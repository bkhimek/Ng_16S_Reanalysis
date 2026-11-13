#!/usr/bin/env bash
# 14_tree.sh - de novo phylogenetic tree of the ASVs, per DADA2 arm.
#
# Needed for phylogeny-aware diversity: Faith's PD (richness weighted by
# branch length) and UniFrac distances. Pipeline (q2-phylogeny 2026.10.0):
#   MAFFT alignment -> mask uninformative/gappy columns -> FastTree (unrooted)
#   -> midpoint rooting. Parameters written out (2026.10.0 defaults):
#   mask max gap frequency 1.0, mask min conservation 0.4, no parttree.
# Check afterwards: the number of tips must equal the number of ASVs.
# Usage (repository root, QIIME 2 environment active):
#   bash scripts/14_tree.sh            (both arms)
#   DRY_RUN=1 bash scripts/14_tree.sh
set -euo pipefail

ARMS="${ARMS:-trunc200_150 trunc218_166}"; THREADS="${THREADS:-6}"
run_cmd() { echo "+ [$(date +%H:%M:%S)] $*"; [[ "${DRY_RUN:-0}" == "1" ]] || "$@"; }
if [[ "${DRY_RUN:-0}" != "1" ]]; then
  command -v qiime >/dev/null || { echo "qiime not found: conda activate rachis-qiime2-2026.10"; exit 1; }
fi
mkdir -p data/qiime/phylogeny
exec 9>"data/qiime/phylogeny/.tree.lock"
flock -n 9 || { echo "Another tree run is active; exiting."; exit 1; }

for arm in $ARMS; do
  echo "=== $arm ==="
  Q="data/qiime/phylogeny/$arm"; R="results/phylogeny/$arm"; mkdir -p "$Q" "$R"
  if [[ -f "$Q/rooted_tree.qza" ]]; then echo "skip: tree exists"; continue; fi
  run_cmd qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "data/qiime/dada2/$arm/humanised_repseqs.qza" \
    --p-n-threads "$THREADS" \
    --p-mask-max-gap-frequency 1.0 --p-mask-min-conservation 0.4 \
    --o-alignment "$Q/aligned_repseqs.qza" \
    --o-masked-alignment "$Q/masked_aligned_repseqs.qza" \
    --o-tree "$Q/unrooted_tree.qza" \
    --o-rooted-tree "$Q/rooted_tree.qza"
  run_cmd qiime tools export --input-path "$Q/rooted_tree.qza" --output-path "$R/rooted_tree_export"
  if [[ "${DRY_RUN:-0}" != "1" ]]; then
    n_asv=$(grep -c "^>" "results/dada2/$arm/repseqs_export/dna-sequences.fasta")
    n_tip=$(python3 -c "
import re,sys
t=open('$R/rooted_tree_export/tree.nwk').read()
print(len(re.findall(r'[(,]\s*([^():,;\s]+)\s*:', t)))")
    echo "ASVs: $n_asv; tree tips: $n_tip"
    [[ "$n_asv" == "$n_tip" ]] && echo "OK: every ASV is a tip" || { echo "MISMATCH: tips != ASVs"; exit 1; }
  fi
done
echo "Done."

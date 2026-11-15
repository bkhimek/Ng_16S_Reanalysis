#!/usr/bin/env bash
# 24_ampliseq_to_qiime.sh - Stage 2: the ampliseq ASV table through the Stage 1
# downstream steps with identical settings (spec docs/STAGE2_COMPARISON_SPEC_v1.md,
# comparisons 4-5; D30).
#   1. import results/ampliseq/dada2/ASV_table.tsv and ASV_seqs.fasta (ASV IDs are
#      MD5 of the sequence, as in QIIME 2)
#   2. taxonomy: D24 SILVA 144 classifier, settings of script 12
#   3. tree: settings of script 14; tips must equal ASVs
#   4. q2-boots core-metrics: settings of script 16 (depth 6,400, n 100, seed
#      20261113, mean / medoid, run-only metadata)
# Outputs: data/qiime/ampliseq/ (artifacts), results/ampliseq_q2/ (exports).
# Usage (repository root, QIIME 2 environment active):
#   bash scripts/24_ampliseq_to_qiime.sh
# Then: python3 scripts/25_stage2_compare.py
set -euo pipefail

AMP="results/ampliseq/dada2"
Q="data/qiime/ampliseq"; R="results/ampliseq_q2"
CLF="data/reference/silva144/SILVA_144_SSURef_NR99_uniform_classifier_V4-515f-806r.qza"
JOBS="${JOBS:-4}"; THREADS="${THREADS:-6}"
DEPTH=6400; N=100; SEED=20261113
run_cmd() { echo "+ [$(date +%H:%M:%S)] $*"; [[ "${DRY_RUN:-0}" == "1" ]] || "$@"; }
if [[ "${DRY_RUN:-0}" != "1" ]]; then
  command -v qiime >/dev/null || { echo "qiime not found: conda activate rachis-qiime2-2026.10"; exit 1; }
  for f in "$AMP/ASV_table.tsv" "$AMP/ASV_seqs.fasta" "$CLF" metadata/run_only_metadata.tsv; do
    [[ -f "$f" ]] || { echo "missing: $f"; exit 1; }
  done
fi
mkdir -p "$Q" "$R"
exec 9>"$Q/.ampliseq_q2.lock"
flock -n 9 || { echo "Another run of this script is active; exiting."; exit 1; }

echo "=== 1. import ==="
if [[ ! -f "$Q/table.qza" ]]; then
  # ASV_table.tsv: ASV_ID <tab> samples...  ->  BIOM-style TSV header
  # not through run_cmd: its echo would be redirected into the table file
  echo "+ [$(date +%H:%M:%S)] awk: $AMP/ASV_table.tsv -> $R/feature-table.tsv"
  awk 'BEGIN{FS=OFS="\t"} NR==1 {print "# Constructed from ampliseq ASV_table.tsv"; $1="#OTU ID"} {print}' \
    "$AMP/ASV_table.tsv" > "$R/feature-table.tsv"
  run_cmd biom convert -i "$R/feature-table.tsv" -o "$R/feature-table.biom" --to-hdf5 --table-type="OTU table"
  run_cmd qiime tools import --type 'FeatureTable[Frequency]' --input-format BIOMV210Format \
    --input-path "$R/feature-table.biom" --output-path "$Q/table.qza"
  run_cmd qiime tools import --type 'FeatureData[Sequence]' \
    --input-path "$AMP/ASV_seqs.fasta" --output-path "$Q/repseqs.qza"
else echo "skip: imported"; fi

echo "=== 2. taxonomy (D24 classifier, script 12 settings) ==="
if [[ ! -f "$Q/taxonomy.qza" ]]; then
  run_cmd qiime feature-classifier classify-sklearn \
    --i-reads "$Q/repseqs.qza" --i-classifier "$CLF" \
    --p-confidence 0.7 --p-read-orientation auto \
    --p-reads-per-batch auto --p-n-jobs "$JOBS" \
    --o-classification "$Q/taxonomy.qza"
  run_cmd qiime tools export --input-path "$Q/taxonomy.qza" --output-path "$R/taxonomy_export"
else echo "skip: classified"; fi

echo "=== 3. tree (script 14 settings) ==="
if [[ ! -f "$Q/rooted_tree.qza" ]]; then
  run_cmd qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "$Q/repseqs.qza" --p-n-threads "$THREADS" \
    --p-mask-max-gap-frequency 1.0 --p-mask-min-conservation 0.4 \
    --o-alignment "$Q/aligned_repseqs.qza" --o-masked-alignment "$Q/masked_aligned_repseqs.qza" \
    --o-tree "$Q/unrooted_tree.qza" --o-rooted-tree "$Q/rooted_tree.qza"
  run_cmd qiime tools export --input-path "$Q/rooted_tree.qza" --output-path "$R/tree_export"
  if [[ "${DRY_RUN:-0}" != "1" ]]; then
    n_asv=$(grep -c "^>" "$AMP/ASV_seqs.fasta")
    n_tip=$(python3 -c "
import re
t=open('$R/tree_export/tree.nwk').read()
print(len(re.findall(r'[(,]\s*([^():,;\s]+)\s*:', t)))")
    echo "ASVs: $n_asv; tree tips: $n_tip"
    [[ "$n_asv" == "$n_tip" ]] || { echo "MISMATCH: tips != ASVs"; exit 1; }
  fi
else echo "skip: tree exists"; fi

echo "=== 4. q2-boots core metrics (script 16 settings) ==="
if [[ ! -d "$Q/diversity/alpha_diversities" ]]; then
  mkdir -p "$Q/diversity" "$R/diversity"
  run_cmd qiime boots core-metrics \
    --i-table "$Q/table.qza" --i-phylogeny "$Q/rooted_tree.qza" \
    --m-metadata-file metadata/run_only_metadata.tsv \
    --p-n "$N" --p-sampling-depth "$DEPTH" --p-no-replacement \
    --p-alpha-average-method mean --p-beta-average-method medoid \
    --p-random-seed "$SEED" \
    --o-resampled-tables "$Q/diversity/resampled_tables" \
    --o-alpha-diversities "$Q/diversity/alpha_diversities" \
    --o-distance-matrices "$Q/diversity/distance_matrices" \
    --o-pcoas "$Q/diversity/pcoas" \
    --o-emperor-plots "$R/diversity/emperor_by_run" \
    --o-scatter-plot "$R/diversity/core_scatter.qzv"
  [[ "${DRY_RUN:-0}" == "1" ]] || { ls "$Q/diversity/alpha_diversities" "$Q/diversity/distance_matrices"; }
else echo "skip: core metrics exist"; fi
echo "Done. Next: python3 scripts/25_stage2_compare.py"

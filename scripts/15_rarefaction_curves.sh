#!/usr/bin/env bash
# 15_rarefaction_curves.sh - technical check of the rarefaction depth (D27).
#
# Alpha-rarefaction curves (observed ASVs, Faith's PD, Shannon) from 400 to
# 19,200 reads in steps of 400 (48 depths; 6,400 and 12,800 are on the grid),
# 10 random draws per depth. The metadata given to QIIME 2 contains ONLY the
# sequencing run, so the curves cannot be grouped by diet or time (D27: no
# diversity result is viewed before the analysis rules are fixed).
# Then: python3 scripts/15_rarefaction_check.py
# Usage (repository root, QIIME 2 environment active):
#   bash scripts/15_rarefaction_curves.sh          (both arms)
set -euo pipefail

ARMS="${ARMS:-trunc200_150 trunc218_166}"
run_cmd() { echo "+ [$(date +%H:%M:%S)] $*"; [[ "${DRY_RUN:-0}" == "1" ]] || "$@"; }
if [[ "${DRY_RUN:-0}" != "1" ]]; then
  command -v qiime >/dev/null || { echo "qiime not found: conda activate rachis-qiime2-2026.10"; exit 1; }
fi
mkdir -p results/diversity
exec 9>"results/diversity/.rarefaction.lock"
flock -n 9 || { echo "Another rarefaction run is active; exiting."; exit 1; }

# run-only metadata (sample-id, run), derived from the full metadata
awk -F'\t' 'NR==1 { for (i=1;i<=NF;i++) if ($i=="run") c=i; print "sample-id\trun"; next }
            $1 ~ /^#q2/ { next }
            { print $1 "\t" $c }' metadata/humanised_sample_metadata.tsv > metadata/run_only_metadata.tsv
echo "run-only metadata: $(($(wc -l < metadata/run_only_metadata.tsv) - 1)) samples, columns: $(head -1 metadata/run_only_metadata.tsv | tr '\t' ' ')"

for arm in $ARMS; do
  echo "=== $arm ==="
  R="results/diversity/$arm"; mkdir -p "$R"
  run_cmd qiime diversity alpha-rarefaction \
    --i-table "data/qiime/dada2/$arm/humanised_table.qza" \
    --i-phylogeny "data/qiime/phylogeny/$arm/rooted_tree.qza" \
    --p-metrics observed_features --p-metrics faith_pd --p-metrics shannon \
    --m-metadata-file metadata/run_only_metadata.tsv \
    --p-min-depth 400 --p-max-depth 19200 --p-steps 48 --p-iterations 10 \
    --o-visualization "$R/alpha_rarefaction_by_run.qzv"
done
echo "Done. Next: python3 scripts/15_rarefaction_check.py $ARMS"

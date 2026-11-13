#!/usr/bin/env bash
# 16_core_metrics.sh - alpha and beta diversity by repeated rarefaction (D27).
#
# q2-boots 2026.10.0 core-metrics: 100 rarefactions without replacement to
# 6,400 reads per sample, fixed random seed (results reproducible exactly);
# alpha diversity averaged by mean, distance matrices by medoid.
# Metadata given = run only (metadata/run_only_metadata.tsv, made by script 15),
# so the plots this action makes cannot show diet or time: biological results
# are examined only after the analysis specification is frozen.
# Outputs (gitignored): data/qiime/diversity/<arm>/ (collections of alpha
# vectors, distance matrices, PCoAs) and results/diversity/<arm>/core_scatter.qzv
# Usage (repository root, QIIME 2 environment active):
#   bash scripts/16_core_metrics.sh            (both arms)
# Check progress with: bash scripts/status.sh
set -euo pipefail

ARMS="${ARMS:-trunc200_150 trunc218_166}"
DEPTH=6400; N=100; SEED=20261113
run_cmd() { echo "+ [$(date +%H:%M:%S)] $*"; [[ "${DRY_RUN:-0}" == "1" ]] || "$@"; }
if [[ "${DRY_RUN:-0}" != "1" ]]; then
  command -v qiime >/dev/null || { echo "qiime not found: conda activate rachis-qiime2-2026.10"; exit 1; }
fi
[[ -f metadata/run_only_metadata.tsv ]] || { echo "metadata/run_only_metadata.tsv missing: run script 15 first"; exit 1; }
mkdir -p data/qiime/diversity
exec 9>"data/qiime/diversity/.core.lock"
flock -n 9 || { echo "Another core-metrics run is active; exiting."; exit 1; }

for arm in $ARMS; do
  echo "=== $arm ==="
  Q="data/qiime/diversity/$arm"; R="results/diversity/$arm"; mkdir -p "$R"
  if [[ -d "$Q/alpha_diversities" ]]; then echo "skip: core metrics exist"; continue; fi
  mkdir -p "$Q"
  run_cmd qiime boots core-metrics \
    --i-table "data/qiime/dada2/$arm/humanised_table.qza" \
    --i-phylogeny "data/qiime/phylogeny/$arm/rooted_tree.qza" \
    --m-metadata-file metadata/run_only_metadata.tsv \
    --p-n "$N" --p-sampling-depth "$DEPTH" --p-no-replacement \
    --p-alpha-average-method mean --p-beta-average-method medoid \
    --p-random-seed "$SEED" \
    --o-resampled-tables "$Q/resampled_tables" \
    --o-alpha-diversities "$Q/alpha_diversities" \
    --o-distance-matrices "$Q/distance_matrices" \
    --o-pcoas "$Q/pcoas" \
    --o-emperor-plots "$R/emperor_by_run" \
    --o-scatter-plot "$R/core_scatter.qzv"
  [[ "${DRY_RUN:-0}" == "1" ]] || { echo "alpha metrics:"; ls "$Q/alpha_diversities"; echo "distance matrices:"; ls "$Q/distance_matrices"; }
done
echo "Done."

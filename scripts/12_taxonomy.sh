#!/usr/bin/env bash
# 12_taxonomy.sh - assign taxonomy with the SILVA 144 V4 classifier (D24).
#
#   TRIAL=1 bash scripts/12_taxonomy.sh   -> classify the first 20 ASVs of the
#                                            primary arm only (minutes)
#   bash scripts/12_taxonomy.sh           -> classify all ASVs of both arms
# The same classifier and settings are used for every arm (D24). Settings are
# written out explicitly: confidence 0.7, read orientation auto (q2-feature-
# classifier 2026.10.0 defaults), reads per batch auto.
# Outputs: data/qiime/taxonomy/<arm>/taxonomy.qza; results/taxonomy/<arm>/
# (taxonomy.qzv and exported taxonomy.tsv). Check with: python3 scripts/13_taxonomy_check.py
# Check progress with: bash scripts/status.sh
set -euo pipefail

CLF="data/reference/silva144/SILVA_144_SSURef_NR99_uniform_classifier_V4-515f-806r.qza"
ARMS="${ARMS:-trunc200_150 trunc218_166}"
JOBS="${JOBS:-4}"
run_cmd() { echo "+ [$(date +%H:%M:%S)] $*"; [[ "${DRY_RUN:-0}" == "1" ]] || "$@"; }
if [[ "${DRY_RUN:-0}" != "1" ]]; then
  command -v qiime >/dev/null || { echo "qiime not found: conda activate rachis-qiime2-2026.10"; exit 1; }
  [[ -f "$CLF" ]] || { echo "classifier missing: run scripts/11_get_classifier.sh"; exit 1; }
fi
mkdir -p data/qiime/taxonomy
exec 9>"data/qiime/taxonomy/.taxonomy.lock"
flock -n 9 || { echo "Another taxonomy run is active; exiting."; exit 1; }

classify() {  # $1 reads.qza, $2 output dir (artifacts), $3 results dir
  mkdir -p "$2" "$3"
  run_cmd qiime feature-classifier classify-sklearn \
    --i-reads "$1" --i-classifier "$CLF" \
    --p-confidence 0.7 --p-read-orientation auto \
    --p-reads-per-batch auto --p-n-jobs "$JOBS" \
    --o-classification "$2/taxonomy.qza"
  run_cmd qiime metadata tabulate --m-input-file "$2/taxonomy.qza" --o-visualization "$3/taxonomy.qzv"
  run_cmd qiime tools export --input-path "$2/taxonomy.qza" --output-path "$3/export"
}

if [[ "${TRIAL:-0}" == "1" ]]; then
  arm=trunc200_150
  T="data/qiime/taxonomy/trial"; mkdir -p "$T" results/taxonomy/trial
  { echo "feature-id"; grep "^>" "results/dada2/$arm/repseqs_export/dna-sequences.fasta" | head -20 | sed 's/^>//; s/ .*//'; } > "$T/trial_ids.tsv"
  run_cmd qiime feature-table filter-seqs \
    --i-data "data/qiime/dada2/$arm/humanised_repseqs.qza" \
    --m-metadata-file "$T/trial_ids.tsv" \
    --o-filtered-data "$T/trial_repseqs.qza"
  classify "$T/trial_repseqs.qza" "$T" results/taxonomy/trial
  echo "Trial taxonomy (feature, taxon, confidence):"
  [[ "${DRY_RUN:-0}" == "1" ]] || cut -f2,3 results/taxonomy/trial/export/taxonomy.tsv | head -21
  echo "Done (trial)."
  exit 0
fi

for arm in $ARMS; do
  echo "=== $arm ==="
  if [[ -f "data/qiime/taxonomy/$arm/taxonomy.qza" ]]; then echo "skip: classified already"; continue; fi
  classify "data/qiime/dada2/$arm/humanised_repseqs.qza" "data/qiime/taxonomy/$arm" "results/taxonomy/$arm"
done
echo "Done. Next: python3 scripts/13_taxonomy_check.py $ARMS"

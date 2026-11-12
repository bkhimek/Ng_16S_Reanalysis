#!/usr/bin/env bash
# 08_dada2_humanised.sh - DADA2 denoising of the 3 humanised runs (QIIME 2).
#
# Decisions applied: D19 (two failed libraries excluded before denoising),
# D20 (truncation 218/166 primary; 200/150 sensitivity; same lengths in all runs).
# Each run is denoised SEPARATELY (its own error model: separate MiSeq runs, D14),
# then the per-run tables and representative sequences are merged. Merging works
# only because every run uses identical truncation: an ASV is an exact sequence.
#
# All DADA2 parameters are written out explicitly (most are q2-dada2 2026.10.1
# defaults) so the record does not depend on future default changes.
#
# Settings by environment variable (defaults = primary arm):
#   TRUNC_F=218 TRUNC_R=166 TAG=trunc218_166 THREADS=6
# Sensitivity arm:
#   TRUNC_F=200 TRUNC_R=150 TAG=trunc200_150 bash scripts/08_dada2_humanised.sh
# DRY_RUN=1 prints the commands without running them.
# Outputs: data/qiime/dada2/<TAG>/ (artifacts), results/dada2/<TAG>/ (qzv, TSV).
set -euo pipefail

TRUNC_F="${TRUNC_F:-218}"; TRUNC_R="${TRUNC_R:-166}"
TAG="${TAG:-trunc${TRUNC_F}_${TRUNC_R}}"; THREADS="${THREADS:-6}"
DEMUX="data/qiime/demux"
Q="data/qiime/dada2/$TAG"; R="results/dada2/$TAG"
EXCL="metadata/excluded_samples.tsv"
mkdir -p "$Q" "$R"

exec 9>"$Q/.dada2.lock"
flock -n 9 || { echo "Another DADA2 run is using $Q; exiting."; exit 1; }

run_cmd() { echo "+ $*"; [[ "${DRY_RUN:-0}" == "1" ]] || "$@"; }
if [[ "${DRY_RUN:-0}" != "1" ]]; then
  command -v qiime >/dev/null || { echo "qiime not found: conda activate rachis-qiime2-2026.10"; exit 1; }
fi
echo "DADA2 arm $TAG: trunc-len-f $TRUNC_F, trunc-len-r $TRUNC_R, threads $THREADS"
printf "arm\t%s\ntrunc_len_f\t%s\ntrunc_len_r\t%s\nstarted\t%s\n" "$TAG" "$TRUNC_F" "$TRUNC_R" "$(date -Is)" > "$R/arm_settings.tsv"

for run in hum_run1 hum_run2 hum_run3; do
  echo "=== $run ==="
  if [[ -f "$Q/${run}_table.qza" ]]; then echo "skip: denoised already"; continue; fi

  input="$DEMUX/${run}_demux.qza"
  # D19: exclude failed libraries, only in runs that contain them
  n_excl=$(comm -12 <(tail -n +2 "$EXCL" | cut -f1 | sort) \
                    <(tail -n +2 "metadata/demux/${run}_barcodes.tsv" | cut -f1 | sort) | wc -l)
  if [[ "$n_excl" -gt 0 ]]; then
    echo "excluding $n_excl failed sample(s) listed in $EXCL"
    run_cmd qiime demux filter-samples \
      --i-demux "$input" \
      --m-metadata-file "$EXCL" \
      --p-exclude-ids \
      --o-filtered-demux "$Q/${run}_demux_filtered.qza"
    input="$Q/${run}_demux_filtered.qza"
  fi

  run_cmd qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "$input" \
    --p-trunc-len-f "$TRUNC_F" --p-trunc-len-r "$TRUNC_R" \
    --p-trim-left-f 0 --p-trim-left-r 0 \
    --p-max-ee-f 2 --p-max-ee-r 2 --p-trunc-q 2 \
    --p-min-overlap 12 --p-max-merge-mismatch 0 \
    --p-pooling-method independent --p-chimera-method consensus \
    --p-min-fold-parent-over-abundance 1 \
    --p-n-reads-learn 1000000 --p-n-threads "$THREADS" \
    --o-table "$Q/${run}_table.qza" \
    --o-representative-sequences "$Q/${run}_repseqs.qza" \
    --o-denoising-stats "$Q/${run}_stats.qza" \
    --o-base-transition-stats "$Q/${run}_base_transition_stats.qza"

  run_cmd qiime tools export --input-path "$Q/${run}_stats.qza" --output-path "$R/${run}_stats"
  run_cmd rm -f "$Q/${run}_demux_filtered.qza"
done

echo "=== merge runs ==="
run_cmd qiime feature-table merge \
  --i-tables "$Q/hum_run1_table.qza" "$Q/hum_run2_table.qza" "$Q/hum_run3_table.qza" \
  --o-merged-table "$Q/humanised_table.qza"
run_cmd qiime feature-table merge-seqs \
  --i-data "$Q/hum_run1_repseqs.qza" "$Q/hum_run2_repseqs.qza" "$Q/hum_run3_repseqs.qza" \
  --o-merged-data "$Q/humanised_repseqs.qza"
run_cmd qiime feature-table summarize \
  --i-table "$Q/humanised_table.qza" \
  --m-sample-metadata-file metadata/humanised_sample_metadata.tsv \
  --o-visualization "$R/humanised_table.qzv"

# plain-text copies for the retention check (09) and later comparisons
run_cmd qiime tools export --input-path "$Q/humanised_table.qza" --output-path "$R/table_export"
run_cmd biom convert -i "$R/table_export/feature-table.biom" -o "$R/table_export/feature-table.tsv" --to-tsv
run_cmd qiime tools export --input-path "$Q/humanised_repseqs.qza" --output-path "$R/repseqs_export"
printf "finished\t%s\n" "$(date -Is)" >> "$R/arm_settings.tsv"
echo "Done. Next: python3 scripts/09_dada2_check.py $TAG"

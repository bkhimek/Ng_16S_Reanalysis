#!/usr/bin/env bash
# 06_demux_humanised.sh - import and demultiplex the 3 humanised runs (QIIME 2).
#
# Per run:
#   1. hard-link the raw files under the names QIIME 2 expects for an
#      EMP paired-end import (forward.fastq.gz, barcodes.fastq.gz,
#      reverse.fastq.gz); hard links take no extra disk space;
#   2. qiime tools import  --type EMPPairedEndSequences;
#   3. qiime demux emp-paired with Golay error correction ON and BOTH
#      reverse-complement options ON (D17): the index reads and the workbook
#      barcodes are both the reverse complement of valid Golay codewords;
#   4. qiime demux summarize (read counts and quality plots, .qzv);
#   5. export the error-correction details to TSV for the assignment check.
# Runs are processed separately because barcodes are reused between runs.
# A run whose outputs already exist is skipped (re-runnable after a failure).
# A lock file stops a second copy from running at the same time.
# DRY_RUN=1 prints the commands without running them.
# RAW, QDIR, RESDIR override the input/output folders (used for the smoke test).
#
# Usage (repository root, QIIME 2 environment active):
#   conda activate rachis-qiime2-2026.10
#   bash scripts/06_demux_humanised.sh
set -euo pipefail

RAW="${RAW:-data/raw/humanised}"     # override for the smoke test
Q="${QDIR:-data/qiime/demux}"
R="${RESDIR:-results/demux}"
mkdir -p "$Q" "$R"

# refuse to start if another copy is already running (same output folder)
exec 9>"$Q/.demux.lock"
flock -n 9 || { echo "Another demultiplexing run is using $Q; exiting."; exit 1; }

run_cmd() { echo "+ $*"; [[ "${DRY_RUN:-0}" == "1" ]] || "$@"; }

if [[ "${DRY_RUN:-0}" != "1" ]]; then
  command -v qiime >/dev/null || { echo "qiime not found: conda activate rachis-qiime2-2026.10"; exit 1; }
fi

declare -A FWD=( [hum_run1]=Figure_3_humanized_double_ciprofloxacin_prep1_R1.fastq.gz
                 [hum_run2]=Figure_3_humanized_double_ciprofloxacin_prep2_R1.fastq.gz
                 [hum_run3]=Figure_3_humanized_double_ciprofloxacin_prep3_R1_001.fastq.gz )
declare -A IDX=( [hum_run1]=Figure_3_humanized_double_ciprofloxacin_prep1_R2.fastq.gz
                 [hum_run2]=Figure_3_humanized_double_ciprofloxacin_prep2_R2.fastq.gz
                 [hum_run3]=Figure_3_humanized_double_ciprofloxacin_prep3_I1_001.fastq.gz )
declare -A REV=( [hum_run1]=Figure_3_humanized_double_ciprofloxacin_prep1_R3.fastq.gz
                 [hum_run2]=Figure_3_humanized_double_ciprofloxacin_prep2_R3.fastq.gz
                 [hum_run3]=Figure_3_humanized_double_ciprofloxacin_prep3_R2_001.fastq.gz )

for run in hum_run1 hum_run2 hum_run3; do
  echo "=== $run ==="
  if [[ -f "$Q/${run}_demux.qza" && -f "$R/${run}_demux.qzv" ]]; then
    echo "skip: outputs exist"; continue
  fi
  imp="$Q/${run}_import"
  run_cmd mkdir -p "$imp"
  run_cmd ln -f "$RAW/${FWD[$run]}" "$imp/forward.fastq.gz"
  run_cmd ln -f "$RAW/${IDX[$run]}" "$imp/barcodes.fastq.gz"
  run_cmd ln -f "$RAW/${REV[$run]}" "$imp/reverse.fastq.gz"

  run_cmd qiime tools import \
    --type EMPPairedEndSequences \
    --input-path "$imp" \
    --output-path "$Q/${run}_multiplexed.qza"

  run_cmd qiime demux emp-paired \
    --i-seqs "$Q/${run}_multiplexed.qza" \
    --m-barcodes-file "metadata/demux/${run}_barcodes.tsv" \
    --m-barcodes-column barcode-sequence \
    --p-golay-error-correction \
    --p-rev-comp-barcodes \
    --p-rev-comp-mapping-barcodes \
    --o-per-sample-sequences "$Q/${run}_demux.qza" \
    --o-error-correction-details "$Q/${run}_ecc_details.qza"

  run_cmd qiime demux summarize \
    --i-data "$Q/${run}_demux.qza" \
    --o-visualization "$R/${run}_demux.qzv"

  run_cmd qiime tools export \
    --input-path "$Q/${run}_ecc_details.qza" \
    --output-path "$R/${run}_ecc_details"

  # the multiplexed artifact is a full copy of the raw reads; remove it once
  # demultiplexing has succeeded (the raw files remain the source of truth)
  run_cmd rm -f "$Q/${run}_multiplexed.qza"
done
echo "Done. Next: python3 scripts/07_demux_check.py"

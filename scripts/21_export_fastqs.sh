#!/usr/bin/env bash
# 21_export_fastqs.sh - per-sample FASTQs for Stage 2 (D9).
#
# Exports the Stage 1 demultiplexed artifacts (data/qiime/demux/<run>_demux.qza)
# to per-sample gzipped FASTQ files, so that QIIME 2 (Stage 1) and
# nf-core/ampliseq (Stage 2) start from identical reads. Then
# scripts/21_fastq_manifest.py checks and records them:
#   - the 451 samples Stage 1 denoised (D19 failures excluded);
#   - forward and reverse read counts equal, and equal to the Stage 1 DADA2
#     'input' count of the same sample;
#   - SHA-256 of every file -> metadata/per_sample_fastq_sha256.tsv (committed);
#   - ampliseq sample sheet -> config/ampliseq/samplesheet.tsv (sample,
#     fastq_1, fastq_2, run; absolute paths on this machine).
# Outputs data/fastq_per_sample/<run>/ (gitignored).
# Usage (repository root, QIIME 2 environment active):
#   bash scripts/21_export_fastqs.sh
set -euo pipefail

OUTD="${OUTD:-data/fastq_per_sample}"
run_cmd() { echo "+ [$(date +%H:%M:%S)] $*"; [[ "${DRY_RUN:-0}" == "1" ]] || "$@"; }
if [[ "${DRY_RUN:-0}" != "1" ]]; then
  command -v qiime >/dev/null || { echo "qiime not found: conda activate rachis-qiime2-2026.10"; exit 1; }
fi
mkdir -p "$OUTD"
exec 9>"$OUTD/.export.lock"
flock -n 9 || { echo "Another export is active; exiting."; exit 1; }

for run in hum_run1 hum_run2 hum_run3; do
  echo "=== $run ==="
  d="$OUTD/$run"
  if [[ -f "$d/.done" ]]; then echo "skip: exported already"; continue; fi
  run_cmd rm -rf "$d.tmp" "$d"
  run_cmd qiime tools export --input-path "data/qiime/demux/${run}_demux.qza" --output-path "$d.tmp"
  run_cmd mv "$d.tmp" "$d"
  run_cmd touch "$d/.done"
done

run_cmd python3 scripts/21_fastq_manifest.py "$OUTD"
echo "Done."

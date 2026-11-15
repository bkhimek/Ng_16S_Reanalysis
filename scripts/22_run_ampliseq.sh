#!/usr/bin/env bash
# 22_run_ampliseq.sh - nf-core/ampliseq 3.1.0, matched Stage 2 run (D30).
#
#   bash scripts/22_run_ampliseq.sh smoke   6 samples of hum_run2 -> results/ampliseq_smoke
#   bash scripts/22_run_ampliseq.sh full    all 451 samples       -> results/ampliseq
# Settings: config/ampliseq/params_matched.yaml (parameters) and
# config/ampliseq/matched.config (module arguments, resources);
# see docs/STAGE2_PARAMETER_MATCH.md. Nextflow version pinned below.
# After the run: Nextflow/pipeline versions and container image digests are
# written to <outdir>/run_versions.txt. Rerunning resumes (-resume).
# Smoke check afterwards: python3 scripts/23_ampliseq_smoke_check.py
set -euo pipefail

MODE="${1:-}"
export NXF_VER="26.10.5"      # stable, 2026-11-02; ampliseq 3.1.0 supports 26.04+ and 26.10
REV="3.1.0"
SHEET_ALL="config/ampliseq/samplesheet.tsv"
case "$MODE" in
  smoke) SHEET="config/ampliseq/samplesheet_smoke.tsv"; OUT="results/ampliseq_smoke"; WORK="work/ampliseq_smoke"
         # one awk, no pipe: a pipe into head under pipefail exits with SIGPIPE (code 141)
         awk -F'\t' 'NR==1 || ($4=="hum_run2" && ++n<=6)' "$SHEET_ALL" > "$SHEET"
         echo "smoke sample sheet: $(($(wc -l < "$SHEET") - 1)) samples (hum_run2)" ;;
  full)  SHEET="$SHEET_ALL"; OUT="results/ampliseq"; WORK="work/ampliseq" ;;
  *)     echo "Usage: bash scripts/22_run_ampliseq.sh smoke|full"; exit 1 ;;
esac
command -v nextflow >/dev/null || { echo "nextflow not found (see HANDOVER: install to ~/bin)"; exit 1; }
docker info >/dev/null 2>&1 || { echo "docker daemon not reachable"; exit 1; }
mkdir -p "$OUT" "$WORK"
exec 9>"$OUT/.ampliseq.lock"
flock -n 9 || { echo "Another ampliseq run with this output is active; exiting."; exit 1; }

echo "+ [$(date +%H:%M:%S)] nextflow run nf-core/ampliseq -r $REV ($MODE)"
nextflow run nf-core/ampliseq -r "$REV" -profile docker \
  -params-file config/ampliseq/params_matched.yaml \
  -c config/ampliseq/matched.config \
  --input "$SHEET" --outdir "$OUT" \
  -work-dir "$WORK" -resume

nf_ver=$(nextflow -version 2>&1); nf_info=$(nextflow info nf-core/ampliseq 2>&1 || true); java_ver=$(java -version 2>&1)
imgs=$(docker images --digests --format '  {{.Repository}}:{{.Tag}} {{.Digest}}')
{
  echo "date: $(date -Iseconds)"
  echo "mode: $MODE"
  echo "nextflow: $(awk '/version/ {print $2, $3, $4; exit}' <<< "$nf_ver")"
  echo "pipeline: nf-core/ampliseq $REV"
  echo "$nf_info" | sed 's/^/  /'
  echo "java: $(awk '/version/ {print; exit}' <<< "$java_ver")"
  echo "container images (repository, tag, digest):"
  awk '!/hello-world/' <<< "$imgs" | sort -u
} > "$OUT/run_versions.txt"
echo "Wrote $OUT/run_versions.txt"
echo "Done ($MODE)."

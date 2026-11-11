#!/usr/bin/env bash
# 04_download_raw.sh <humanised|conventional>
# Download the full Figure 3 FASTQ files of one experiment and verify them.
#
# The repository publishes no checksums, so each file must pass three checks:
#   1. exact byte size equal to the inventory (metadata/sdr_inventory_figure3.tsv,
#      taken from the server's Content-Length on 2026-11-11);
#   2. gzip -t (the whole compressed stream decodes without error);
#   3. per run: forward, index and reverse files hold the same number of reads.
# A file is downloaded to <name>.part and renamed only after checks 1-2 pass,
# so a file with its final name is always complete. No resume (curl -C -) is
# used: resumed transfers have silently corrupted files in an earlier project.
# Re-running the script skips files already verified (idempotent), and a lock
# prevents two copies running at once.
# A sha256 is computed for every verified file and written to
# metadata/raw_sha256_<experiment>.txt, so later copies can be checked exactly.
#
# Usage (repository root):  bash scripts/04_download_raw.sh humanised
set -euo pipefail

EXP="${1:?usage: $0 <humanised|conventional>}"
BASE="${SDR_BASE:-https://stacks.stanford.edu/file/rm696bw4911}"
INV="metadata/sdr_inventory_figure3.tsv"
OUT="data/raw/${EXP}"
LOG="data/raw/download_log_${EXP}.tsv"
SHA="metadata/raw_sha256_${EXP}.txt"
mkdir -p "$OUT"

exec 9>"data/raw/.download_${EXP}.lock"
flock -n 9 || { echo "Another download for ${EXP} is running; exiting."; exit 1; }

[[ -f "$LOG" ]] || printf "time\tfile\texpected_bytes\tbytes\tgzip_test\tstatus\n" > "$LOG"
ZCAT=$(command -v pigz >/dev/null && echo "pigz -dc" || echo "gzip -dc")

verify() {  # $1 file, $2 expected bytes -> 0 if size and gzip ok
  local f="$1" exp="$2" got
  got=$(stat -c %s "$f")
  [[ "$got" == "$exp" ]] || { echo "size $got != $exp"; return 1; }
  gzip -t "$f" 2>/dev/null || { echo "gzip -t failed"; return 1; }
  return 0
}

mapfile -t ROWS < <(awk -F'\t' -v e="$EXP" 'NR>1 && $1==e {print $4"\t"$5"\t"$2"\t"$3}' "$INV")
[[ ${#ROWS[@]} -eq 9 ]] || { echo "Expected 9 inventory rows for ${EXP}, found ${#ROWS[@]}"; exit 1; }

for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r name bytes run role <<< "$row"
  dest="$OUT/$name"
  if [[ -f "$dest" && -f "$dest.verified" ]]; then
    echo "skip (verified earlier): $name"; continue
  fi
  for attempt in 1 2; do
    echo "downloading ($attempt/2): $name  [$bytes bytes]"
    rm -f "$dest.part"
    code=$(curl -sS --fail -o "$dest.part" -w "%{http_code}" "$BASE/$name" || echo "curl_error")
    if [[ "$code" == "200" ]] && msg=$(verify "$dest.part" "$bytes"); then
      mv "$dest.part" "$dest"; touch "$dest.verified"
      printf "%s\t%s\t%s\t%s\tok\tOK\n" "$(date -Is)" "$name" "$bytes" "$(stat -c %s "$dest")" >> "$LOG"
      break
    fi
    printf "%s\t%s\t%s\t%s\t%s\tFAIL(http %s)\n" "$(date -Is)" "$name" "$bytes" \
      "$(stat -c %s "$dest.part" 2>/dev/null || echo 0)" "${msg:-n/a}" "$code" >> "$LOG"
    [[ $attempt -eq 2 ]] && { echo "FAILED twice: $name (see $LOG)"; exit 1; }
  done
done

echo; echo "Read counts per run (forward / index / reverse must be equal):"
fail=0
for run in $(printf "%s\n" "${ROWS[@]}" | cut -f3 | sort -u); do
  declare -A n=()
  for row in "${ROWS[@]}"; do
    IFS=$'\t' read -r name bytes r role <<< "$row"
    [[ "$r" == "$run" ]] || continue
    n[$role]=$(( $($ZCAT "$OUT/$name" | wc -l) / 4 ))
  done
  status=$([[ "${n[forward]}" == "${n[index]}" && "${n[index]}" == "${n[reverse]}" ]] && echo OK || echo MISMATCH)
  [[ "$status" == OK ]] || fail=1
  printf "  %s %s: forward %s, index %s, reverse %s  %s\n" "$EXP" "$run" "${n[forward]}" "${n[index]}" "${n[reverse]}" "$status" | tee -a "$LOG"
  unset n
done

( cd "$OUT" && sha256sum *.fastq.gz ) > "$SHA"
echo; echo "sha256 written to $SHA"
[[ $fail -eq 0 ]] && echo "All ${EXP} files verified." || { echo "Read-count MISMATCH: do not use these files."; exit 1; }

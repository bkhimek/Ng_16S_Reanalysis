#!/usr/bin/env bash
# 00_fetch_read_heads.sh
# Download only the first 4 MiB of each of the 18 Figure 3 FASTQ files
# (forward, index and reverse reads for 6 MiSeq runs). A gzip stream can be
# decompressed from its start, so each partial file yields a few thousand
# complete reads: enough to check read lengths, read roles, barcodes and
# primers before committing to the 24 GB download.
#
# Each piece is accepted only if the server answers HTTP 206 (partial content)
# and exactly the requested number of bytes arrived. A 200 means the server
# ignored the range; an HTML error page can also arrive with a success code.
# Usage: bash scripts/00_fetch_read_heads.sh   (run from the repository root)
set -euo pipefail

BASE="https://stacks.stanford.edu/file/rm696bw4911"
OUT="data/readcheck"
BYTES=4194304   # 4 MiB
mkdir -p "$OUT"
LOG="$OUT/fetch_log.tsv"
printf "file\thttp_code\tbytes_received\tfirst_two_bytes\tstatus\n" > "$LOG"

FILES=(
  Figure_3_humanized_double_ciprofloxacin_prep1_R1.fastq.gz
  Figure_3_humanized_double_ciprofloxacin_prep1_R2.fastq.gz
  Figure_3_humanized_double_ciprofloxacin_prep1_R3.fastq.gz
  Figure_3_humanized_double_ciprofloxacin_prep2_R1.fastq.gz
  Figure_3_humanized_double_ciprofloxacin_prep2_R2.fastq.gz
  Figure_3_humanized_double_ciprofloxacin_prep2_R3.fastq.gz
  Figure_3_humanized_double_ciprofloxacin_prep3_R1_001.fastq.gz
  Figure_3_humanized_double_ciprofloxacin_prep3_I1_001.fastq.gz
  Figure_3_humanized_double_ciprofloxacin_prep3_R2_001.fastq.gz
  Figure_3_conventional_ciprofloxacin_prep_1_R1.fastq.gz
  Figure_3_conventional_ciprofloxacin_prep_1_R2.fastq.gz
  Figure_3_conventional_ciprofloxacin_prep_1_R3.fastq.gz
  Figure_3_conventional_ciprofloxacin_prep_2_R1.fastq.gz
  Figure_3_conventional_ciprofloxacin_prep_2_R2.fastq.gz
  Figure_3_conventional_ciprofloxacin_prep_2_R3.fastq.gz
  Figure_3_conventional_ciprofloxacin_prep_3_R1.fastq.gz
  Figure_3_conventional_ciprofloxacin_prep_3_R2.fastq.gz
  Figure_3_conventional_ciprofloxacin_prep_3_R3.fastq.gz
)

fail=0
for f in "${FILES[@]}"; do
  dest="$OUT/${f%.fastq.gz}.head.fastq.gz"
  code=$(curl -sS --retry 2 -r "0-$((BYTES-1))" -o "$dest" -w "%{http_code}" "$BASE/$f" || echo "curl_error")
  got=$(stat -c %s "$dest" 2>/dev/null || echo 0)
  magic=$(head -c 2 "$dest" 2>/dev/null | od -An -tx1 | tr -d ' \n')
  if [[ "$code" == "206" && "$got" == "$BYTES" && "$magic" == "1f8b" ]]; then
    status="OK"
  else
    status="FAIL"; fail=1
  fi
  printf "%s\t%s\t%s\t%s\t%s\n" "$f" "$code" "$got" "$magic" "$status" | tee -a "$LOG"
done

echo
if [[ $fail -eq 0 ]]; then
  echo "All 18 pieces OK (HTTP 206, ${BYTES} bytes, gzip magic 1f8b)."
else
  echo "One or more pieces FAILED; see $LOG. Do not run the read check on failed pieces."
  exit 1
fi

#!/usr/bin/env bash
# 11_get_classifier.sh - download and verify the SILVA 144 V4 classifier (D24).
#
# Source: SILVA official QIIME 2 classifiers,
#   https://www.arb-silva.de/fileadmin/silva_databases/current/QIIME2/2026.7/SSU/V4-515f-806r/uniform/
# Values recorded on 2026-11-13 from SILVA's own files:
#   classifier size 79,326,205 bytes; md5 683a71dcd81ce90d724add22d2ad4d0f;
#   DOI 10.82364/144/SSU/Ref-NR99/QIIME2/2026.7/V4-515f-806r/uniform;
#   README primers GTGYCAGCMGCCGCGGTAA / GGACTACNVGGGTWTCTAAT; trained with QIIME 2 2026.7.0.
#   Signing key: SILVA Team <contact@arb-silva.de>, RSA 4096,
#   fingerprint 01877A6A3933ED067EA5B8C3F3A26EDEDB5A4596 (key text copied from
#   https://www.arb-silva.de/documentation/classifiers/qiime-2 into config/silva/).
# Checks: HTTP success, exact size, md5 against SILVA's .md5 file AND the value
# above, key fingerprint, then (QIIME 2 environment active) the GPG signature
# 'SILVA-Pipeline-Official' with `qiime tools signature-verify`.
# Writes metadata/database_manifest.tsv (one row).
# Usage (repository root):  bash scripts/11_get_classifier.sh
set -euo pipefail

BASE="https://www.arb-silva.de/fileadmin/silva_databases/current/QIIME2/2026.7/SSU/V4-515f-806r/uniform"
QZA="SILVA_144_SSURef_NR99_uniform_classifier_V4-515f-806r.qza"
OUT="data/reference/silva144"
KEY="config/silva/SILVA_PUBLIC_KEY.asc"
EXP_BYTES=79326205
EXP_MD5=683a71dcd81ce90d724add22d2ad4d0f
EXP_FPR=01877A6A3933ED067EA5B8C3F3A26EDEDB5A4596
mkdir -p "$OUT"

fail() { echo "FAILED: $*"; exit 1; }

for f in README "$QZA.md5" "$QZA.doi"; do
  curl -sS --fail -o "$OUT/$f" "$BASE/$f" || fail "download $f"
done
if [[ -f "$OUT/$QZA" && "$(stat -c %s "$OUT/$QZA")" == "$EXP_BYTES" ]]; then
  echo "classifier already present with the expected size; re-checking"
else
  rm -f "$OUT/$QZA.part"
  curl -sS --fail -o "$OUT/$QZA.part" "$BASE/$QZA" || fail "download $QZA"
  mv "$OUT/$QZA.part" "$OUT/$QZA"
fi

got=$(stat -c %s "$OUT/$QZA");      [[ "$got" == "$EXP_BYTES" ]] || fail "size $got != $EXP_BYTES"
md5=$(md5sum "$OUT/$QZA" | cut -d' ' -f1)
[[ "$md5" == "$(cut -d' ' -f1 "$OUT/$QZA.md5")" ]] || fail "md5 does not match SILVA's .md5 file"
[[ "$md5" == "$EXP_MD5" ]] || fail "md5 differs from the value recorded on 2026-11-13 (file changed upstream?)"
echo "OK size $got bytes, md5 $md5"
echo "DOI: $(cat "$OUT/$QZA.doi")"
grep -q "GTGYCAGCMGCCGCGGTAA-GGACTACNVGGGTWTCTAAT" "$OUT/README" && echo "OK README primer set GTGYCAGCMGCCGCGGTAA / GGACTACNVGGGTWTCTAAT" \
  || fail "README primer set not as recorded"

[[ -f "$KEY" ]] || fail "$KEY missing"
fpr=$(gpg --batch --with-colons --show-keys "$KEY" | awk -F: '$1=="fpr"{print $10; exit}')
[[ "$fpr" == "$EXP_FPR" ]] || fail "key fingerprint $fpr != $EXP_FPR"
gpg --batch --import "$KEY" 2>&1 | tail -1
echo "OK key fingerprint $fpr"

sig="not checked (QIIME 2 environment not active)"
if command -v qiime >/dev/null; then
  qiime tools signature-verify --input-path "$OUT/$QZA" --name SILVA-Pipeline-Official && sig="verified"
  qiime tools peek "$OUT/$QZA" | head -5
fi
echo "Signature: $sig"

mkdir -p metadata
printf "name\trelease\tsource\tfile\tbytes\tmd5\tdoi\tprimers\ttrained_with\tsignature\tsigning_key_fingerprint\trecorded\n" > metadata/database_manifest.tsv
printf "SILVA SSU Ref NR99 uniform classifier V4-515f-806r\t144\t%s\t%s\t%s\t%s\t%s\tGTGYCAGCMGCCGCGGTAA/GGACTACNVGGGTWTCTAAT\tQIIME 2 2026.7.0\t%s\t%s\t%s\n" \
  "$BASE" "$QZA" "$got" "$md5" "$(cat "$OUT/$QZA.doi")" "$sig" "$fpr" "$(date -I)" >> metadata/database_manifest.tsv
echo "Written: metadata/database_manifest.tsv"

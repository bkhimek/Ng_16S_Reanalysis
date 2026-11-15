#!/usr/bin/env python3
"""21_fastq_manifest.py - check and record the per-sample FASTQs (D9, Stage 2).

Called by 21_export_fastqs.sh. For each run directory exported by QIIME 2
(MANIFEST: sample-id,filename,direction):
  - keeps the samples Stage 1 denoised: all samples in the run minus
    metadata/excluded_samples.tsv (D19); expected 451 in total;
  - sample IDs must suit nf-core/ampliseq (start with a letter; letters,
    digits, underscores only);
  - counts reads in both files (must be equal) and compares with the Stage 1
    DADA2 'input' count of the sample (results/dada2/trunc200_150/<run>_stats);
  - SHA-256 of every file.
Writes metadata/per_sample_fastq_sha256.tsv and config/ampliseq/samplesheet.tsv.
Exit code 1 if any check fails.
Usage: python3 scripts/21_fastq_manifest.py data/fastq_per_sample
"""
import csv
import gzip
import hashlib
import re
import sys
from pathlib import Path

RUNS = ["hum_run1", "hum_run2", "hum_run3"]
EXPECTED = 451
EXCL = Path("metadata/excluded_samples.tsv")
STATS = Path("results/dada2/trunc200_150")
OUT_SHA = Path("metadata/per_sample_fastq_sha256.tsv")
OUT_SHEET = Path("config/ampliseq/samplesheet.tsv")
ID_OK = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


def reads_and_sha(path):
    h = hashlib.sha256(path.read_bytes()).hexdigest()
    n = 0
    with gzip.open(path, "rb") as fh:
        for n, _ in enumerate(fh, 1):
            pass
    if n % 4:
        raise ValueError(f"{path}: line count {n} not a multiple of 4")
    return n // 4, h


def stage1_input(run):
    p = STATS / f"{run}_stats" / "stats.tsv"
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["sample-id"].startswith("#"):
                continue
            out[r["sample-id"]] = int(float(r["input"]))
    return out


def main(outd):
    ok = True
    with open(EXCL) as fh:
        excl = {r["sample-id"] for r in csv.DictReader(fh, delimiter="\t")}
    rows, sheet = [], []
    for run in RUNS:
        d = Path(outd) / run
        files = {}
        with open(d / "MANIFEST") as fh:
            for r in csv.DictReader(l for l in fh if not l.startswith("#")):
                files.setdefault(r["sample-id"], {})[r["direction"]] = d / r["filename"]
        s1 = stage1_input(run)
        keep = sorted(s for s in files if s not in excl)
        n_bad = 0
        for s in keep:
            if not ID_OK.match(s):
                print(f"FAIL {run} {s}: sample ID not valid for ampliseq"); ok = False
            f, r = files[s].get("forward"), files[s].get("reverse")
            if f is None or r is None:
                print(f"FAIL {run} {s}: forward or reverse file missing"); ok = False; continue
            nf, hf = reads_and_sha(f)
            nr, hr = reads_and_sha(r)
            if nf != nr:
                print(f"FAIL {run} {s}: forward {nf} != reverse {nr} reads"); ok = False
            if s not in s1:
                print(f"FAIL {run} {s}: not in Stage 1 DADA2 stats"); ok = False
            elif s1[s] != nf:
                n_bad += 1
                print(f"FAIL {run} {s}: {nf} read pairs exported, Stage 1 DADA2 input {s1[s]}"); ok = False
            rows += [(s, run, "forward", str(f), nf, hf), (s, run, "reverse", str(r), nr, hr)]
            sheet.append((s, str(f.resolve()), str(r.resolve()), run))
        missing = set(s1) - set(keep)
        if missing:
            print(f"FAIL {run}: {len(missing)} Stage 1 samples not exported: {sorted(missing)[:5]}"); ok = False
        print(f"{run}: {len(keep)} samples exported (excluded here: {len(set(files) & excl)}); "
              f"read pairs {sum(r[4] for r in rows if r[1] == run and r[2] == 'forward'):,}; "
              f"mismatches with Stage 1 input: {n_bad}")
    if len(sheet) != EXPECTED:
        print(f"FAIL: {len(sheet)} samples, expected {EXPECTED}"); ok = False
    OUT_SHA.parent.mkdir(exist_ok=True)
    with open(OUT_SHA, "w") as fh:
        fh.write("sample\trun\tdirection\tpath\treads\tsha256\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    OUT_SHEET.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_SHEET, "w") as fh:
        fh.write("sample\tfastq_1\tfastq_2\trun\n")
        for r in sheet:
            fh.write("\t".join(r) + "\n")
    print(f"Total: {len(sheet)} samples, {sum(r[4] for r in rows if r[2] == 'forward'):,} read pairs")
    print(f"Wrote {OUT_SHA} and {OUT_SHEET}")
    print("ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])

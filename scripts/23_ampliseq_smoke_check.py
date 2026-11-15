#!/usr/bin/env python3
"""23_ampliseq_smoke_check.py - technical checks of an ampliseq run (D30).

Not a comparison of results (those follow docs/STAGE2_COMPARISON_SPEC_v1.md
after the full run). Checks:
  1. DADA2 arguments actually used (published <out>/dada2/args/*.args.txt)
     contain the matched settings of docs/STAGE2_PARAMETER_MATCH.md.
  2. DADA2 input read pairs per sample equal the exported FASTQ counts
     (metadata/per_sample_fastq_sha256.tsv), i.e. ampliseq read every read.
  3. ASV sequences: none empty, none containing N (a concatenated pair would
     carry a run of N), lengths within the Stage 1 range.
  4. DADA2 version used by ampliseq (pipeline_info versions file).
Usage (repository root):
  python3 scripts/23_ampliseq_smoke_check.py [results/ampliseq_smoke]
"""
import csv
import re
import statistics
import sys
from pathlib import Path

REQUIRED = {
    "filterAndTrim": ["truncLen = c(200, 150)", "maxEE = c(2, 2)", "truncQ = 2", "rm.phix = TRUE", "maxN = 0", "trimLeft = 0"],
    "learnErrors": ["nreads = 1e6", "randomize = FALSE", "nbases = 1e8"],
    "dada": ["pool = FALSE"],
    "mergePairs": ["minOverlap = 12", "maxMismatch = 0", "trimOverhang = FALSE"],
    "removeBimeraDenovo": ['method = "consensus"', "minFoldParentOverAbundance = 1", "allowOneOff = FALSE"],
}
STAGE1_FASTA = Path("results/dada2/trunc200_150/repseqs_export/dna-sequences.fasta")


def norm(s):
    return re.sub(r"\s+", " ", s.replace("\\", "")).strip()


def fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = ""
        elif name is not None:
            seqs[name] += line
    return seqs


def main(out):
    out = Path(out)
    ok = True
    # 1. arguments
    args_dir = out / "dada2" / "args"
    for step, needles in REQUIRED.items():
        f = args_dir / f"{step}.args.txt"
        if not f.exists():
            print(f"FAIL args: {f} missing"); ok = False; continue
        text = norm(f.read_text())
        missing = [n for n in needles if norm(n) not in text]
        if missing:
            print(f"FAIL args {step}: not found {missing}\n     file says: {text[:300]}"); ok = False
        else:
            print(f"OK   args {step}: {', '.join(needles)}")
    # 2. input reads
    exported = {}
    with open("metadata/per_sample_fastq_sha256.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["direction"] == "forward":
                exported[r["sample"]] = int(r["reads"])
    if not (out / "dada2" / "DADA2_stats.tsv").exists():
        print(f"FAIL {out}/dada2/DADA2_stats.tsv missing: the run did not finish (see the run log)")
        print("CHECKS FAILED"); sys.exit(1)
    with open(out / "dada2" / "DADA2_stats.tsv") as fh:
        stats = list(csv.DictReader(fh, delimiter="\t"))
    bad = [(r["sample"], r["DADA2_input"], exported.get(r["sample"])) for r in stats
           if exported.get(r["sample"]) != int(float(r["DADA2_input"]))]
    if bad:
        print(f"FAIL input reads differ from exported FASTQ counts: {bad[:5]}"); ok = False
    else:
        print(f"OK   input reads: {len(stats)} samples, DADA2 input = exported read pairs for every sample")
    cols = [c for c in stats[0] if c != "sample"]
    for r in stats[:3]:
        print("     " + r["sample"] + ": " + ", ".join(f"{c} {r[c]}" for c in cols))
    # 3. sequences
    seqs = fasta(out / "dada2" / "ASV_seqs.fasta")
    lens = [len(s) for s in seqs.values()]
    s1 = [len(s) for s in fasta(STAGE1_FASTA).values()] if STAGE1_FASTA.exists() else []
    empty = sum(1 for s in seqs.values() if not s)
    with_n = sum(1 for s in seqs.values() if "N" in s.upper())
    print(f"     ASVs {len(seqs)}; length min {min(lens)}, median {statistics.median(lens):.0f}, max {max(lens)}"
          + (f"; Stage 1 range {min(s1)}-{max(s1)} (median {statistics.median(s1):.0f})" if s1 else ""))
    if empty or with_n:
        print(f"FAIL sequences: {empty} empty, {with_n} containing N"); ok = False
    elif s1 and (min(lens) < min(s1) - 5 or max(lens) > max(s1) + 5):
        print("FAIL sequences: lengths outside the Stage 1 range (+-5)"); ok = False
    else:
        print("OK   sequences: none empty, none with N, lengths within the Stage 1 range")
    # 4. versions
    vf = sorted((out / "pipeline_info").glob("*versions*.yml"))
    if vf:
        lines = [l.strip() for l in vf[0].read_text().splitlines() if re.search(r"dada2|^ *R:|nextflow", l, re.I)]
        print("     versions (" + vf[0].name + "): " + "; ".join(sorted(set(lines))[:6]))
    print("ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    if len(sys.argv) > 2:
        sys.exit(__doc__)
    main(sys.argv[1] if len(sys.argv) == 2 else "results/ampliseq_smoke")

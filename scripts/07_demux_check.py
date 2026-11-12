#!/usr/bin/env python3
"""07_demux_check.py - how did demultiplexing go, per run and per sample?

Reads the exported error-correction details of each run (one row per read pair:
assigned sample, raw and corrected barcode, number of corrected errors) and the
humanised sample metadata. Reports:
  - per run: read pairs, share assigned to a workbook sample, and how many of
    the assigned reads needed Golay correction (0 / 1 / 2 / 3 bit errors);
  - unassigned reads split into "valid codeword, not one of our samples"
    (expected: these runs carried other samples) and "not decodable";
  - per sample: read pairs, summarised by run, diet and source, and the list of
    samples below 1,000 and 5,000 read pairs (inputs to the depth decision).
Writes results/demux/per_sample_read_pairs.tsv and demux_check_report.txt.

Usage (repository root):
  python3 scripts/07_demux_check.py [results/demux] [metadata/humanised_sample_metadata.tsv]
"""
import csv
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

RUNS = ["hum_run1", "hum_run2", "hum_run3"]


def main(res_dir="results/demux", meta_path="metadata/humanised_sample_metadata.tsv"):
    res = Path(res_dir)
    meta = {}
    with open(meta_path) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if not r["sample-id"].startswith("#q2")]
    for r in rows:
        meta[r["sample-id"]] = r

    counts = Counter()
    rep = []
    for run in RUNS:
        tsvs = sorted((res / f"{run}_ecc_details").glob("*.tsv"))
        if not tsvs:
            rep.append(f"=== {run}: no error-correction details found, skipped ===\n")
            continue
        total = assigned = 0
        err = Counter()
        unassigned_valid = unassigned_undecodable = 0
        with open(tsvs[0]) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                total += 1
                if r["sample"] not in ("None", "", None):
                    assigned += 1
                    counts[r["sample"]] += 1
                    err[r["barcode-errors"]] += 1
                elif r["barcode-corrected"] in ("None", "", None):
                    unassigned_undecodable += 1
                else:
                    unassigned_valid += 1
        rep.append(f"=== {run} ===")
        rep.append(f"read pairs: {total:,}; assigned to a workbook sample: {assigned:,} ({assigned/total:.3f})")
        rep.append("  assigned reads by corrected barcode errors (bits): "
                   + ", ".join(f"{k}: {v:,} ({v/assigned:.3f})" for k, v in sorted(err.items())))
        rep.append(f"  unassigned, valid codeword not in this run's list: {unassigned_valid:,} ({unassigned_valid/total:.3f})")
        rep.append(f"  unassigned, barcode not decodable:                {unassigned_undecodable:,} ({unassigned_undecodable/total:.3f})")
        rep.append("")

    # per-sample summary
    by_group = defaultdict(list)
    out = ["sample-id\trun\tdiet\tsource\tday\tmouse_id\tcage_id\tread_pairs"]
    for sid, r in meta.items():
        n = counts.get(sid, 0)
        by_group[(r["run"], r["diet"], r["source"])].append(n)
        out.append("\t".join([sid, r["run"], r["diet"], r["source"], r["day"], r["mouse_id"], r["cage_id"], str(n)]))
    (res / "per_sample_read_pairs.tsv").write_text("\n".join(out) + "\n")

    rep.append("=== read pairs per sample (min / median / max, n) ===")
    for k in sorted(by_group):
        v = by_group[k]
        rep.append(f"  {k[0]} {k[1]:3s} {k[2]:6s}: {min(v):>7,} / {int(statistics.median(v)):>7,} / {max(v):>7,}  (n={len(v)})")
    allv = [counts.get(s, 0) for s in meta]
    for thr in (1000, 5000):
        low = sorted((counts.get(s, 0), s) for s in meta if counts.get(s, 0) < thr)
        rep.append(f"  samples below {thr:,} read pairs: {len(low)} of {len(allv)}"
                   + (": " + ", ".join(f"{s}({n})" for n, s in low[:25]) + (" ..." if len(low) > 25 else "") if low else ""))
    extra = set(counts) - set(meta)
    if extra:
        rep.append(f"  WARNING: reads assigned to IDs not in metadata: {sorted(extra)[:10]}")

    (res / "demux_check_report.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))
    print(f"Written: {res / 'demux_check_report.txt'} and per_sample_read_pairs.tsv")


if __name__ == "__main__":
    main(*sys.argv[1:])

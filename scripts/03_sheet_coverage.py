#!/usr/bin/env python3
"""03_sheet_coverage.py - which barcode sheet describes each run?

If a sheet truly describes a run, nearly every barcode listed in it should
appear in that run's index reads (every listed sample was sequenced). A sheet
that merely shares some barcodes with the run will show only the shared ones.

Uses ALL complete index reads in each partial index file (not only 5,000).

For each Figure 3 run x each of the 13 sheets:
  barcodes listed, barcodes observed (>= MIN_READS reads), fraction observed,
  share of index reads explained.
Also: share of index reads containing N, and a per-barcode table for the run's
own sheet (written to TSV) so missing samples can be identified.

Usage (repository root):
  python3 scripts/03_sheet_coverage.py data/readcheck /path/to/mapping_files.xlsx \
      results/readcheck
"""
import sys
from collections import Counter
from importlib import import_module
from pathlib import Path

import openpyxl

sys.path.insert(0, str(Path(__file__).resolve().parent))
rc_mod = import_module("01_read_check")
read_partial_gz, RUNS = rc_mod.read_partial_gz, rc_mod.RUNS

MIN_READS = 3
L = 12


def sheets_with_ids(xlsx):
    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    out = {}
    for ws in wb.worksheets:
        rows = ws.iter_rows(values_only=True)
        try:
            hdr = [str(c) for c in next(rows)]
        except StopIteration:
            continue
        if "BarcodeSequence" not in hdr:
            continue
        bi, si = hdr.index("BarcodeSequence"), hdr.index("#SampleID")
        out[ws.title] = {str(r[bi]).strip().upper(): str(r[si]) for r in rows if r[bi]}
    return out


def main(rc_dir, xlsx, out_dir):
    rc_dir, out_dir = Path(rc_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sheets = sheets_with_ids(xlsx)
    rep, per_bc = [], ["run\tsheet\tSampleID\tbarcode\treads"]

    for run, (prefix, suffixes, own) in RUNS.items():
        idx = read_partial_gz(rc_dir / f"{prefix}_{suffixes[1]}.head.fastq.gz")
        n = len(idx)
        counts = Counter(s[:L] for _, s in idx)
        with_n = sum(c for k, c in counts.items() if "N" in k) / n

        rep.append(f"=== {run}  (index reads analysed: {n}; containing N: {with_n:.3f}) ===")
        rep.append(f"    {'sheet':32s} {'listed':>6s} {'observed':>8s} {'frac_obs':>8s} {'reads_expl':>10s}")
        res = []
        for sh, bcs in sheets.items():
            obs = sum(1 for b in bcs if counts[b] >= MIN_READS)
            expl = sum(counts[b] for b in bcs) / n
            res.append((sh, len(bcs), obs, obs / len(bcs), expl))
        for sh, nl, ob, fo, ex in sorted(res, key=lambda r: -r[3]):
            mark = "  <- own sheet" if sh == own else ""
            rep.append(f"    {sh:32s} {nl:6d} {ob:8d} {fo:8.3f} {ex:10.3f}{mark}")
        own_counts = sorted(((counts[b], sid, b) for b, sid in sheets[own].items()), reverse=True)
        zero = sum(1 for c, _, _ in own_counts if c < MIN_READS)
        rep.append(f"    own sheet: {zero} of {len(own_counts)} listed samples have < {MIN_READS} reads; "
                   f"median reads per listed sample: {own_counts[len(own_counts)//2][0]}")
        rep.append("")
        for c, sid, b in own_counts:
            per_bc.append(f"{run}\t{own}\t{sid}\t{b}\t{c}")

    (out_dir / "sheet_coverage_report.txt").write_text("\n".join(rep))
    (out_dir / "own_sheet_barcode_counts.tsv").write_text("\n".join(per_bc) + "\n")
    print("\n".join(rep))
    print(f"Written: {out_dir / 'sheet_coverage_report.txt'} and own_sheet_barcode_counts.tsv")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])

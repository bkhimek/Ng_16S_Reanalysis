#!/usr/bin/env python3
"""02_barcode_census.py - follow-up to 01_read_check.py.

Uses the same partial .head.fastq.gz files (no new download).

Per Figure 3 run it reports:
  1. Flowcell  - instrument, run number and flowcell ID from the Illumina read
                 headers. One MiSeq run = one flowcell, so this shows which
                 files came from the same physical run.
  2. Census    - share of index reads matching ANY barcode in ALL 13 sheets of
                 the workbook, hits per sheet (all experiments, not only
                 Figure 3), and the 15 most common index sequences with the
                 sheets that list them.
  3. V4 check  - most common first 10 bases of forward and reverse reads (a
                 positive check that reads are 16S V4 when no primer is seen).
  4. Read-through - share of reads containing the reverse complement of the
                 opposite primer (i.e. the read ran past the end of the
                 amplicon), and the median position where it starts.

Usage (repository root):
  python3 scripts/02_barcode_census.py data/readcheck /path/to/mapping_files.xlsx \
      results/readcheck
"""
import re
import statistics
import sys
from collections import Counter
from pathlib import Path

import openpyxl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module  # noqa: E402

rc_mod = import_module("01_read_check")
read_partial_gz, rc, RUNS, primer_regex = (rc_mod.read_partial_gz, rc_mod.rc,
                                           rc_mod.RUNS, rc_mod.primer_regex)

N_READS = 5000
F515 = ["GTGCCAGCMGCCGCGGTAA", "GTGYCAGCMGCCGCGGTAA"]
R806 = ["GGACTACVSGGGTATCTAAT", "GGACTACNVGGGTWTCTAAT"]


def rc_iupac(p):
    comp = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
    return p.translate(comp)[::-1]


def all_sheet_barcodes(xlsx):
    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    out = {}
    for ws in wb.worksheets:
        rows = ws.iter_rows(values_only=True)
        try:
            hdr = [str(c) for c in next(rows)]
        except StopIteration:
            continue
        if "BarcodeSequence" not in hdr:
            continue  # README sheet
        bi = hdr.index("BarcodeSequence")
        out[ws.title] = {str(r[bi]).strip().upper() for r in rows if r[bi]}
    return out


def header_fields(h):
    # Illumina: @instrument:run:flowcell:lane:tile:x:y  [read:filtered:control:index]
    f = h[1:].split()[0].split(":")
    return tuple(f[:3]) if len(f) >= 7 else ("?", "?", "?")


def readthrough(recs, primers_rc):
    rxs = [primer_regex(p) for p in primers_rc]
    pos = []
    for _, s in recs:
        for rx in rxs:
            m = rx.search(s)
            if m:
                pos.append(m.start())
                break
    frac = len(pos) / len(recs)
    return frac, (int(statistics.median(pos)) if pos else None)


def main(rc_dir, xlsx, out_dir):
    rc_dir, out_dir = Path(rc_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sheets = all_sheet_barcodes(xlsx)
    union = set().union(*sheets.values())
    where = {}
    for sh, bcs in sheets.items():
        for b in bcs:
            where.setdefault(b, []).append(sh)
    rep = [f"Workbook: {len(sheets)} barcode sheets, {len(union)} distinct barcodes.", ""]

    for run, (prefix, suffixes, own) in RUNS.items():
        files = [rc_dir / f"{prefix}_{s}.head.fastq.gz" for s in suffixes]
        recs = [read_partial_gz(f) for f in files]
        n = min(N_READS, *(len(r) for r in recs))
        fwd, idx, rev = (r[:n] for r in recs)
        L = 12
        keys = [s[:L] for _, s in idx]

        flow = Counter(header_fields(h) for _, h in [(None, r[0]) for r in fwd])
        any_hit = sum(1 for k in keys if k in union) / n
        per_sheet = {sh: sum(1 for k in keys if k in bcs) for sh, bcs in sheets.items()}
        top = Counter(keys).most_common(15)
        f10 = Counter(s[:10] for _, s in fwd).most_common(5)
        r10 = Counter(s[:10] for _, s in rev).most_common(5)
        rt_f = readthrough(fwd, [rc_iupac(p) for p in R806])   # forward read runs into rc(806R)
        rt_r = readthrough(rev, [rc_iupac(p) for p in F515])   # reverse read runs into rc(515F)

        rep.append(f"=== {run}  ({prefix}; own sheet {own}) ===")
        rep.append("flowcell (instrument, run no., flowcell ID): "
                   + "; ".join(f"{':'.join(k)} x{v}" for k, v in flow.items()))
        rep.append(f"index reads matching ANY barcode in the workbook: {any_hit:.3f}")
        rep.append("hits per sheet (all experiments):")
        for sh, c in sorted(per_sheet.items(), key=lambda kv: -kv[1]):
            mark = "  <- own sheet" if sh == own else ""
            rep.append(f"    {sh:32s} {c/n:.3f}{mark}")
        rep.append("15 most common index sequences (count, share, sheets listing it):")
        for k, c in top:
            rep.append(f"    {k}  {c:5d}  {c/n:.3f}  {', '.join(where.get(k, ['(in no sheet)']))}")
        rep.append("most common first 10 bases, forward: " + ", ".join(f"{k} {c/n:.2f}" for k, c in f10))
        rep.append("most common first 10 bases, reverse: " + ", ".join(f"{k} {c/n:.2f}" for k, c in r10))
        rep.append(f"read-through: forward reads containing rc(806R) {rt_f[0]:.3f} (median start {rt_f[1]}); "
                   f"reverse reads containing rc(515F) {rt_r[0]:.3f} (median start {rt_r[1]})")
        rep.append("")

    (out_dir / "barcode_census_report.txt").write_text("\n".join(rep))
    print("\n".join(rep))
    print(f"Written: {out_dir / 'barcode_census_report.txt'}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])

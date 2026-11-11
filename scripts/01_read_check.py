#!/usr/bin/env python3
"""01_read_check.py - inspect the first reads of each Figure 3 run.

Input: the partial .head.fastq.gz files made by 00_fetch_read_heads.sh and the
depositors' mapping_files.xlsx. Nothing here trusts the file names: each check
is made from the reads themselves.

For each of the 6 runs it reports:
  1. Pairing   - do the read IDs agree, record by record, across the 3 files?
  2. Read role - read lengths per file; the short file is the index read.
  3. Barcodes  - fraction of index reads matching the sheet barcodes, as given
                 and reverse-complemented, against all 6 run sheets (tests the
                 assumption that sheet runN = file prepN).
  4. Primers   - fraction of forward/reverse reads that start with (or contain
                 near the start) the 515F and 806R primer variants.

Usage (repository root):
  python3 scripts/01_read_check.py data/readcheck /path/to/mapping_files.xlsx \
      results/readcheck
"""
import re
import statistics
import sys
import zlib
from pathlib import Path

import openpyxl

N_READS = 5000  # records compared per run (fewer if a piece holds fewer)

RUNS = {
    # run label: (file prefix, (forward, index, reverse) suffixes per README, sheet)
    "hum_run1": ("Figure_3_humanized_double_ciprofloxacin_prep1", ("R1", "R2", "R3"), "Fig_3_hum_double_cipro_run1"),
    "hum_run2": ("Figure_3_humanized_double_ciprofloxacin_prep2", ("R1", "R2", "R3"), "Fig_3_hum_double_cipro_run2"),
    "hum_run3": ("Figure_3_humanized_double_ciprofloxacin_prep3", ("R1_001", "I1_001", "R2_001"), "Fig_3_hum_double_cipro_run3"),
    "conv_run1": ("Figure_3_conventional_ciprofloxacin_prep_1", ("R1", "R2", "R3"), "Fig_3_conv_cipro_run1"),
    "conv_run2": ("Figure_3_conventional_ciprofloxacin_prep_2", ("R1", "R2", "R3"), "Fig_3_conv_cipro_run2"),
    "conv_run3": ("Figure_3_conventional_ciprofloxacin_prep_3", ("R1", "R2", "R3"), "Fig_3_conv_cipro_run3"),
}

# 16S V4 primer variants: original (Caporaso 2011) and modified (Parada 2016 / Apprill 2015)
PRIMERS = {
    "515F_orig": "GTGCCAGCMGCCGCGGTAA",
    "515F_Parada": "GTGYCAGCMGCCGCGGTAA",
    "806R_orig": "GGACTACVSGGGTATCTAAT",
    "806R_Apprill": "GGACTACNVGGGTWTCTAAT",
}
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "S": "[CG]",
         "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]",
         "V": "[ACG]", "N": "[ACGT]"}
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def primer_regex(p):
    return re.compile("".join(IUPAC[c] for c in p))


def read_partial_gz(path):
    """Decompress a truncated gzip file; return complete FASTQ records only."""
    raw = Path(path).read_bytes()
    text = bytearray()
    while raw:
        d = zlib.decompressobj(wbits=31)
        try:
            text += d.decompress(raw)
        except zlib.error:
            break  # truncated mid-block: keep what was decoded
        if not d.eof:
            break  # stream cut off (expected for a partial download)
        raw = d.unused_data  # next gzip member, if the file has several
    lines = text.decode("ascii", errors="replace").split("\n")
    lines = lines[: (len(lines) - 1)]  # last line may be incomplete
    recs = []
    for i in range(0, len(lines) - 3, 4):
        h, s, plus, q = lines[i: i + 4]
        if not h.startswith("@") or not plus.startswith("+") or len(s) != len(q):
            break  # malformed: stop rather than misread
        recs.append((h, s))
    return recs


def read_id(header):
    tok = header[1:].split()[0]
    return re.sub(r"/[123]$", "", tok)


def load_barcodes(xlsx):
    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    out = {}
    for _, (_, _, sheet) in RUNS.items():
        ws = wb[sheet]
        rows = ws.iter_rows(values_only=True)
        hdr = [str(c) for c in next(rows)]
        bi, li = hdr.index("BarcodeSequence"), hdr.index("LinkerPrimerSequence")
        bcs, linkers = set(), set()
        for r in rows:
            if r[bi]:
                bcs.add(str(r[bi]).strip().upper())
                linkers.add(str(r[li]).strip().upper())
        out[sheet] = (bcs, linkers)
    return out


def length_summary(recs):
    ls = [len(s) for _, s in recs]
    return f"{min(ls)}/{int(statistics.median(ls))}/{max(ls)}"


def primer_hits(recs):
    res = {}
    for name, p in PRIMERS.items():
        rx = primer_regex(p)
        at0 = sum(1 for _, s in recs if rx.match(s))
        near = sum(1 for _, s in recs if rx.search(s[:40]))
        res[name] = (at0 / len(recs), near / len(recs))
    return res


def main(rc_dir, xlsx, out_dir):
    rc_dir, out_dir = Path(rc_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    barcodes = load_barcodes(xlsx)
    report, rows = [], []

    for run, (prefix, suffixes, sheet) in RUNS.items():
        files = [rc_dir / f"{prefix}_{s}.head.fastq.gz" for s in suffixes]
        recs = [read_partial_gz(f) for f in files]
        n = min(N_READS, *(len(r) for r in recs))
        recs = [r[:n] for r in recs]
        fwd, idx, rev = recs

        ids_agree = sum(1 for a, b, c in zip(*recs) if read_id(a[0]) == read_id(b[0]) == read_id(c[0]))
        lens = [length_summary(r) for r in recs]
        med = [statistics.median(len(s) for _, s in r) for r in recs]
        shortest = suffixes[med.index(min(med))]

        L = len(next(iter(barcodes[sheet][0])))
        bc_rates = {}
        for sh, (bcs, _) in barcodes.items():
            asis = sum(1 for _, s in idx if s[:L] in bcs) / n
            rcomp = sum(1 for _, s in idx if rc(s)[:L] in bcs or rc(s[:L]) in bcs) / n
            bc_rates[sh] = (asis, rcomp)
        best_sheet = max(bc_rates, key=lambda k: max(bc_rates[k]))
        # Barcodes are reused between runs, so a high match rate alone cannot
        # identify the run. Count reads whose barcode occurs in ONE sheet only.
        # orientation chosen from the best match over ALL sheets, so that a
        # wrongly labelled run still gets its barcodes read the right way round
        orient_rc = max(b for _, b in bc_rates.values()) > max(a for a, _ in bc_rates.values())
        keys = [(rc(s[:L]) if orient_rc else s[:L]) for _, s in idx]
        unique_hits = {}
        for sh, (bcs, _) in barcodes.items():
            others = set().union(*(b for o, (b, _) in barcodes.items() if o != sh))
            uniq = bcs - others
            unique_hits[sh] = (sum(1 for k in keys if k in uniq), len(uniq))

        pf, pr = primer_hits(fwd), primer_hits(rev)

        report.append(f"=== {run}  ({prefix}; sheet {sheet}) ===")
        report.append(f"records compared: {n}; read IDs identical across 3 files: {ids_agree}/{n}")
        report.append(f"length min/median/max  {suffixes[0]}: {lens[0]}  {suffixes[1]}: {lens[1]}  {suffixes[2]}: {lens[2]}")
        report.append(f"shortest file (expected index): {shortest}  [README index: {suffixes[1]}]")
        report.append(f"barcode length in sheet: {L}; sheet linker(s): {', '.join(sorted(barcodes[sheet][1]))}")
        own = bc_rates[sheet]
        report.append(f"index reads matching OWN sheet barcodes: as-is {own[0]:.3f}, reverse-complement {own[1]:.3f}")
        report.append("  per sheet: match rate as-is / revcomp | reads with a barcode found ONLY in that sheet (n unique barcodes):")
        for sh, (a, b) in bc_rates.items():
            mark = "  <- own sheet" if sh == sheet else ""
            u, nu = unique_hits[sh]
            report.append(f"    {sh:30s} {a:.3f} / {b:.3f} | {u:5d} ({nu}){mark}")
        top_unique = max(unique_hits, key=lambda k: unique_hits[k][0])
        own_hits, own_n = unique_hits[sheet]
        if own_n == 0:
            verdict = "NOT TESTABLE (own sheet has no sheet-unique barcodes; rely on match rates)"
        elif top_unique == sheet and own_hits > 0:
            verdict = "SUPPORTED"
        else:
            verdict = "NOT SUPPORTED"
        report.append(f"  sheet with most sheet-unique barcode hits: {top_unique} -> run mapping {verdict}")
        report.append("primer at read start / within first 40 bp:")
        for name in PRIMERS:
            report.append(f"    {name:14s} forward {pf[name][0]:.3f} / {pf[name][1]:.3f}   reverse {pr[name][0]:.3f} / {pr[name][1]:.3f}")
        report.append("")
        rows.append([run, sheet, n, ids_agree, *lens, shortest, L, f"{own[0]:.3f}", f"{own[1]:.3f}",
                     unique_hits[sheet][0], top_unique, verdict.split(" (")[0]]
                    + [f"{pf[k][1]:.3f}" for k in PRIMERS] + [f"{pr[k][1]:.3f}" for k in PRIMERS])

    hdr = ["run", "sheet", "n_records", "ids_agree", "len_fwd", "len_index", "len_rev", "shortest_file",
           "barcode_len", "bc_match_asis", "bc_match_revcomp", "own_unique_hits", "top_unique_sheet", "run_mapping"] \
        + [f"fwd_{k}_first40" for k in PRIMERS] + [f"rev_{k}_first40" for k in PRIMERS]
    (out_dir / "read_check_summary.tsv").write_text(
        "\t".join(hdr) + "\n" + "\n".join("\t".join(map(str, r)) for r in rows) + "\n")
    (out_dir / "read_check_report.txt").write_text("\n".join(report))
    print("\n".join(report))
    print(f"Written: {out_dir / 'read_check_summary.tsv'} and read_check_report.txt")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])

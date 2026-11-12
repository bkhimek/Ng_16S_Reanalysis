#!/usr/bin/env python3
"""05_make_sample_metadata.py - QIIME 2 metadata for the humanised experiment.

Reads the three humanised Figure 3 sheets of the depositors' workbook and writes:

  metadata/humanised_sample_metadata.tsv
      one row per sample (all 3 runs), QIIME 2 metadata format. Sample IDs are
      experiment-qualified ("H" + original ID) so they can never collide with
      conventional samples and never look numeric to R or QIIME 2. Mouse and
      cage IDs are qualified the same way ("H_m1", "H_c1") so they stay
      categorical. The original workbook ID is kept in its own column.
  metadata/demux/hum_run{1,2,3}_barcodes.tsv
      sample-id + barcode-sequence for ONE run. Barcodes are reused between
      runs, so each run must be demultiplexed with its own barcode list only.

Checks (the script stops if any fails): 12-nt A/C/G/T barcodes, no duplicate
barcode within a run, no duplicate sample ID, expected sample counts.

Usage (repository root):
  python3 scripts/05_make_sample_metadata.py /path/to/mapping_files.xlsx
"""
import re
import sys
from collections import Counter
from pathlib import Path

import openpyxl

SHEETS = {"hum_run1": "Fig_3_hum_double_cipro_run1",
          "hum_run2": "Fig_3_hum_double_cipro_run2",
          "hum_run3": "Fig_3_hum_double_cipro_run3"}
EXPECTED = {"hum_run1": 62, "hum_run2": 364, "hum_run3": 27}


def fmt_day(v):
    f = float(v)
    return str(int(f)) if f.is_integer() else str(f)


def main(xlsx):
    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    out_rows, problems = [], []
    Path("metadata/demux").mkdir(parents=True, exist_ok=True)

    for run, sheet in SHEETS.items():
        it = wb[sheet].iter_rows(values_only=True)
        hdr = [str(c) for c in next(it)]
        col = {name: hdr.index(name) for name in
               ["#SampleID", "BarcodeSequence", "day", "mouse", "chow", "cage", "source", "antibiotic"]}
        rows = [r for r in it if r[0] is not None]
        if len(rows) != EXPECTED[run]:
            problems.append(f"{run}: {len(rows)} rows, expected {EXPECTED[run]}")
        bcs = [str(r[col["BarcodeSequence"]]).strip().upper() for r in rows]
        bad = [b for b in bcs if not re.fullmatch(r"[ACGT]{12}", b)]
        dups = [b for b, c in Counter(bcs).items() if c > 1]
        if bad:
            problems.append(f"{run}: invalid barcodes {bad[:5]}")
        if dups:
            problems.append(f"{run}: duplicate barcodes within run {dups[:5]}")

        demux_lines = ["sample-id\tbarcode-sequence"]
        for r, bc in zip(rows, bcs):
            orig = str(r[col["#SampleID"]])
            sid = f"H{orig}"
            src = str(r[col["source"]])
            out_rows.append({
                "sample-id": sid,
                "barcode-sequence": bc,
                "original_sample_id": orig,
                "experiment": "humanised",
                "run": run,
                "mouse_id": f"H_m{r[col['mouse']]}",
                "cage_id": f"H_c{r[col['cage']]}",
                "diet": str(r[col["chow"]]),
                "day": fmt_day(r[col["day"]]),
                "source": "faeces" if src == "Feces" else "cecal",
                "antibiotic": str(r[col["antibiotic"]]),
            })
            demux_lines.append(f"{sid}\t{bc}")
        Path(f"metadata/demux/{run}_barcodes.tsv").write_text("\n".join(demux_lines) + "\n")

    ids = [r["sample-id"] for r in out_rows]
    if len(ids) != len(set(ids)):
        problems.append("duplicate sample IDs across runs")

    # flag repeated mouse-day-source entries (not assumed to be technical replicates)
    key = Counter((r["mouse_id"], r["day"], r["source"]) for r in out_rows)
    for r in out_rows:
        r["repeat_mouse_day"] = "yes" if key[(r["mouse_id"], r["day"], r["source"])] > 1 else "no"

    if problems:
        sys.exit("STOP - metadata problems:\n  " + "\n  ".join(problems))

    cols = list(out_rows[0].keys())
    lines = ["\t".join(cols),
             "#q2:types\t" + "\t".join("numeric" if c == "day" else "categorical" for c in cols[1:])]
    lines += ["\t".join(r[c] for c in cols) for r in out_rows]
    Path("metadata/humanised_sample_metadata.tsv").write_text("\n".join(lines) + "\n")

    print(f"Samples: {len(out_rows)}  ({', '.join(f'{k} {v}' for k, v in Counter(r['run'] for r in out_rows).items())})")
    print("By diet x source:", dict(Counter((r["diet"], r["source"]) for r in out_rows)))
    print("Mice:", len({r["mouse_id"] for r in out_rows}), " Cages:", len({r["cage_id"] for r in out_rows}))
    print("Rows flagged repeat_mouse_day:", sum(r["repeat_mouse_day"] == "yes" for r in out_rows))
    print("Written: metadata/humanised_sample_metadata.tsv, metadata/demux/hum_run{1,2,3}_barcodes.tsv")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])

#!/usr/bin/env python3
"""15_rarefaction_check.py - is 6,400 reads near the plateau? (D27, technical)

Reads the per-sample rarefaction values stored inside the .qzv made by
15_rarefaction_curves.sh (one CSV per metric; columns depth-<d>_iter-<i>).
For each sample the 10 draws at a depth are averaged. Reported per run (no
diet, no time):
  - median value at 6,400 as a share of the value at 12,800, among samples
    with at least 12,800 reads (how much more would be found at double depth);
  - median relative gain from 5,600 to 6,400 reads (local slope near the
    chosen depth; close to 0 = flat).
Interpretation: Shannon and Faith's PD usually plateau early; observed ASVs
keeps rising slowly. A share >= ~0.90 at 6,400 means the depth captures most of
the detectable richness; the remaining gap is the same for every sample
because all are rarefied to the same depth.

Usage (repository root):
  python3 scripts/15_rarefaction_check.py trunc200_150 [trunc218_166]
"""
import csv
import io
import statistics
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

METRICS = ["observed_features", "faith_pd", "shannon"]
D_MAIN, D_PREV, D_DOUBLE = 6400, 5600, 12800


def per_sample(zf, metric):
    name = next(n for n in zf.namelist() if n.endswith(f"/data/{metric}.csv"))
    rows = list(csv.DictReader(io.TextIOWrapper(zf.open(name), encoding="utf-8")))
    out = {}
    for r in rows:
        vals = defaultdict(list)
        for k, v in r.items():
            if k and k.startswith("depth-") and v not in ("", "nan", "NaN", None):
                d = int(k.split("_")[0].split("-")[1])
                vals[d].append(float(v))
        out[r["sample-id"]] = ({d: statistics.mean(v) for d, v in vals.items()}, r.get("run", "?"))
    return out


def main(arms):
    for arm in arms:
        qzv = Path("results/diversity") / arm / "alpha_rarefaction_by_run.qzv"
        rep = [f"##### Rarefaction check {arm} (depth {D_MAIN}; grouped by run only)"]
        with zipfile.ZipFile(qzv) as zf:
            for m in METRICS:
                data = per_sample(zf, m)
                share, gain = defaultdict(list), defaultdict(list)
                for s, (v, run) in data.items():
                    if D_MAIN in v and D_DOUBLE in v and v[D_DOUBLE] > 0:
                        share[run].append(v[D_MAIN] / v[D_DOUBLE])
                    if D_MAIN in v and D_PREV in v and v[D_PREV] > 0:
                        gain[run].append(v[D_MAIN] / v[D_PREV] - 1)
                rep.append(f"  {m}:")
                for run in sorted(set(share) | set(gain)):
                    sh, g = share.get(run, []), gain.get(run, [])
                    rep.append(f"    {run}: value at {D_MAIN} / value at {D_DOUBLE}: "
                               + (f"median {statistics.median(sh):.3f} (min {min(sh):.3f}, n={len(sh)})" if sh else "n/a")
                               + f"; gain {D_PREV}->{D_MAIN}: "
                               + (f"median {statistics.median(g):+.3f} (n={len(g)})" if g else "n/a"))
        out = Path("results/diversity") / arm / "rarefaction_check_report.txt"
        out.write_text("\n".join(rep) + "\n")
        print("\n".join(rep) + "\n")


if __name__ == "__main__":
    if not 1 <= len(sys.argv[1:]) <= 2:
        sys.exit(__doc__)
    main(sys.argv[1:])

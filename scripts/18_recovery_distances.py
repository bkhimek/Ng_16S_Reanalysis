#!/usr/bin/env python3
"""18_recovery_distances.py - distance-from-baseline analyses of the frozen
specification (docs/ANALYSIS_SPEC_FROZEN_v1.md, tag analysis-frozen-v1).

Samples: metadata/primary_analysis_samples.tsv (script 17).
For every analysis: per mouse, mean distance of its window samples from its own
day-0 sample; per cage, mean of its mice; per diet, unweighted mean of its cages;
delta = mean(PD cages) - mean(STD cages). Positive delta = PD further from
baseline (read as less recovery). No p-value (spec section 3).

  PRIMARY  Bray-Curtis, 200/150 rarefied medoid matrix, days 8-12
  SEC-wU   weighted UniFrac, 200/150 rarefied medoid matrix, days 8-12
  S1       Bray-Curtis, 218/166 rarefied medoid matrix (D21; both deltas and difference)
  S2       Bray-Curtis on relative abundances, unrarefied 200/150 table
  S3       Aitchison: ASVs with >=1 read in >=10% of the 102 primary samples,
           pseudocount 1, CLR per sample, Euclidean distance
  S4       primary matrix, window days 8-13
Supplement: delta under all 10 allocations of the 5 cages into 3 STD + 2 PD
(descriptive position only, not a test).

Before computing, every input artifact's MD5 is checked against
docs/frozen_artifacts.tsv.
Outputs: results/recovery/ (report, mouse and cage values, allocations, figure).
Usage (repository root, QIIME 2 environment active for numpy/matplotlib):
  python3 scripts/18_recovery_distances.py
"""
import csv
import hashlib
import io
import itertools
import math
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

PRIMARY, SENS = "trunc200_150", "trunc218_166"
OUT = Path("results/recovery")
SAMPLES = Path("metadata/primary_analysis_samples.tsv")
FROZEN = Path("docs/frozen_artifacts.tsv")
PREV = 0.10  # S3 prevalence filter
DIETS = ("STD", "PD")


# ---------- inputs ----------
def load_samples():
    with open(SAMPLES) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    mice = {}
    for r in rows:
        m = mice.setdefault(r["mouse_id"], {"cage": r["cage_id"], "diet": r["diet"], "days": {}})
        m["days"][float(r["day"])] = r["sample-id"]
    return mice


def frozen_md5():
    with open(FROZEN) as fh:
        return {r["path"]: r["md5"] for r in csv.DictReader(fh, delimiter="\t")}


def check_md5(path, frozen):
    want = frozen.get(str(path))
    got = hashlib.md5(Path(path).read_bytes()).hexdigest()
    if want is None:
        sys.exit(f"STOP: {path} is not listed in {FROZEN}")
    if got != want:
        sys.exit(f"STOP: {path} md5 {got} differs from frozen {want}")
    return got


def find_dm(arm, names, frozen):
    cands = [Path(p) for p in frozen if p.startswith(f"data/qiime/diversity/{arm}/distance_matrices/")
             and Path(p).stem in names]
    if len(cands) != 1:
        sys.exit(f"STOP: expected one frozen {names} matrix for {arm}, found {cands}")
    check_md5(cands[0], frozen)
    return cands[0]


def read_dm(path, keep):
    """Distance matrix rows/columns for the samples in keep only."""
    with zipfile.ZipFile(path) as zf:
        name = next(n for n in zf.namelist() if n.endswith("/data/distance-matrix.tsv"))
        text = io.TextIOWrapper(zf.open(name), encoding="utf-8")
        header = text.readline().rstrip("\n").split("\t")[1:]
        idx = {s: i for i, s in enumerate(header)}
        missing = keep - set(idx)
        if missing:
            sys.exit(f"STOP: {path} lacks {sorted(missing)[:5]}")
        dm = {}
        for line in text:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in keep:
                dm[parts[0]] = {s: float(parts[1 + idx[s]]) for s in keep}
    return lambda a, b: dm[a][b]


def read_table(arm, keep):
    path = Path("results/dada2") / arm / "table_export" / "feature-table.tsv"
    with open(path) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("# Constructed")]
    cols = lines[0][1:]
    idx = [i for i, s in enumerate(cols) if s in keep]
    if len(idx) != len(keep):
        sys.exit(f"STOP: {path} lacks some samples")
    names = [cols[i] for i in idx]
    counts = {s: [] for s in names}
    for row in lines[1:]:
        for i, s in zip(idx, names):
            counts[s].append(float(row[1 + i]))
    return counts   # sample -> list of counts (same ASV order)


# ---------- distances from the table ----------
def bc_relative(counts):
    rel = {s: [x / sum(v) for x in v] for s, v in counts.items()}
    return lambda a, b: 0.5 * sum(abs(x - y) for x, y in zip(rel[a], rel[b]))


def aitchison(counts, base_ids):
    n = len(base_ids)
    n_asv = len(next(iter(counts.values())))
    keep = [j for j in range(n_asv)
            if sum(counts[s][j] > 0 for s in base_ids) >= math.ceil(PREV * n)]
    clr = {}
    for s, v in counts.items():
        logs = [math.log(v[j] + 1) for j in keep]
        g = sum(logs) / len(logs)
        clr[s] = [x - g for x in logs]
    f = lambda a, b: math.sqrt(sum((x - y) ** 2 for x, y in zip(clr[a], clr[b])))
    return f, len(keep), n_asv


# ---------- the frozen estimate ----------
def estimate(mice, dist, days):
    mouse = {}
    for m, d in mice.items():
        base = d["days"][0.0]
        vals = [dist(d["days"][float(x)], base) for x in days]
        mouse[m] = sum(vals) / len(vals)
    cage = defaultdict(list)
    for m, v in mouse.items():
        cage[mice[m]["cage"]].append(v)
    cage_mean = {c: sum(v) / len(v) for c, v in cage.items()}
    cage_n = {c: len(v) for c, v in cage.items()}
    cage_diet = {mice[m]["cage"]: mice[m]["diet"] for m in mice}
    diet_mean = {g: _mean([cage_mean[c] for c in cage_mean if cage_diet[c] == g]) for g in DIETS}
    return {"mouse": mouse, "cage": cage_mean, "cage_n": cage_n, "cage_diet": cage_diet,
            "diet": diet_mean, "delta": diet_mean["PD"] - diet_mean["STD"]}


def _mean(x):
    return sum(x) / len(x)


def allocations(res):
    cages = sorted(res["cage"])
    n_pd = sum(1 for c in cages if res["cage_diet"][c] == "PD")
    out = []
    for pd in itertools.combinations(cages, n_pd):
        d = _mean([res["cage"][c] for c in pd]) - _mean([res["cage"][c] for c in cages if c not in pd])
        observed = set(pd) == {c for c in cages if res["cage_diet"][c] == "PD"}
        out.append((d, pd, observed))
    return sorted(out, key=lambda x: -x[0])


# ---------- main ----------
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mice = load_samples()
    frozen = frozen_md5()
    win, win13 = [8, 9, 10, 11, 12], [8, 9, 10, 11, 12, 13]
    if any(len(m["days"]) < 6 for m in mice.values()):
        sys.exit("STOP: a mouse in the sample list lacks day 0 or a window day")
    ids_prim = {s for m in mice.values() for d, s in m["days"].items() if d in [0.0] + [float(x) for x in win]}
    ids_all = {s for m in mice.values() for s in m["days"].values()}
    have13 = all(13.0 in m["days"] for m in mice.values())

    for arm in (PRIMARY, SENS):
        check_md5(Path("data/qiime/dada2") / arm / "humanised_table.qza", frozen)

    analyses = []
    BC = ("braycurtis", "bray_curtis", "bray_curtis_distance_matrix")
    bc = read_dm(find_dm(PRIMARY, BC, frozen), ids_all)
    analyses.append(("PRIMARY", "Bray-Curtis, 200/150 rarefied medoid, days 8-12", estimate(mice, bc, win)))
    wu = read_dm(find_dm(PRIMARY, ("weighted_unifrac", "weighted_unifrac_distance_matrix"), frozen), ids_prim)
    analyses.append(("SEC-wU", "weighted UniFrac, 200/150 rarefied medoid, days 8-12", estimate(mice, wu, win)))
    bc2 = read_dm(find_dm(SENS, BC, frozen), ids_prim)
    analyses.append(("S1", "Bray-Curtis, 218/166 rarefied medoid, days 8-12", estimate(mice, bc2, win)))
    counts = read_table(PRIMARY, ids_prim)
    analyses.append(("S2", "Bray-Curtis on relative abundances, unrarefied 200/150, days 8-12",
                     estimate(mice, bc_relative(counts), win)))
    ait, n_keep, n_asv = aitchison(counts, sorted(ids_prim))
    analyses.append(("S3", f"Aitchison (CLR, pseudocount 1; {n_keep} of {n_asv} ASVs pass the 10% filter), unrarefied 200/150, days 8-12",
                     estimate(mice, ait, win)))
    if have13:
        analyses.append(("S4", "Bray-Curtis, 200/150 rarefied medoid, days 8-13", estimate(mice, bc, win13)))

    prim = analyses[0][2]
    cages = sorted(prim["cage"], key=lambda c: (DIETS.index(prim["cage_diet"][c]), c))
    rep = ["##### Recovery: distance from each mouse's own day-0 sample (frozen spec v1)",
           f"{len(mice)} mice, {len(cages)} cages; delta = mean(PD cages) - mean(STD cages); positive = PD further from baseline",
           "", "PRIMARY: " + analyses[0][1],
           f"  delta = {prim['delta']:+.4f}   (STD mean {prim['diet']['STD']:.4f}, PD mean {prim['diet']['PD']:.4f})",
           "  cage values (mice):"]
    for c in cages:
        rep.append(f"    {prim['cage_diet'][c]:3s} {c}: {prim['cage'][c]:.4f}  (n={prim['cage_n'][c]})")
    rep += ["", "All distance-from-baseline analyses (same mice, days, cages):",
            f"  {'ID':7s} {'delta':>8s}  {'STD':>7s}  {'PD':>7s}  " + "  ".join(f"{c:>7s}" for c in cages) + "   definition"]
    for aid, desc, r in analyses:
        rep.append(f"  {aid:7s} {r['delta']:+8.4f}  {r['diet']['STD']:7.4f}  {r['diet']['PD']:7.4f}  "
                   + "  ".join(f"{r['cage'][c]:7.4f}" for c in cages) + f"   {desc}")
    s1 = next(r for a, _, r in analyses if a == "S1")
    rep += ["", f"D21 (S1): delta 200/150 {prim['delta']:+.4f}; delta 218/166 {s1['delta']:+.4f}; "
            f"difference (218/166 - 200/150) {s1['delta'] - prim['delta']:+.4f}; "
            f"same direction: {'yes' if prim['delta'] * s1['delta'] > 0 else 'NO'}"]
    rep += ["", "Supplement: delta under all 10 cage allocations (3 STD + 2 PD), largest first; descriptive, not a test",
            "  (cage labels are not exchangeable: STD baselines are in run 1, all other samples in run 2)"]
    allo = allocations(prim)
    for i, (d, pd, obs) in enumerate(allo, 1):
        rep.append(f"  {i:2d}. {d:+.4f}  'PD' = {', '.join(pd)}" + ("   <- observed" if obs else ""))
    rep += ["", "Limit: STD distances cross sequencing runs (run 1 -> run 2); PD distances stay in run 2.",
            "Delta describes the cage-level pattern; attributing it to diet is limited by this confounding."]
    (OUT / "recovery_report.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))

    with open(OUT / "mouse_values.tsv", "w") as fh:
        fh.write("analysis\tmouse_id\tcage_id\tdiet\tmean_distance_from_day0\n")
        for aid, _, r in analyses:
            for m in sorted(r["mouse"], key=lambda x: int(x.split("m")[-1])):
                fh.write(f"{aid}\t{m}\t{mice[m]['cage']}\t{mice[m]['diet']}\t{r['mouse'][m]:.6f}\n")
    with open(OUT / "cage_values.tsv", "w") as fh:
        fh.write("analysis\tcage_id\tdiet\tn_mice\tcage_mean\tdelta_PD_minus_STD\n")
        for aid, _, r in analyses:
            for c in cages:
                fh.write(f"{aid}\t{c}\t{r['cage_diet'][c]}\t{r['cage_n'][c]}\t{r['cage'][c]:.6f}\t{r['delta']:.6f}\n")
    with open(OUT / "allocations_primary.tsv", "w") as fh:
        fh.write("rank\tdelta\tcages_labelled_PD\tobserved\n")
        for i, (d, pd, obs) in enumerate(allo, 1):
            fh.write(f"{i}\t{d:.6f}\t{','.join(pd)}\t{'yes' if obs else 'no'}\n")
    plot(prim, mice, cages)
    print(f"\nWrote {OUT}/recovery_report.txt, mouse_values.tsv, cage_values.tsv, allocations_primary.tsv, primary_cages.png")


def plot(r, mice, cages):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available: figure skipped")
        return
    col = {"STD": "#3E5C76", "PD": "#C0703A"}
    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=150)
    for i, c in enumerate(cages):
        ys = [v for m, v in r["mouse"].items() if mice[m]["cage"] == c]
        xs = [i + 0.12 * (k - (len(ys) - 1) / 2) for k in range(len(ys))]
        ax.scatter(xs, ys, s=22, color=col[r["cage_diet"][c]], alpha=0.55, zorder=2)
        ax.hlines(r["cage"][c], i - 0.3, i + 0.3, color=col[r["cage_diet"][c]], lw=2.5, zorder=3)
    for g in DIETS:
        ax.axhline(r["diet"][g], color=col[g], lw=1, ls="--", zorder=1, label=f"{g} mean of cages")
    ax.set_xticks(range(len(cages)))
    ax.set_xticklabels([f"{c}\n{r['cage_diet'][c]} (n={r['cage_n'][c]})" for c in cages])
    ax.set_ylabel("Bray-Curtis from own day 0\n(mean of days 8-12)")
    ax.set_title(f"Primary: PD - STD = {r['delta']:+.3f} (cage means)", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "primary_cages.png")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        sys.exit(__doc__)
    main()

#!/usr/bin/env python3
"""20_exploratory_window.py - EXPLORATORY analyses after the frozen results (D29).
Not part of the frozen specification (tag analysis-frozen-v1); results are
labelled exploratory. Same 17 mice and cage-equal weighting as scripts 18-19.

A. Window levels (days 8-12; every window sample of both diets is in run 2, so
   this removes the run difference between baselines; cohort concurrency and
   sample preparation remain unverified). Per mouse mean over the 5 days ->
   cage mean -> diet mean of cages. Observed ASVs and Faith's PD (rarefied,
   both truncation arms), and family shares (unrarefied 200/150).
B. Muribaculaceae detection per sample: reads > 0 in any Muribaculaceae ASV,
   counted over the eligible window samples and day-0 samples of each diet.
C. Family decomposition of the unrarefied relative-abundance Bray-Curtis from
   day 0 (S2). For one pair, BC = 0.5 * sum over ASVs |p_window - p_day0|, so
   it splits exactly into family parts (0.5 * sum over that family's ASVs).
   Parts are averaged mouse -> cage -> diet like S2; per family the PD - STD
   difference of its part; the parts add up to the S2 delta (checked).
D. Window samples only, rarefied medoid Bray-Curtis (200/150): mean distance
   between cages (pairs of different mice) and a PCoA plot coloured by diet,
   marker by cage.

Inputs checked against docs/frozen_artifacts.tsv where listed (alpha vectors,
distance matrix).
Outputs: results/exploratory/ (report, tables, figure)
Usage (repository root, QIIME 2 environment active for numpy/matplotlib):
  python3 scripts/20_exploratory_window.py
"""
import csv
import hashlib
import io
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

PRIMARY, SENS = "trunc200_150", "trunc218_166"
OUT = Path("results/exploratory")
SAMPLES = Path("metadata/primary_analysis_samples.tsv")
FROZEN = Path("docs/frozen_artifacts.tsv")
WIN = [8, 9, 10, 11, 12]
FAMILIES = ["Muribaculaceae", "Bacteroidaceae", "Akkermansiaceae"]
DIETS = ("STD", "PD")
BC_NAMES = ("braycurtis", "bray_curtis", "bray_curtis_distance_matrix")


def mean(x):
    return sum(x) / len(x)


# ---------- inputs ----------
def load_mice():
    mice = {}
    with open(SAMPLES) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m = mice.setdefault(r["mouse_id"], {"cage": r["cage_id"], "diet": r["diet"], "days": {}})
            if r["role"] in ("baseline", "window"):
                m["days"][float(r["day"])] = r["sample-id"]
    if any(len(d["days"]) != 6 for d in mice.values()):
        sys.exit("STOP: sample list does not give day 0 + days 8-12 for every mouse")
    return mice


def frozen_path(prefix, stems):
    with open(FROZEN) as fh:
        rows = {r["path"]: r["md5"] for r in csv.DictReader(fh, delimiter="\t")}
    c = [p for p in rows if p.startswith(prefix) and Path(p).stem in stems]
    if len(c) != 1:
        sys.exit(f"STOP: expected one frozen file {prefix}{stems}, found {c}")
    if hashlib.md5(Path(c[0]).read_bytes()).hexdigest() != rows[c[0]]:
        sys.exit(f"STOP: {c[0]} differs from the frozen MD5")
    return Path(c[0])


def read_alpha(arm, metric):
    p = frozen_path(f"data/qiime/diversity/{arm}/alpha_diversities/", (metric,))
    with zipfile.ZipFile(p) as zf:
        name = next(n for n in zf.namelist() if n.endswith("/data/alpha-diversity.tsv"))
        text = io.TextIOWrapper(zf.open(name), encoding="utf-8")
        next(text)
        return {l.split("\t")[0]: float(l.rstrip("\n").split("\t")[1]) for l in text if l.strip()}


def read_dm(arm, keep):
    p = frozen_path(f"data/qiime/diversity/{arm}/distance_matrices/", BC_NAMES)
    with zipfile.ZipFile(p) as zf:
        name = next(n for n in zf.namelist() if n.endswith("/data/distance-matrix.tsv"))
        text = io.TextIOWrapper(zf.open(name), encoding="utf-8")
        header = text.readline().rstrip("\n").split("\t")[1:]
        idx = {s: i for i, s in enumerate(header)}
        dm = {}
        for line in text:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in keep:
                dm[parts[0]] = {s: float(parts[1 + idx[s]]) for s in keep}
    return dm


def read_table_and_families(arm):
    fam = {}
    with open(Path("results/taxonomy") / arm / "export" / "taxonomy.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            f = ""
            for part in r["Taxon"].split(";"):
                part = part.strip()
                if part.startswith("f__"):
                    f = part[3:].strip()
            fam[r["Feature ID"]] = f or "(no family)"
    with open(Path("results/dada2") / arm / "table_export" / "feature-table.tsv") as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("# Constructed")]
    cols = lines[0][1:]
    asvs = [r[0] for r in lines[1:]]
    counts = {s: [float(r[1 + i]) for r in lines[1:]] for i, s in enumerate(cols)}
    return asvs, [fam.get(a, "(no taxonomy)") for a in asvs], counts


# ---------- weighting ----------
def cage_weighted(mice, mouse_value):
    cages = sorted({d["cage"] for d in mice.values()},
                   key=lambda c: (DIETS.index(next(d["diet"] for d in mice.values() if d["cage"] == c)), c))
    cage = {c: mean([mouse_value[m] for m in mice if mice[m]["cage"] == c]) for c in cages}
    n = {c: sum(1 for m in mice if mice[m]["cage"] == c) for c in cages}
    diet_of = {c: next(d["diet"] for d in mice.values() if d["cage"] == c) for c in cages}
    diet = {g: mean([cage[c] for c in cages if diet_of[c] == g]) for g in DIETS}
    return {"cages": cages, "cage": cage, "n": n, "diet_of": diet_of, "diet": diet,
            "delta": diet["PD"] - diet["STD"]}


def separation(r):
    pd = [r["cage"][c] for c in r["cages"] if r["diet_of"][c] == "PD"]
    st = [r["cage"][c] for c in r["cages"] if r["diet_of"][c] == "STD"]
    if max(pd) < min(st):
        return "all PD cages below all STD cages"
    if min(pd) > max(st):
        return "all PD cages above all STD cages"
    return "cages overlap between diets"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mice = load_mice()
    win_ids = {m: [d["days"][float(x)] for x in WIN] for m, d in mice.items()}
    base_id = {m: d["days"][0.0] for m, d in mice.items()}
    rep = ["##### EXPLORATORY (not in the frozen specification; D29) - window samples and decomposition",
           f"{len(mice)} mice; cage-equal weighting as in scripts 18-19", ""]
    tsv = []

    # A. window levels
    rep.append("A. Window levels, days 8-12 (all window samples in run 2). Per cage: mean of mouse means.")
    alpha = {}
    for metric, label, fmt in [("observed_features", "Observed ASVs", ".1f"), ("faith_pd", "Faith's PD", ".2f")]:
        for arm in (PRIMARY, SENS):
            alpha[(metric, arm)] = read_alpha(arm, metric)
            r = cage_weighted(mice, {m: mean([alpha[(metric, arm)][s] for s in win_ids[m]]) for m in mice})
            tag = "200/150" if arm == PRIMARY else "218/166"
            rep.append(f"  {label} [{tag}]: STD {r['diet']['STD']:{fmt}}, PD {r['diet']['PD']:{fmt}}, "
                       f"PD - STD {r['delta']:+{fmt}}; {separation(r)}")
            rep.append("     " + "; ".join(f"{r['diet_of'][c]} {c} (n={r['n'][c]}) {r['cage'][c]:{fmt}}" for c in r["cages"]))
            tsv += [("A_window_level", f"{metric}_{tag}", c, r["diet_of"][c], r["n"][c], r["cage"][c], r["delta"]) for c in r["cages"]]
    asvs, fams, counts = read_table_and_families(PRIMARY)
    tot = {s: sum(v) for s, v in counts.items()}

    def share(s, f):
        return sum(x for x, g in zip(counts[s], fams) if g == f) / tot[s]
    for f in FAMILIES:
        r = cage_weighted(mice, {m: mean([share(s, f) for s in win_ids[m]]) for m in mice})
        rep.append(f"  {f} share [200/150, unrarefied]: STD {r['diet']['STD']:.3f}, PD {r['diet']['PD']:.3f}, "
                   f"PD - STD {r['delta']:+.3f}; {separation(r)}")
        rep.append("     " + "; ".join(f"{r['diet_of'][c]} {c} {r['cage'][c]:.3f}" for c in r["cages"]))
        tsv += [("A_window_level", f"{f}_share", c, r["diet_of"][c], r["n"][c], r["cage"][c], r["delta"]) for c in r["cages"]]
    rep.append("")

    # B. Muribaculaceae detection
    rep.append("B. Muribaculaceae detection per sample (reads > 0 in any Muribaculaceae ASV; unrarefied 200/150)")
    for g in DIETS:
        for label, ids in [("day 0", [base_id[m] for m in mice if mice[m]["diet"] == g]),
                           ("days 8-12", [s for m in mice if mice[m]["diet"] == g for s in win_ids[m]])]:
            reads = [sum(x for x, f in zip(counts[s], fams) if f == "Muribaculaceae") for s in ids]
            det = sum(1 for x in reads if x > 0)
            rep.append(f"  {g} {label}: detected in {det} of {len(ids)} samples; "
                       f"Muribaculaceae reads total {int(sum(reads))}, max per sample {int(max(reads))}; "
                       f"median sample depth {sorted(tot[s] for s in ids)[len(ids) // 2]:.0f}")
    rep.append("")

    # C. family decomposition of S2
    rel = {s: [x / tot[s] for x in v] for s, v in counts.items()}
    fam_names = sorted(set(fams))
    parts = defaultdict(dict)   # family -> mouse -> mean part
    for m in mice:
        b = rel[base_id[m]]
        per_day = []
        for s in win_ids[m]:
            w = rel[s]
            d = defaultdict(float)
            for x, y, f in zip(w, b, fams):
                d[f] += 0.5 * abs(x - y)
            per_day.append(d)
        for f in fam_names:
            parts[f][m] = mean([d[f] for d in per_day])
    total_mouse = {m: sum(parts[f][m] for f in fam_names) for m in mice}
    total = cage_weighted(mice, total_mouse)
    fam_res = {f: cage_weighted(mice, parts[f]) for f in fam_names}
    check = sum(r["delta"] for r in fam_res.values())
    rep.append("C. Family parts of the unrarefied relative-abundance Bray-Curtis from day 0 (S2), days 8-12")
    rep.append(f"  total: STD {total['diet']['STD']:.4f}, PD {total['diet']['PD']:.4f}, delta {total['delta']:+.4f} "
               f"(sum of family deltas {check:+.4f}; should equal the S2 delta of script 18)")
    rep.append(f"  {'family':32s} {'STD part':>9s} {'PD part':>9s} {'PD-STD':>9s}  {'share of STD distance':>22s}")
    order = sorted(fam_names, key=lambda f: -abs(fam_res[f]["delta"]))
    for f in order[:10]:
        r = fam_res[f]
        rep.append(f"  {f[:32]:32s} {r['diet']['STD']:9.4f} {r['diet']['PD']:9.4f} {r['delta']:+9.4f}  "
                   f"{r['diet']['STD'] / total['diet']['STD']:22.1%}")
    rest = [f for f in order[10:]]
    if rest:
        rep.append(f"  {'(all other families: ' + str(len(rest)) + ')':32s} "
                   f"{sum(fam_res[f]['diet']['STD'] for f in rest):9.4f} {sum(fam_res[f]['diet']['PD'] for f in rest):9.4f} "
                   f"{sum(fam_res[f]['delta'] for f in rest):+9.4f}")
    for f in fam_names:
        r = fam_res[f]
        tsv += [("C_family_part", f, c, r["diet_of"][c], r["n"][c], r["cage"][c], r["delta"]) for c in r["cages"]]
    rep.append("")

    # D. between-cage distances among window samples
    all_win = {s: m for m in mice for s in win_ids[m]}
    dm = read_dm(PRIMARY, set(all_win))
    cages = total["cages"]
    rep.append("D. Window samples, rarefied medoid Bray-Curtis (200/150): mean distance between cages")
    rep.append("   (pairs of samples from different mice; diagonal = between mice of the same cage)")
    rep.append("   " + " " * 12 + "".join(f"{c:>9s}" for c in cages))
    cmat = {}
    for c1 in cages:
        row = []
        for c2 in cages:
            ds = [dm[a][b] for a, ma in all_win.items() for b, mb in all_win.items()
                  if ma != mb and mice[ma]["cage"] == c1 and mice[mb]["cage"] == c2]
            cmat[(c1, c2)] = mean(ds)
            row.append(f"{cmat[(c1, c2)]:9.3f}")
        rep.append(f"   {total['diet_of'][c1]:3s} {c1:8s}" + "".join(row))
    within_diet = mean([cmat[(a, b)] for a in cages for b in cages if a != b and total["diet_of"][a] == total["diet_of"][b]])
    between = mean([cmat[(a, b)] for a in cages for b in cages if total["diet_of"][a] != total["diet_of"][b]])
    rep.append(f"   mean of cage-pair means: different cages, same diet {within_diet:.3f}; different diets {between:.3f}")
    rep += ["", "Limits: exploratory; window levels avoid the run difference between baselines but compare",
            "dietary groups whose cohort concurrency and sample preparation are unverified; no p-values."]
    (OUT / "exploratory_window_report.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))
    with open(OUT / "exploratory_cage_values.tsv", "w") as fh:
        fh.write("section\tmeasure\tcage_id\tdiet\tn_mice\tcage_value\tdelta_PD_minus_STD\n")
        for t in tsv:
            fh.write("\t".join(str(x) if not isinstance(x, float) else f"{x:.6f}" for x in t) + "\n")
    pcoa_plot(dm, all_win, mice)
    print(f"\nWrote {OUT}/exploratory_window_report.txt, exploratory_cage_values.tsv, window_pcoa.png")


def pcoa_plot(dm, all_win, mice):
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("numpy/matplotlib not available: figure skipped")
        return
    ids = sorted(all_win)
    D = np.array([[dm[a][b] for b in ids] for a in ids])
    n = len(ids)
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    pos = vals[vals > 0].sum()
    coords = vecs[:, :2] * np.sqrt(np.maximum(vals[:2], 0))
    col = {"STD": "#3E5C76", "PD": "#C0703A"}
    marks = ["o", "s", "^", "D", "v"]
    cages = sorted({mice[m]["cage"] for m in mice})
    fig, ax = plt.subplots(figsize=(6.2, 4.8), dpi=150)
    for k, c in enumerate(cages):
        idx = [i for i, s in enumerate(ids) if mice[all_win[s]]["cage"] == c]
        diet = mice[all_win[ids[idx[0]]]]["diet"]
        ax.scatter(coords[idx, 0], coords[idx, 1], marker=marks[k % len(marks)], s=30,
                   facecolor=col[diet], edgecolor="white", linewidth=0.6, alpha=0.85, label=f"{c} {diet}")
    ax.set_xlabel(f"PCoA 1 ({vals[0] / pos:.0%})")
    ax.set_ylabel(f"PCoA 2 ({vals[1] / pos:.0%})")
    ax.set_title("EXPLORATORY: window samples (days 8-12, all run 2), Bray-Curtis 200/150", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "window_pcoa.png")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        sys.exit(__doc__)
    main()

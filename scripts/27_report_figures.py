#!/usr/bin/env python3
"""27_report_figures.py - figures for docs/STAGE3_REPORT.md (committed in docs/figures/).

Built only from result tables already written by scripts 18, 20 and 25:
  results/recovery/mouse_values.tsv, results/recovery/cage_values.tsv   (script 18)
  results/exploratory/exploratory_cage_values.tsv                        (script 20)
  results/stage2/cage_values.tsv                                         (script 25)
  data/qiime/diversity/trunc200_150/distance_matrices/braycurtis.qza    (script 16; frozen
  rarefied medoid matrix) with metadata/primary_analysis_samples.tsv     (script 17)
Figures:
  fig1_primary.png   A: frozen primary measure, mouse values and cage means;
                     B: PD - STD contrast in every distance-from-baseline analysis,
                        as % of the STD mean (scales differ between metrics)
  fig2_window_levels.png  exploratory window levels (days 8-12, all run 2), per cage
  fig3_window_pcoa.png    PCoA of the 85 window samples (same matrix and method as script 20)
  fig4_stage2.png    cage values, QIIME 2 vs nf-core/ampliseq
Colours: STD #2a78d6 (circles), PD #eb6834 (squares); validated palette (dataviz).
Usage (repository root, QIIME 2 environment active for matplotlib):
  python3 scripts/27_report_figures.py
"""
import csv
import io
import zipfile
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("docs/figures")
COL = {"STD": "#2a78d6", "PD": "#eb6834"}
MARK = {"STD": "o", "PD": "s"}
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "font.family": "DejaVu Sans",   # bundled with matplotlib: identical rendering on every machine
    "text.color": INK, "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
})


def tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def mean(x):
    return sum(x) / len(x)


def cage_order(rows):
    cages = {(r["diet"], r["cage_id"]) for r in rows}
    return sorted(cages, key=lambda t: (0 if t[0] == "STD" else 1, t[1]))


def fig1():
    mv = [r for r in tsv("results/recovery/mouse_values.tsv") if r["analysis"] == "PRIMARY"]
    cv = tsv("results/recovery/cage_values.tsv")
    prim = [r for r in cv if r["analysis"] == "PRIMARY"]
    order = cage_order(prim)
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.0, 4.3), gridspec_kw={"width_ratios": [1.15, 1]})
    for i, (diet, cage) in enumerate(order):
        ys = [float(r["mean_distance_from_day0"]) for r in mv if r["cage_id"] == cage]
        xs = [i + 0.09 * (k - (len(ys) - 1) / 2) for k in range(len(ys))]
        a.scatter(xs, ys, s=26, marker=MARK[diet], facecolor=COL[diet], edgecolor=SURF, linewidth=0.8, alpha=0.75, zorder=3)
        cm = next(float(r["cage_mean"]) for r in prim if r["cage_id"] == cage)
        a.hlines(cm, i - 0.28, i + 0.28, color=COL[diet], lw=2.2, zorder=4)
    for diet in ("STD", "PD"):
        dm = mean([float(r["cage_mean"]) for r in prim if r["diet"] == diet])
        a.axhline(dm, color=COL[diet], lw=1, ls=(0, (4, 3)), zorder=2)
        # label in the right margin, outside the data area (no collision with cages)
        a.annotate(f"{diet} mean\nof cages\n{dm:.3f}", xy=(1.0, dm), xycoords=("axes fraction", "data"),
                   xytext=(4, 0), textcoords="offset points", va="center", ha="left", fontsize=7.5, color=INK2,
                   annotation_clip=False)
    n = {r["cage_id"]: r["n_mice"] for r in prim}
    a.set_xticks(range(len(order)))
    a.set_xticklabels([f"{c.replace('H_', '')}\n{d}\nn={n[c]}" for d, c in order])
    a.set_ylabel("Bray–Curtis from own day 0\n(mean of days 8–12, per mouse)")
    delta = float(prim[0]["delta_PD_minus_STD"])
    a.set_title(f"A  Frozen primary measure: PD − STD = {delta:+.4f}".replace("-", "−"), loc="left")
    a.grid(axis="x", visible=False)
    # B: contrast as % of STD mean, every analysis
    labels = {"PRIMARY": "Bray–Curtis (primary)", "S1": "Bray–Curtis, 218/166", "S2": "Bray–Curtis, unrarefied",
              "S4": "Bray–Curtis, days 8–13", "SEC-wU": "Weighted UniFrac", "S3": "Aitchison"}
    rows = []
    for aid in ["PRIMARY", "S1", "S2", "S4", "SEC-wU", "S3"]:
        rs = [r for r in cv if r["analysis"] == aid]
        if not rs:
            continue
        std = mean([float(r["cage_mean"]) for r in rs if r["diet"] == "STD"])
        pd = mean([float(r["cage_mean"]) for r in rs if r["diet"] == "PD"])
        rows.append((labels[aid], 100 * (pd - std) / std, pd - std))
    y = list(range(len(rows)))[::-1]
    b.axvline(0, color=INK2, lw=0.9)
    for yi, (lab, pct, d) in zip(y, rows):
        b.hlines(yi, min(0, pct), max(0, pct), color=INK2, lw=1.2, alpha=0.5)
        b.scatter([pct], [yi], s=46, color=INK, zorder=3)
        # absolute contrast printed right of the zero line, clear of points and axis
        b.text(0.8, yi, (f"{d:+.4f}" if abs(d) < 1 else f"{d:+.2f}").replace("-", "−"), ha="left", va="center",
               fontsize=7.5, color=INK2)
    b.set_yticks(y)
    b.set_yticklabels([r[0] for r in rows])
    b.set_xlabel("PD − STD, % of the STD mean (cage-equal)")
    b.set_title("B  Contrast in every distance analysis", loc="left")
    b.grid(axis="y", visible=False)
    lo = min(r[1] for r in rows)
    b.set_xlim(min(-20, lo * 1.15), 9)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_primary.png", dpi=200)
    plt.close(fig)


def fig2():
    rows = [r for r in tsv("results/exploratory/exploratory_cage_values.tsv") if r["section"] == "A_window_level"]
    panels = [("observed_features_200/150", "Observed ASVs", 1, "{:.0f}"), ("faith_pd_200/150", "Faith’s PD", 1, "{:.1f}"),
              ("Muribaculaceae_share", "Muribaculaceae (%)", 100, "{:.0f}"),
              ("Bacteroidaceae_share", "Bacteroidaceae (%)", 100, "{:.0f}"),
              ("Akkermansiaceae_share", "Akkermansiaceae (%)", 100, "{:.0f}")]
    fig, axes = plt.subplots(1, 5, figsize=(12.5, 3.1))
    for ax, (meas, title, sc, fmt) in zip(axes, panels):
        rs = [r for r in rows if r["measure"] == meas]
        order = cage_order(rs)
        for i, (diet, cage) in enumerate(order):
            v = sc * float(next(r["cage_value"] for r in rs if r["cage_id"] == cage))
            ax.scatter([i], [v], s=46, marker=MARK[diet], color=COL[diet], zorder=3, clip_on=False)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([c.replace("H_", "") for _, c in order], fontsize=8)
        ax.set_title(title, fontsize=9, loc="left")
        ax.set_xlim(-0.6, len(order) - 0.4)
        ax.grid(axis="x", visible=False)
        top = max(sc * float(r["cage_value"]) for r in rs)
        ax.set_ylim(-0.05 * top if top > 0 else 0, top * 1.18 if top > 0 else 1)
    handles = [plt.Line2D([], [], marker=MARK[d], color=COL[d], ls="", markersize=7, label=d) for d in ("STD", "PD")]
    fig.legend(handles=handles, loc="upper right", ncol=2, frameon=False, fontsize=9)
    fig.suptitle("Exploratory: window levels, days 8–12 (all samples in sequencing run 2); one point per cage, cage-equal",
                 x=0.01, ha="left", fontsize=10, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT / "fig2_window_levels.png", dpi=200)
    plt.close(fig)


def fig4():
    rs = tsv("results/stage2/cage_values.tsv")
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    xs = [float(r["qiime2"]) for r in rs]
    ys = [float(r["ampliseq"]) for r in rs]
    lo, hi = min(xs + ys) - 0.03, max(xs + ys) + 0.03
    ax.plot([lo, hi], [lo, hi], color=INK2, lw=0.9, ls=(0, (4, 3)), zorder=1)
    for r, x, y in zip(rs, xs, ys):
        ax.scatter([x], [y], s=52, marker=MARK[r["diet"]], color=COL[r["diet"]], edgecolor=SURF, linewidth=0.8, zorder=3)
        ax.text(x + 0.006, y - 0.004, r["cage_id"].replace("H_", ""), fontsize=8, color=INK2, va="top")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("QIIME 2 cage value")
    ax.set_ylabel("nf-core/ampliseq cage value")
    def delta(key):
        return mean([float(r[key]) for r in rs if r["diet"] == "PD"]) - mean([float(r[key]) for r in rs if r["diet"] == "STD"])
    ax.set_title(f"Stage 2: PD − STD {delta('qiime2'):+.4f} vs {delta('ampliseq'):+.4f}".replace("-", "−").replace("nf−core", "nf-core"),
                 loc="left", fontsize=9.5)
    handles = [plt.Line2D([], [], marker=MARK[d], color=COL[d], ls="", markersize=7, label=d) for d in ("STD", "PD")]
    handles.append(plt.Line2D([], [], color=INK2, lw=0.9, ls=(0, (4, 3)), label="same value in both"))
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig4_stage2.png", dpi=200)
    plt.close(fig)


def read_dm(path, keep):
    with zipfile.ZipFile(path) as zf:
        name = next(n for n in zf.namelist() if n.endswith("/data/distance-matrix.tsv"))
        text = io.TextIOWrapper(zf.open(name), encoding="utf-8")
        header = text.readline().rstrip("\n").split("\t")[1:]
        idx = {s: i for i, s in enumerate(header)}
        dm = {}
        for line in text:
            p = line.rstrip("\n").split("\t")
            if p[0] in keep:
                dm[p[0]] = {s: float(p[1 + idx[s]]) for s in keep}
    return dm


def fig3():
    import numpy as np
    win = {r["sample-id"]: r for r in tsv("metadata/primary_analysis_samples.tsv") if r["role"] == "window"}
    dm = read_dm(Path("data/qiime/diversity/trunc200_150/distance_matrices/braycurtis.qza"), set(win))
    ids = sorted(win)
    D = np.array([[dm[a][b] for b in ids] for a in ids])
    n = len(ids)
    J = np.eye(n) - np.ones((n, n)) / n
    vals, vecs = np.linalg.eigh(-0.5 * J @ (D ** 2) @ J)
    o = np.argsort(vals)[::-1]
    vals, vecs = vals[o], vecs[:, o]
    xy = vecs[:, :2] * np.sqrt(np.maximum(vals[:2], 0))
    pos = vals[vals > 0].sum()
    order = sorted({(win[s]["diet"], win[s]["cage_id"]) for s in ids}, key=lambda t: (0 if t[0] == "STD" else 1, t[1]))
    shapes = {"STD": ["o", "^", "v"], "PD": ["s", "D"]}
    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    k = {"STD": 0, "PD": 0}
    for diet, cage in order:
        sel = [i for i, s in enumerate(ids) if win[s]["cage_id"] == cage]
        ax.scatter(xy[sel, 0], xy[sel, 1], s=30, marker=shapes[diet][k[diet] % len(shapes[diet])], color=COL[diet],
                   edgecolor=SURF, linewidth=0.6, alpha=0.85, label=f"{cage.replace('H_', '')} {diet}", zorder=3)
        k[diet] += 1
    ax.set_xlabel(f"PCoA 1 ({vals[0] / pos:.0%} of variation)")
    ax.set_ylabel(f"PCoA 2 ({vals[1] / pos:.0%})")
    ax.set_title("Exploratory: window samples, days 8–12 (all run 2)\nBray–Curtis, rarefied medoid, 200/150", loc="left",
                 fontsize=9.5)
    ax.legend(frameon=False, fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(OUT / "fig3_window_pcoa.png", dpi=200)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig1(); fig2(); fig3(); fig4()
    for p in sorted(OUT.glob("*.png")):
        print(f"wrote {p} ({p.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        sys.exit(__doc__)
    main()

#!/usr/bin/env python3
"""28_pcoa_split_check.py - EXPLORATORY: the two PD groups in the window PCoA (D34).

Report figure 3 (script 27) shows the PD window samples (days 8-12, all run 2)
in two groups along PCoA axis 2, each containing both PD cages. This check asks
whether the groups follow sampling day, individual mice, or a processing proxy.
Same ordination as figure 3: 85 window samples of the frozen set, rarefied medoid
Bray-Curtis (200/150), classical PCoA.
  1. PD groups: split at the largest gap in PD axis-2 scores.
  2. Group membership by day, by mouse and by cage; purity = share of samples
     in the majority group within each day / mouse (1.0 = fully determined).
  3. Processing proxies: original sample ID (numbering order), DADA2 final reads,
     demultiplexed read pairs. Sequencing run is constant (run 2).
  4. Composition: families and ASVs differing most between the groups (mean
     relative abundance, unrarefied 200/150).
  5. Figure: A all window samples coloured by day; B PD samples labelled by mouse.
Usage (repository root, QIIME 2 environment active for numpy/matplotlib):
  python3 scripts/28_pcoa_split_check.py
Outputs: results/exploratory/pcoa_split_report.txt, pcoa_split.png
"""
import csv
import io
import statistics
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path("results/exploratory")
ARM = "trunc200_150"
DAYS = [8, 9, 10, 11, 12]
DAYCOL = dict(zip(DAYS, ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]))   # ordinal ramp, validated
MARK = {"STD": "o", "PD": "s"}
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
})


def tsv(path):
    with open(path) as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if not list(r.values())[0].startswith("#")]


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


def table(arm):
    with open(Path("results/dada2") / arm / "table_export" / "feature-table.tsv") as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("# Constructed")]
    cols = lines[0][1:]
    out = {s: {} for s in cols}
    for row in lines[1:]:
        for s, x in zip(cols, row[1:]):
            if float(x):
                out[s][row[0]] = float(x)
    return out


def taxonomy(arm):
    out = {}
    for r in tsv(Path("results/taxonomy") / arm / "export" / "taxonomy.tsv"):
        ranks = {p.strip()[0]: p.strip()[3:] for p in r["Taxon"].split(";") if len(p.strip()) > 3 and p.strip()[1:3] == "__"}
        out[r["Feature ID"]] = (ranks.get("f", "") or "(no family)", ranks.get("g", "") or "(no genus)")
    return out


def purity(groups, key):
    by = defaultdict(Counter)
    for s, g in groups.items():
        by[key(s)][g] += 1
    return sum(c.most_common(1)[0][1] for c in by.values()) / len(groups), by


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    win = {r["sample-id"]: r for r in tsv("metadata/primary_analysis_samples.tsv") if r["role"] == "window"}
    meta = {r["sample-id"]: r for r in tsv("metadata/humanised_sample_metadata.tsv")}
    ids = sorted(win)
    dm = read_dm(Path(f"data/qiime/diversity/{ARM}/distance_matrices/braycurtis.qza"), set(ids))
    D = np.array([[dm[a][b] for b in ids] for a in ids])
    n = len(ids)
    J = np.eye(n) - np.ones((n, n)) / n
    vals, vecs = np.linalg.eigh(-0.5 * J @ (D ** 2) @ J)
    o = np.argsort(vals)[::-1]
    vals, vecs = vals[o], vecs[:, o]
    xy = vecs[:, :2] * np.sqrt(np.maximum(vals[:2], 0))
    pos = {s: (xy[i, 0], xy[i, 1]) for i, s in enumerate(ids)}
    pd_ids = [s for s in ids if win[s]["diet"] == "PD"]
    # 'upper'/'lower' as plotted in report figure 3 (same computation, same axis sign)
    y = sorted((pos[s][1], s) for s in pd_ids)
    gaps = [(y[i + 1][0] - y[i][0], i) for i in range(len(y) - 1)]
    gap, cut = max(gaps)
    threshold = (y[cut][0] + y[cut + 1][0]) / 2
    groups = {s: ("upper" if pos[s][1] > threshold else "lower") for s in pd_ids}
    sizes = Counter(groups.values())
    rep = ["##### EXPLORATORY (D34): the two PD groups in the window PCoA (report figure 3)",
           f"85 window samples; PCoA axis 1 {vals[0] / vals[vals > 0].sum():.0%}, axis 2 {vals[1] / vals[vals > 0].sum():.0%}; "
           f"PD samples split at the largest gap on axis 2 ({gap:.3f}; next largest {sorted(g for g, _ in gaps)[-2]:.3f}): "
           f"upper {sizes['upper']}, lower {sizes['lower']}", ""]

    # 2. membership
    day = lambda s: int(float(win[s]["day"]))
    mouse = lambda s: win[s]["mouse_id"]
    p_day, by_day = purity(groups, day)
    p_mouse, by_mouse = purity(groups, mouse)
    rep.append("2. Group membership (PD window samples)")
    rep.append("   by day:   " + "; ".join(f"day {d}: upper {by_day[d]['upper']}, lower {by_day[d]['lower']}" for d in DAYS))
    rep.append("   by mouse: " + "; ".join(f"{m} ({win[next(s for s in pd_ids if mouse(s) == m)]['cage_id']}): "
                                           f"{by_mouse[m]['upper']}/{by_mouse[m]['lower']}"
                                           for m in sorted(by_mouse, key=lambda x: int(x.split('m')[-1]))) + "  (upper/lower)")
    p_cage, by_cage = purity(groups, lambda s: win[s]["cage_id"])
    rep.append("   by cage:  " + "; ".join(f"{c}: {by_cage[c]['upper']}/{by_cage[c]['lower']}" for c in sorted(by_cage)))
    rep.append(f"   purity (share of samples in the majority group within each level): day {p_day:.2f}, mouse {p_mouse:.2f}, "
               f"cage {p_cage:.2f}  (1.00 = groups fully determined by that variable; "
               f"{max(sizes.values()) / len(pd_ids):.2f} = no information)")
    rep.append("")

    # 3. processing proxies
    tab = table(ARM)
    depth = {s: sum(tab[s].values()) for s in pd_ids}
    demux = {}
    pf = Path("results/demux/per_sample_read_pairs.tsv")
    if pf.exists():
        for r in tsv(pf):
            k = next((c for c in r if c.lower().startswith("sample")), None)
            v = next((c for c in r if "read" in c.lower() or "count" in c.lower()), None)
            if k and v:
                try:
                    demux[r[k]] = float(r[v])
                except ValueError:
                    pass
    rep.append("3. Processing proxies (all PD window samples are in run 2)")
    for label, f in [("original sample ID", lambda s: float(meta[s]["original_sample_id"])),
                     ("DADA2 final reads", lambda s: depth[s]),
                     ("demultiplexed read pairs", lambda s: demux.get(s))]:
        line = []
        for g in ("upper", "lower"):
            v = [f(s) for s in pd_ids if groups[s] == g and f(s) is not None]
            if v:
                line.append(f"{g} median {statistics.median(v):,.0f} (range {min(v):,.0f}-{max(v):,.0f})")
        rep.append(f"   {label:24s} " + ("; ".join(line) if line else "n/a (input not found)"))
    rep.append("   note: original sample IDs may follow collection order, so an ID difference can mirror a day difference")
    rep.append("")

    # 4. composition
    tax = taxonomy(ARM)
    rel = {s: {a: c / sum(tab[s].values()) for a, c in tab[s].items()} for s in pd_ids}

    def gmean(g, key):
        ss = [s for s in pd_ids if groups[s] == g]
        acc = Counter()
        for s in ss:
            for a, r in rel[s].items():
                acc[key(a)] += r
        return {k: v / len(ss) for k, v in acc.items()}
    rep.append("4. Largest differences in mean relative abundance between the PD groups (unrarefied)")
    for level, key in [("family", lambda a: tax.get(a, ("?", "?"))[0]), ("ASV", lambda a: a)]:
        up, lo = gmean("upper", key), gmean("lower", key)
        top = sorted(set(up) | set(lo), key=lambda k: -abs(up.get(k, 0) - lo.get(k, 0)))[:8]
        rep.append(f"   {level}:")
        for k in top:
            name = k if level == "family" else f"{k[:10]}... {tax.get(k, ('?', '?'))[0]} / {tax.get(k, ('?', '?'))[1]}"
            rep.append(f"     {name[:60]:60s} upper {up.get(k, 0):.3f}  lower {lo.get(k, 0):.3f}  "
                       f"difference {up.get(k, 0) - lo.get(k, 0):+.3f}")
    rep += ["", "Interpretation guide: purity near 1.00 for day -> temporal shift; for mouse -> persistent",
            "differences between mice; neither -> look at processing proxies and composition; otherwise unresolved."]
    (OUT / "pcoa_split_report.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))

    # 5. figure
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.0, 4.6))
    for s in ids:
        x_, y_ = pos[s]
        a.scatter([x_], [y_], s=30, marker=MARK[win[s]["diet"]], color=DAYCOL[day(s)], edgecolor=SURF, linewidth=0.6, zorder=3)
    hd = [plt.Line2D([], [], marker="o", color=DAYCOL[d], ls="", markersize=7, label=f"day {d}") for d in DAYS]
    hd += [plt.Line2D([], [], marker=MARK[g], color=INK2, ls="", markersize=7, label=g) for g in ("STD", "PD")]
    a.legend(handles=hd, frameon=False, fontsize=8, loc="best", ncol=2)
    a.set_xlabel(f"PCoA 1 ({vals[0] / vals[vals > 0].sum():.0%})")
    a.set_ylabel(f"PCoA 2 ({vals[1] / vals[vals > 0].sum():.0%})")
    a.set_title("A  Window samples by sampling day", loc="left")
    for s in pd_ids:
        x_, y_ = pos[s]
        b.scatter([x_], [y_], s=34, marker="s", color=DAYCOL[day(s)], edgecolor=SURF, linewidth=0.6, zorder=3)
        b.annotate(f"{mouse(s).replace('H_', '')}·d{day(s)}", (x_, y_), xytext=(4, 2), textcoords="offset points",
                   fontsize=6.5, color=INK2)
    b.axhline(threshold, color=INK2, lw=0.8, ls=(0, (4, 3)))
    b.set_xlabel("PCoA 1")
    b.set_ylabel("PCoA 2")
    b.set_title("B  PD samples, labelled mouse·day", loc="left")
    fig.suptitle("Exploratory (D34): structure of the PD window samples", x=0.01, ha="left", fontsize=10, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "pcoa_split.png", dpi=200)
    print(f"\nWrote {OUT}/pcoa_split_report.txt and pcoa_split.png")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        sys.exit(__doc__)
    main()

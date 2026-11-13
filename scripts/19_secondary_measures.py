#!/usr/bin/env python3
"""19_secondary_measures.py - secondary measures of the frozen specification
(docs/ANALYSIS_SPEC_FROZEN_v1.md section 4, D21; tag analysis-frozen-v1).

Same 17 mice and weighting as script 18 (metadata/primary_analysis_samples.tsv):
per mouse, mean over days 8-12 minus the day-0 value; per cage, mean of mice;
per diet, unweighted mean of cages; delta = mean(PD cages) - mean(STD cages).
Descriptive; no p-value.

  Observed ASVs, Faith's PD : rarefied means (q2-boots, depth 6,400, n = 100),
                              200/150; D21 check: same in 218/166, both deltas
                              reported with their difference.
  Families                  : share of reads per sample, unrarefied 200/150 table,
                              SILVA 144 family (f__) of each ASV: Muribaculaceae,
                              Bacteroidaceae, Akkermansiaceae.
Day-0 and window levels are reported per cage as well, because a change can
come from a different starting level (baselines sit in different runs).
Illustration: cage-mean trajectories over all faecal days (same 17 mice).

Alpha inputs are checked against docs/frozen_artifacts.tsv (MD5).
Outputs: results/recovery/secondary_report.txt, secondary_cage_values.tsv,
secondary_trajectories.png
Usage (repository root, QIIME 2 environment active):
  python3 scripts/19_secondary_measures.py
"""
import csv
import hashlib
import io
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

PRIMARY, SENS = "trunc200_150", "trunc218_166"
OUT = Path("results/recovery")
SAMPLES = Path("metadata/primary_analysis_samples.tsv")
META = Path("metadata/humanised_sample_metadata.tsv")
EXCL = Path("metadata/excluded_samples.tsv")
D23 = {"H3058", "H3072"}
FROZEN = Path("docs/frozen_artifacts.tsv")
WIN = [8, 9, 10, 11, 12]
FAMILIES = ["Muribaculaceae", "Bacteroidaceae", "Akkermansiaceae"]
DIETS = ("STD", "PD")


def mean(x):
    return sum(x) / len(x)


def load_mice():
    mice = {}
    with open(SAMPLES) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m = mice.setdefault(r["mouse_id"], {"cage": r["cage_id"], "diet": r["diet"], "days": {}})
            if r["role"] in ("baseline", "window"):
                m["days"][float(r["day"])] = r["sample-id"]
    return mice


def all_faecal(mice):
    """All faecal samples of the 17 mice (for trajectories only), minus exclusions."""
    with open(EXCL) as fh:
        excl = {r["sample-id"] for r in csv.DictReader(fh, delimiter="\t")} | D23
    out = []
    with open(META) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["sample-id"].startswith("#") or r["mouse_id"] not in mice:
                continue
            if r["source"] == "faeces" and r["sample-id"] not in excl:
                out.append((r["sample-id"], r["mouse_id"], float(r["day"])))
    return out


def frozen_alpha(arm, metric):
    with open(FROZEN) as fh:
        rows = {r["path"]: r["md5"] for r in csv.DictReader(fh, delimiter="\t")}
    p = Path(f"data/qiime/diversity/{arm}/alpha_diversities/{metric}.qza")
    if str(p) not in rows:
        sys.exit(f"STOP: {p} not in {FROZEN}")
    if hashlib.md5(p.read_bytes()).hexdigest() != rows[str(p)]:
        sys.exit(f"STOP: {p} differs from the frozen MD5")
    with zipfile.ZipFile(p) as zf:
        name = next(n for n in zf.namelist() if n.endswith("/data/alpha-diversity.tsv"))
        text = io.TextIOWrapper(zf.open(name), encoding="utf-8")
        next(text)
        return {l.split("\t")[0]: float(l.rstrip("\n").split("\t")[1]) for l in text if l.strip()}


def family_shares(arm):
    with open(Path("results/taxonomy") / arm / "export" / "taxonomy.tsv") as fh:
        fam = {}
        for r in csv.DictReader(fh, delimiter="\t"):
            f = ""
            for part in r["Taxon"].split(";"):
                part = part.strip()
                if part.startswith("f__"):
                    f = part[3:].strip()
            fam[r["Feature ID"]] = f
    with open(Path("results/dada2") / arm / "table_export" / "feature-table.tsv") as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("# Constructed")]
    cols = lines[0][1:]
    tot = [0.0] * len(cols)
    per = {f: [0.0] * len(cols) for f in FAMILIES}
    for row in lines[1:]:
        f = fam.get(row[0], "")
        for i, x in enumerate(row[1:]):
            v = float(x)
            tot[i] += v
            if f in per:
                per[f][i] += v
    return {f: {s: per[f][i] / tot[i] for i, s in enumerate(cols)} for f in FAMILIES}


def estimate(mice, value):
    """value: sample -> number. Change = window mean - day 0, per mouse."""
    mouse, base, win = {}, {}, {}
    for m, d in mice.items():
        base[m] = value[d["days"][0.0]]
        win[m] = mean([value[d["days"][float(x)]] for x in WIN])
        mouse[m] = win[m] - base[m]
    cages = sorted({d["cage"] for d in mice.values()}, key=lambda c: (DIETS.index(next(d["diet"] for d in mice.values() if d["cage"] == c)), c))
    res = {"cages": cages, "cage": {}, "base": {}, "win": {}, "n": {}, "diet_of": {}}
    for c in cages:
        ms = [m for m in mice if mice[m]["cage"] == c]
        res["cage"][c] = mean([mouse[m] for m in ms])
        res["base"][c] = mean([base[m] for m in ms])
        res["win"][c] = mean([win[m] for m in ms])
        res["n"][c] = len(ms)
        res["diet_of"][c] = mice[ms[0]]["diet"]
    res["diet"] = {g: mean([res["cage"][c] for c in cages if res["diet_of"][c] == g]) for g in DIETS}
    res["diet_base"] = {g: mean([res["base"][c] for c in cages if res["diet_of"][c] == g]) for g in DIETS}
    res["diet_win"] = {g: mean([res["win"][c] for c in cages if res["diet_of"][c] == g]) for g in DIETS}
    res["delta"] = res["diet"]["PD"] - res["diet"]["STD"]
    return res


def block(name, r, fmt):
    cs = r["cages"]
    sep = ("all PD cages below all STD cages" if max(r["cage"][c] for c in cs if r["diet_of"][c] == "PD")
           < min(r["cage"][c] for c in cs if r["diet_of"][c] == "STD") else
           "all PD cages above all STD cages" if min(r["cage"][c] for c in cs if r["diet_of"][c] == "PD")
           > max(r["cage"][c] for c in cs if r["diet_of"][c] == "STD") else "cages overlap between diets")
    out = [f"{name}: delta (PD - STD, change from day 0) = {r['delta']:+{fmt}}; {sep}",
           f"  STD: day 0 {r['diet_base']['STD']:{fmt}} -> window {r['diet_win']['STD']:{fmt}} (change {r['diet']['STD']:+{fmt}})",
           f"  PD : day 0 {r['diet_base']['PD']:{fmt}} -> window {r['diet_win']['PD']:{fmt}} (change {r['diet']['PD']:+{fmt}})"]
    for c in cs:
        out.append(f"    {r['diet_of'][c]:3s} {c} (n={r['n'][c]}): day 0 {r['base'][c]:{fmt}}, window {r['win'][c]:{fmt}}, change {r['cage'][c]:+{fmt}}")
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mice = load_mice()
    if any(len(d["days"]) != 6 for d in mice.values()):
        sys.exit("STOP: sample list does not give day 0 + days 8-12 for every mouse")
    rep = ["##### Secondary measures (frozen spec v1, section 4): change from own day 0 to days 8-12 mean",
           f"{len(mice)} mice; cage-equal weighting; delta = mean(PD cages) - mean(STD cages); descriptive", ""]
    rows, results = [], {}
    alpha = {}
    for metric, label, fmt in [("observed_features", "Observed ASVs", ".1f"), ("faith_pd", "Faith's PD", ".2f")]:
        for arm in (PRIMARY, SENS):
            alpha[(metric, arm)] = frozen_alpha(arm, metric)
            results[(metric, arm)] = estimate(mice, alpha[(metric, arm)])
        rep += block(f"{label} [200/150]", results[(metric, PRIMARY)], fmt)
        a, b = results[(metric, PRIMARY)]["delta"], results[(metric, SENS)]["delta"]
        rep.append(f"  D21: delta 200/150 {a:+{fmt}}, 218/166 {b:+{fmt}}, difference {b - a:+{fmt}}, "
                   f"same direction: {'yes' if a * b > 0 else 'NO'}"
                   + (f", ratio 218/166 / 200/150 {b / a:.2f}" if a else ""))
        rep.append("")
    shares = family_shares(PRIMARY)
    for f in FAMILIES:
        results[(f, PRIMARY)] = estimate(mice, shares[f])
        rep += block(f"{f} share of reads [200/150, unrarefied]", results[(f, PRIMARY)], ".3f")
        rep.append("")
    rep += ["Limit: day-0 values of STD come from run 1, all other values from run 2; a run difference in",
            "richness or composition enters STD changes only. Compare the day-0 levels before reading a change."]
    (OUT / "secondary_report.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))

    with open(OUT / "secondary_cage_values.tsv", "w") as fh:
        fh.write("measure\tarm\tcage_id\tdiet\tn_mice\tday0\twindow\tchange\tdelta_PD_minus_STD\n")
        for (meas, arm), r in results.items():
            for c in r["cages"]:
                fh.write(f"{meas}\t{arm}\t{c}\t{r['diet_of'][c]}\t{r['n'][c]}\t{r['base'][c]:.6f}\t"
                         f"{r['win'][c]:.6f}\t{r['cage'][c]:.6f}\t{r['delta']:.6f}\n")
    plot(mice, all_faecal(mice), alpha, shares)
    print(f"\nWrote {OUT}/secondary_report.txt, secondary_cage_values.tsv, secondary_trajectories.png")


def plot(mice, samples, alpha, shares):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available: figure skipped")
        return
    panels = [("Observed ASVs", alpha[("observed_features", PRIMARY)]),
              ("Faith's PD", alpha[("faith_pd", PRIMARY)])] + [(f + " share", shares[f]) for f in FAMILIES]
    col = {"STD": "#3E5C76", "PD": "#C0703A"}
    fig, axes = plt.subplots(len(panels), 1, figsize=(7, 2.1 * len(panels)), dpi=140, sharex=True)
    cages = sorted({d["cage"] for d in mice.values()})
    for ax, (title, val) in zip(axes, panels):
        for c in cages:
            by_day = defaultdict(list)
            for s, m, day in samples:
                if mice[m]["cage"] == c and s in val:
                    by_day[day].append(val[s])
            days = sorted(by_day)
            diet = next(d["diet"] for d in mice.values() if d["cage"] == c)
            ax.plot(days, [mean(by_day[d]) for d in days], color=col[diet], lw=1.2, alpha=0.8, marker="o", ms=2.5,
                    label=f"{c} {diet}")
        ax.axvspan(min(WIN), max(WIN), color="#999999", alpha=0.12, lw=0)
        ax.axvline(0, color="#777777", lw=0.8, ls=":")
        ax.set_ylabel(title, fontsize=8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].legend(frameon=False, fontsize=7, ncol=5, loc="upper center", bbox_to_anchor=(0.5, 1.35))
    axes[-1].set_xlabel("day (shaded: days 8-12; dotted: day 0)")
    fig.suptitle("Cage means over time (illustration; 200/150)", fontsize=10, y=1.0)
    fig.tight_layout()
    fig.savefig(OUT / "secondary_trajectories.png", bbox_inches="tight")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        sys.exit(__doc__)
    main()

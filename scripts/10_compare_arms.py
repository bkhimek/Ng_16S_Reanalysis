#!/usr/bin/env python3
"""10_compare_arms.py - technical comparison of two DADA2 arms, per sample (D22).

(a) Retention by DADA2 step, per sample. For each step the pass rate is
    filtered/input, denoised/filtered, merged/denoised, non-chimeric/merged.
    Per run: median pass rate in each arm and the difference, so the step that
    explains a retention gain is visible; and the spread of the per-sample gain
    in final reads (is it broad, or concentrated in a few samples?).
(b) Same sample, two arms: Bray-Curtis dissimilarity between the sample's
    relative-abundance profile in arm A and in arm B, with ASVs aligned by
    exact sequence. 0 = identical profiles. Summarised by run and diet; the 10
    largest are listed. As a technical yardstick, the Bray-Curtis between the
    one repeated mouse-day pair (two libraries from the same mouse and day) is
    shown within arm A.

No taxonomy, diversity or diet contrast is computed here.

Usage (repository root):
  python3 scripts/10_compare_arms.py trunc200_150 trunc218_166
Writes results/dada2/compare_arms_<A>_vs_<B>/ (report + two TSVs).
"""
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

RUNS = ["hum_run1", "hum_run2", "hum_run3"]
STEPS = [("filter", "input", "filtered"), ("denoise", "filtered", "denoised"),
         ("merge", "denoised", "merged"), ("chimera", "merged", "non-chimeric")]
META = "metadata/humanised_sample_metadata.tsv"


def read_tsv(path):
    with open(path) as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t")
                if not list(r.values())[0].startswith("#q2")]


def stats(arm):
    out = {}
    for run in RUNS:
        for r in read_tsv(next((Path("results/dada2") / arm / f"{run}_stats").glob("*.tsv"))):
            out[r["sample-id"]] = {k: float(r[k]) for k in ["input", "filtered", "denoised", "merged", "non-chimeric"]}
    return out


def profiles(arm):
    """sample -> {sequence: relative abundance}"""
    res = Path("results/dada2") / arm
    seqs, name = {}, None
    for line in open(next((res / "repseqs_export").glob("*.fasta"))):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = ""
        elif name:
            seqs[name] += line
    with open(res / "table_export" / "feature-table.tsv") as fh:
        rows = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("# Constructed")]
    samples = rows[0][1:]
    prof = {s: defaultdict(float) for s in samples}
    for row in rows[1:]:
        sq = seqs[row[0]]
        for s, v in zip(samples, row[1:]):
            if float(v) > 0:
                prof[s][sq] += float(v)
    for s in samples:
        tot = sum(prof[s].values())
        prof[s] = {k: v / tot for k, v in prof[s].items()} if tot else {}
    return prof


def bray_curtis(p, q):
    keys = set(p) | set(q)
    return 0.5 * sum(abs(p.get(k, 0.0) - q.get(k, 0.0)) for k in keys)


def q(v, f):
    v = sorted(v)
    return v[min(len(v) - 1, int(round(f * (len(v) - 1))))]


def main(a, b):
    meta = {r["sample-id"]: r for r in read_tsv(META)}
    out = Path("results/dada2") / f"compare_arms_{a}_vs_{b}"
    out.mkdir(parents=True, exist_ok=True)
    sa, sb = stats(a), stats(b)
    common = sorted(set(sa) & set(sb))
    rep = [f"##### Arm comparison {a} (A) vs {b} (B); samples in both: {len(common)}", "",
           "(a) Median per-sample pass rate by DADA2 step (A / B / A-B)"]
    rows_a = ["sample-id\trun\tdiet\t" + "\t".join(f"{s}_{x}" for s, _, _ in STEPS for x in (a, b)) + f"\tfinal_ratio_{a}_over_{b}"]
    by_run = defaultdict(list)
    for s in common:
        by_run[meta[s]["run"]].append(s)
        line = [s, meta[s]["run"], meta[s]["diet"]]
        for _, x, y in STEPS:
            for arm_stats in (sa, sb):
                d = arm_stats[s]
                line.append(f"{d[y]/d[x]:.4f}" if d[x] else "NA")
        line.append(f"{sa[s]['non-chimeric']/sb[s]['non-chimeric']:.4f}" if sb[s]["non-chimeric"] else "NA")
        rows_a.append("\t".join(line))
    for run in RUNS:
        ss = by_run[run]
        parts = []
        for name, x, y in STEPS:
            ra = statistics.median(sa[s][y] / sa[s][x] for s in ss if sa[s][x])
            rb = statistics.median(sb[s][y] / sb[s][x] for s in ss if sb[s][x])
            parts.append(f"{name} {ra:.3f}/{rb:.3f}/{ra-rb:+.3f}")
        rep.append(f"  {run} (n={len(ss)}): " + "; ".join(parts))
        g = [sa[s]["non-chimeric"] / sb[s]["non-chimeric"] for s in ss if sb[s]["non-chimeric"]]
        rep.append(f"      final reads A/B per sample: 5th {q(g,.05):.3f}, median {statistics.median(g):.3f}, "
                   f"95th {q(g,.95):.3f}, min {min(g):.3f}, max {max(g):.3f}; samples with A/B < 1: {sum(x < 1 for x in g)}")
    (out / "retention_by_step_per_sample.tsv").write_text("\n".join(rows_a) + "\n")

    pa, pb = profiles(a), profiles(b)
    bc = {s: bray_curtis(pa[s], pb[s]) for s in common if s in pa and s in pb}
    rep += ["", "(b) Bray-Curtis between the two arms' profiles of the same sample (0 = identical)"]
    grp = defaultdict(list)
    for s, v in bc.items():
        grp[(meta[s]["run"], meta[s]["diet"])].append(v)
    for k in sorted(grp):
        v = grp[k]
        rep.append(f"  {k[0]} {k[1]:3s}: median {statistics.median(v):.4f}, 95th {q(v,.95):.4f}, max {max(v):.4f} (n={len(v)})")
    allv = list(bc.values())
    rep.append(f"  all samples: median {statistics.median(allv):.4f}, 95th {q(allv,.95):.4f}, max {max(allv):.4f}")
    rep.append("  10 largest: " + ", ".join(f"{s}({v:.4f}, {meta[s]['run']}, {meta[s]['diet']}, day {meta[s]['day']})"
                                          for s, v in sorted(bc.items(), key=lambda kv: -kv[1])[:10]))
    rep_pairs = defaultdict(list)
    for s, m in meta.items():
        if m.get("repeat_mouse_day") == "yes" and s in pa:
            rep_pairs[(m["mouse_id"], m["day"], m["source"])].append(s)
    for k, ids in rep_pairs.items():
        if len(ids) == 2:
            rep.append(f"  yardstick, repeated library pair {ids[0]}/{ids[1]} ({k[0]}, day {k[1]}) within arm {a}: "
                       f"{bray_curtis(pa[ids[0]], pa[ids[1]]):.4f}")
    (out / "bray_curtis_between_arms_per_sample.tsv").write_text(
        "sample-id\trun\tdiet\tday\tbray_curtis\n" +
        "".join(f"{s}\t{meta[s]['run']}\t{meta[s]['diet']}\t{meta[s]['day']}\t{v:.5f}\n" for s, v in sorted(bc.items())))
    (out / "report.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))
    print(f"\nWritten: {out}/report.txt and two per-sample TSVs")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])

#!/usr/bin/env python3
"""09_dada2_check.py - technical summary of one or two DADA2 arms.

Judges the denoising on TECHNICAL criteria only (no diversity or diet results):
  1. Read retention per run through each DADA2 stage
     (input -> filtered -> denoised -> merged -> non-chimeric), totals and
     the per-sample median, plus the samples with the lowest final counts.
  2. ASVs: number, length distribution (V4 between the primers should be
     ~253 bp; other lengths flag off-target products), and share of reads in
     ASVs of unexpected length.
  3. Final reads per sample by run, diet and source.
  4. With two arms: ASVs shared by exact sequence, share of reads in shared
     ASVs (each arm), and the per-sample ratio of final counts.

Usage (repository root):
  python3 scripts/09_dada2_check.py trunc218_166 [trunc200_150]
Writes results/dada2/<arm>/dada2_check_report.txt (and a comparison file).
"""
import csv
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

RUNS = ["hum_run1", "hum_run2", "hum_run3"]
STAGES = ["input", "filtered", "denoised", "merged", "non-chimeric"]
META = "metadata/humanised_sample_metadata.tsv"


def read_tsv(path):
    with open(path) as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t")
                if not list(r.values())[0].startswith("#q2")]


def load_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif name:
            seqs[name] += line
    return seqs


def load_table(path):
    """biom TSV: '# Constructed from biom file', then '#OTU ID<TAB>samples...'."""
    with open(path) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("# Constructed")]
    samples = lines[0][1:]
    tab = {}
    for row in lines[1:]:
        tab[row[0]] = dict(zip(samples, (float(x) for x in row[1:])))
    return samples, tab


def arm_summary(arm, meta):
    res = Path("results/dada2") / arm
    rep = [f"##### DADA2 arm {arm}"]
    settings = dict(l.rstrip("\n").split("\t", 1) for l in open(res / "arm_settings.tsv"))
    rep.append("settings: " + ", ".join(f"{k}={v}" for k, v in settings.items()))
    rep.append("")
    rep.append("1. Retention by stage (sum over samples; share of input)")
    for run in RUNS:
        rows = read_tsv(next((res / f"{run}_stats").glob("*.tsv")))
        tot = {s: sum(float(r[s]) for r in rows) for s in STAGES}
        rep.append(f"  {run}: " + " -> ".join(f"{s} {int(tot[s]):,} ({tot[s]/tot['input']:.3f})" for s in STAGES))
        per = [float(r["non-chimeric"]) / float(r["input"]) for r in rows if float(r["input"]) > 0]
        mer = [float(r["merged"]) / float(r["denoised"]) for r in rows if float(r["denoised"]) > 0]
        rep.append(f"      per sample: final/input median {statistics.median(per):.3f} "
                   f"(min {min(per):.3f}); merged/denoised median {statistics.median(mer):.3f} (min {min(mer):.3f})")
    rep.append("")

    seqs = load_fasta(next((res / "repseqs_export").glob("*.fasta")))
    samples, tab = load_table(res / "table_export" / "feature-table.tsv")
    reads_by_asv = {a: sum(v.values()) for a, v in tab.items()}
    total = sum(reads_by_asv.values())
    lens = Counter(len(seqs[a]) for a in tab)
    reads_len = Counter()
    for a, n in reads_by_asv.items():
        reads_len[len(seqs[a])] += n
    rep.append(f"2. ASVs: {len(tab):,}; reads in table: {int(total):,}")
    rep.append("  ASV length (bp): count of ASVs / share of reads, most common first:")
    for L, c in lens.most_common(8):
        rep.append(f"    {L}: {c:,} ASVs / {reads_len[L]/total:.4f} of reads")
    odd = sum(n for L, n in reads_len.items() if not 250 <= L <= 256)
    rep.append(f"  reads in ASVs outside 250-256 bp: {odd/total:.4f}")
    rep.append("")

    per_sample = {s: sum(tab[a].get(s, 0) for a in tab) for s in samples}
    grp = defaultdict(list)
    for s, n in per_sample.items():
        m = meta.get(s, {})
        grp[(m.get("run", "?"), m.get("diet", "?"), m.get("source", "?"))].append(n)
    rep.append("3. Final reads per sample (min / median / max, n)")
    for k in sorted(grp):
        v = grp[k]
        rep.append(f"  {k[0]} {k[1]:3s} {k[2]:6s}: {int(min(v)):>7,} / {int(statistics.median(v)):>7,} / {int(max(v)):>7,}  (n={len(v)})")
    low = sorted((n, s) for s, n in per_sample.items())[:10]
    rep.append("  10 lowest: " + ", ".join(f"{s}({int(n)})" for n, s in low))
    rep.append(f"  samples in table: {len(samples)}")
    return rep, seqs, tab, per_sample


def main(arms):
    meta = {r["sample-id"]: r for r in read_tsv(META)}
    out = {}
    for arm in arms:
        rep, seqs, tab, per_sample = arm_summary(arm, meta)
        out[arm] = (seqs, tab, per_sample)
        (Path("results/dada2") / arm / "dada2_check_report.txt").write_text("\n".join(rep) + "\n")
        print("\n".join(rep) + "\n")

    if len(arms) == 2:
        a, b = arms
        sa = {out[a][0][x]: sum(out[a][1][x].values()) for x in out[a][1]}
        sb = {out[b][0][x]: sum(out[b][1][x].values()) for x in out[b][1]}
        shared = set(sa) & set(sb)
        rep = [f"##### Comparison {a} vs {b} (ASVs matched by exact sequence)",
               f"  ASVs: {a} {len(sa):,}, {b} {len(sb):,}, shared {len(shared):,}",
               f"  share of reads in shared ASVs: {a} {sum(sa[s] for s in shared)/sum(sa.values()):.4f}, "
               f"{b} {sum(sb[s] for s in shared)/sum(sb.values()):.4f}"]
        pa, pb = out[a][2], out[b][2]
        ratios = sorted(pb[s] / pa[s] for s in pa if s in pb and pa[s] > 0)
        rep.append(f"  per-sample final reads {b}/{a}: median {statistics.median(ratios):.3f}, "
                   f"range {ratios[0]:.3f}-{ratios[-1]:.3f} (n={len(ratios)})")
        (Path("results/dada2") / f"compare_{a}_vs_{b}.txt").write_text("\n".join(rep) + "\n")
        print("\n".join(rep))


if __name__ == "__main__":
    if not 1 <= len(sys.argv[1:]) <= 2:
        sys.exit(__doc__)
    main(sys.argv[1:])

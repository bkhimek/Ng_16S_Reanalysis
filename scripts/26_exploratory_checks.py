#!/usr/bin/env python3
"""26_exploratory_checks.py - two EXPLORATORY checks (outside the frozen specs).

A. Muribaculaceae traces in the 35 eligible PD window samples (D29). For every
   Muribaculaceae ASV with reads there (QIIME 2, 200/150, unrarefied):
     - its reads and samples in the PD window;
     - its rank among Muribaculaceae ASVs, and share of Muribaculaceae reads,
       in STD samples of the same run (run 2);
     - whether PD mice carried it before the antibiotic (faecal days <= 0, all
       PD mice; day -14 is run 1, days -3..0 run 2).
   Cross-talk yardstick: abundant run-2 STD ASVs (top 30 by reads) with no read
   in any pre-antibiotic PD sample. Their reads in PD window samples estimate
   leakage between samples of the same run (or other transfer). A trace ASV
   that behaves like these fits cross-talk; one carried by PD mice before the
   antibiotic also fits survival. Neither is proof.
B. Samples where QIIME 2 and ampliseq differ most (same-sample Bray-Curtis >
   0.01, results of script 25): the ASVs making up the difference, and for
   workflow-specific ASVs the nearest ASV in the other workflow (Hamming
   distance, same length) and whether that neighbour gained the reads, i.e.
   reads assigned to a different ASV rather than lost.
Usage (repository root): python3 scripts/26_exploratory_checks.py
Output: results/exploratory/checks_report.txt
"""
import csv
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

OUT = Path("results/exploratory")
S1 = "trunc200_150"
WIN = [8, 9, 10, 11, 12]
D23 = {"H3058", "H3072"}


def fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = ""
        elif name is not None:
            seqs[name] += line.upper()
    return seqs


def table_tsv(path):
    with open(path) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("# Constructed")]
    cols = lines[0][1:]
    out = {s: {} for s in cols}
    for row in lines[1:]:
        for s, x in zip(cols, row[1:]):
            if float(x):
                out[s][row[0]] = float(x)
    return out


def families(path):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            f = ""
            for part in r["Taxon"].split(";"):
                part = part.strip()
                if part.startswith("f__"):
                    f = part[3:].strip()
            out[r["Feature ID"]] = f or "(no family)"
    return out


def by_seq(tab, seqs):
    return {s: {seqs[f]: v for f, v in d.items()} for s, d in tab.items()}


def ham(a, b):
    return sum(x != y for x, y in zip(a, b)) if len(a) == len(b) else None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    meta = {}
    with open("metadata/humanised_sample_metadata.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r["sample-id"].startswith("#"):
                meta[r["sample-id"]] = r
    with open("metadata/excluded_samples.tsv") as fh:
        excl = {r["sample-id"] for r in csv.DictReader(fh, delimiter="\t")} | D23
    prim = {}
    with open("metadata/primary_analysis_samples.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            prim[r["sample-id"]] = r
    s_seqs = fasta(next((Path("results/dada2") / S1 / "repseqs_export").glob("*.fasta")))
    S = by_seq(table_tsv(Path("results/dada2") / S1 / "table_export" / "feature-table.tsv"), s_seqs)
    fam = {s_seqs[i]: f for i, f in families(Path("results/taxonomy") / S1 / "export" / "taxonomy.tsv").items() if i in s_seqs}
    rel = {s: {x: n / sum(d.values()) for x, n in d.items()} for s, d in S.items() if d}

    def faecal(s):
        return s in rel and s not in excl and meta[s]["source"] == "faeces"
    pd_win = [s for s, r in prim.items() if r["diet"] == "PD" and r["role"] == "window"]
    std_r2 = [s for s in rel if faecal(s) and meta[s]["diet"] == "STD" and meta[s]["run"] == "hum_run2"]
    pd_pre = [s for s in rel if faecal(s) and meta[s]["diet"] == "PD" and float(meta[s]["day"]) <= 0]
    pd_win_all = [s for s in rel if faecal(s) and meta[s]["diet"] == "PD" and float(meta[s]["day"]) in WIN]

    rep = ["##### EXPLORATORY checks (outside the frozen specifications)", ""]
    rep.append(f"A. Muribaculaceae traces: {len(pd_win)} eligible PD window samples; {len(std_r2)} STD faecal samples in run 2; "
               f"{len(pd_pre)} PD faecal samples before the antibiotic (days <= 0; runs "
               f"{', '.join(sorted(Counter(meta[s]['run'] for s in pd_pre)))})")
    mur_std = Counter()
    for s in std_r2:
        for x, n in S[s].items():
            if fam.get(x) == "Muribaculaceae":
                mur_std[x] += n
    mur_rank = {x: i + 1 for i, (x, _) in enumerate(mur_std.most_common())}
    mur_tot = sum(mur_std.values())
    trace = Counter()
    for s in pd_win:
        for x, n in S[s].items():
            if fam.get(x) == "Muribaculaceae":
                trace[x] += n
    for x, n in trace.most_common():
        nwin = sum(1 for s in pd_win if S[s].get(x))
        pre = [s for s in pd_pre if S[s].get(x)]
        pre_mice = sorted({meta[s]["mouse_id"] for s in pre})
        rep.append(f"   ASV {len(x)} bp ({x[:12]}...): {int(n)} reads in {nwin} of {len(pd_win)} PD window samples; "
                   f"in STD run 2: rank {mur_rank.get(x, 'absent')} of {len(mur_std)} Muribaculaceae ASVs, "
                   f"{mur_std[x] / mur_tot:.3f} of their reads, in {sum(1 for s in std_r2 if S[s].get(x))} of {len(std_r2)} samples, "
                   f"mean rel. abundance {statistics.mean(rel[s].get(x, 0) for s in std_r2):.4f}; "
                   f"PD before antibiotic: {len(pre)} of {len(pd_pre)} samples"
                   + (f" ({len(pre_mice)} mice; days {sorted({meta[s]['day'] for s in pre}, key=float)}; "
                      f"mean rel. abundance where present {statistics.mean(rel[s][x] for s in pre):.4f})" if pre else ""))
    # yardstick: abundant STD run-2 ASVs never seen in PD before the antibiotic
    std_reads = Counter()
    for s in std_r2:
        std_reads.update(S[s])
    absent_pd = [x for x, _ in std_reads.most_common() if not any(S[s].get(x) for s in pd_pre)][:30]
    rows = []
    for x in absent_pd:
        mean_std = statistics.mean(rel[s].get(x, 0) for s in std_r2)
        det = sum(1 for s in pd_win if S[s].get(x))
        reads = sum(S[s].get(x, 0) for s in pd_win)
        rows.append((x, mean_std, det, reads, fam.get(x, "?")))
    rep.append(f"   yardstick: the {len(rows)} most abundant STD run-2 ASVs with no read in any pre-antibiotic PD sample")
    rep.append(f"   {'family':26s} {'STD run 2 mean rel.':>19s} {'PD window samples with reads':>29s} {'reads':>6s} {'reads per 1% STD abundance':>27s}")
    for x, m, det, reads, f in rows[:12]:
        rep.append(f"   {f[:26]:26s} {m:19.4f} {det:>21d} of {len(pd_win):<5d} {int(reads):6d} {reads / (100 * m) if m else 0:27.1f}")
    if rows:
        det_frac = [r[2] / len(pd_win) for r in rows]
        per1 = [r[3] / (100 * r[1]) for r in rows if r[1] > 0]
        rep.append(f"   over all {len(rows)}: median share of PD window samples with reads {statistics.median(det_frac):.2f} "
                   f"(range {min(det_frac):.2f}-{max(det_frac):.2f}); median reads per 1% mean STD abundance {statistics.median(per1):.1f}")
    for x, n in trace.most_common(1):
        m = statistics.mean(rel[s].get(x, 0) for s in std_r2)
        rep.append(f"   main trace ASV for comparison: {sum(1 for s in pd_win if S[s].get(x)) / len(pd_win):.2f} of PD window samples, "
                   f"{n / (100 * m) if m else float('nan'):.1f} reads per 1% mean STD abundance")
    rep.append("   also checked: the main trace ASV in all PD faecal window samples (not only eligible mice): "
               + (f"{sum(1 for s in pd_win_all if S[s].get(trace.most_common(1)[0][0]))} of {len(pd_win_all)}" if trace else "n/a"))
    rep.append("")

    # ---- B ----
    a_seqs = fasta(Path("results/ampliseq/dada2/ASV_seqs.fasta"))
    A = by_seq(table_tsv(Path("results/ampliseq/dada2/ASV_table.tsv")), a_seqs)
    a_fam = {a_seqs[i]: f for i, f in families(Path("results/ampliseq_q2/taxonomy_export/taxonomy.tsv")).items() if i in a_seqs}
    sset, aset = set(s_seqs.values()), set(a_seqs.values())
    per = {}
    with open("results/stage2/per_sample.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            per[r["sample"]] = float(r["bray_curtis_workflows"])
    top = sorted((s for s in per if per[s] > 0.01), key=lambda s: -per[s])
    rep.append(f"B. Samples with QIIME 2 vs ampliseq Bray-Curtis > 0.01: {len(top)}")
    for s in top:
        r = meta[s]
        in_prim = f"in the frozen set ({prim[s]['role']}, {prim[s]['mouse_id']}, {prim[s]['cage_id']})" if s in prim else "not in the frozen set"
        rep.append(f"   {s}: {r['run']}, {r['diet']}, mouse {r['mouse_id']}, day {r['day']}; Bray-Curtis {per[s]:.4f}; {in_prim}")
        ts, ta = sum(S[s].values()), sum(A[s].values())
        contrib = []
        for x in set(S[s]) | set(A[s]):
            d = 0.5 * abs(S[s].get(x, 0) / ts - A[s].get(x, 0) / ta)
            contrib.append((d, x))
        contrib.sort(reverse=True)
        for d, x in contrib[:3]:
            where = "both" if x in sset and x in aset else "QIIME 2 only" if x in sset else "ampliseq only"
            f = fam.get(x) or a_fam.get(x, "?")
            line = (f"     {d:.4f} of the distance: {len(x)} bp, {f}, {where}; reads QIIME 2 {int(S[s].get(x, 0))}, "
                    f"ampliseq {int(A[s].get(x, 0))}")
            if where != "both":
                other = aset if where == "QIIME 2 only" else sset
                near = min(((ham(x, y), y) for y in other if ham(x, y) is not None), default=(None, None))
                if near[0] is not None:
                    y = near[1]
                    line += (f"; nearest ASV in the other workflow: {near[0]} base(s) different, "
                             f"reads QIIME 2 {int(S[s].get(y, 0))} / ampliseq {int(A[s].get(y, 0))}")
                total_here = sum(t.get(x, 0) for t in (S if where == "QIIME 2 only" else A).values())
                line += f"; this ASV's reads over all samples: {int(total_here)}"
            rep.append(line)
    rep += ["", "Limits: exploratory; no blanks or mock community, so cross-talk cannot be measured directly."]
    (OUT / "checks_report.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))
    print(f"\nWrote {OUT}/checks_report.txt")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        sys.exit(__doc__)
    main()

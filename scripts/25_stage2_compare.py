#!/usr/bin/env python3
"""25_stage2_compare.py - Stage 2 comparison, QIIME 2 (Stage 1, 200/150) vs
nf-core/ampliseq 3.1.0, as fixed in docs/STAGE2_COMPARISON_SPEC_v1.md
(tag stage2-spec-v1; D30).

  1. Read retention at equivalent stages, per sample and run; dropped samples;
     the frozen 17 mice / 102 samples complete?
  2. Exact ASV agreement by full sequence (length distributions, orientation
     check, shared ASVs and the share of reads they carry).
  3. Same-sample Bray-Curtis between workflows (unrarefied relative
     abundances, union of ASVs); Stage 1 reference recomputed the same way
     (200/150 vs 218/166).
  4. Family composition with the same SILVA 144 classifier; Muribaculaceae
     traces in the 35 PD window samples in both workflows.
  5. Frozen biological analysis on the ampliseq outputs (script 24): primary
     contrast and cage values, window levels; headline rule:
     same direction and |delta(ampliseq) - (-0.1454)| <= 0.02; pattern check:
     both PD cages below all three STD cages.
Outputs: results/stage2/ (report and tables).
Usage (repository root, after script 24):
  python3 scripts/25_stage2_compare.py
"""
import csv
import io
import math
import statistics
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

FROZEN_DELTA, TOL = -0.1454, 0.02
S1, S1B = "trunc200_150", "trunc218_166"
RUNS = ["hum_run1", "hum_run2", "hum_run3"]
FAMS = ["Muribaculaceae", "Bacteroidaceae", "Akkermansiaceae"]
DIETS = ("STD", "PD")
WIN = [8, 9, 10, 11, 12]
OUT = Path("results/stage2")
AMP = Path("results/ampliseq/dada2")
AQ = Path("data/qiime/ampliseq/diversity")
BC_NAMES = ("braycurtis", "bray_curtis", "bray_curtis_distance_matrix")


def mean(x):
    return sum(x) / len(x)


def q(v, f):
    v = sorted(v)
    return v[min(len(v) - 1, int(round(f * (len(v) - 1))))]


# ---------------- loading ----------------
def fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = ""
        elif name is not None:
            seqs[name] += line.upper()
    return seqs


def table_tsv(path, id_col_names=("#OTU ID", "ASV_ID")):
    """feature table TSV -> {sample: {feature_id: count}}"""
    with open(path) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("# Constructed")]
    if lines[0][0] not in id_col_names:
        sys.exit(f"STOP: unexpected header in {path}: {lines[0][:2]}")
    cols = lines[0][1:]
    out = {s: {} for s in cols}
    for row in lines[1:]:
        for s, x in zip(cols, row[1:]):
            v = float(x)
            if v:
                out[s][row[0]] = v
    return out


def by_sequence(tab, seqs):
    return {s: {seqs[f]: v for f, v in d.items()} for s, d in tab.items()}


def stage1(arm):
    seqs = fasta(next((Path("results/dada2") / arm / "repseqs_export").glob("*.fasta")))
    tab = table_tsv(Path("results/dada2") / arm / "table_export" / "feature-table.tsv")
    return tab, seqs


def families(path):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            f = ""
            for part in r["Taxon"].split(";"):
                part = part.strip()
                if part.startswith("f__"):
                    f = part[3:].strip()
            out[r["Feature ID"]] = (f or "(no family)", r["Taxon"])
    return out


def read_dm(path, keep):
    with zipfile.ZipFile(path) as zf:
        name = next(n for n in zf.namelist() if n.endswith("/data/distance-matrix.tsv"))
        text = io.TextIOWrapper(zf.open(name), encoding="utf-8")
        header = text.readline().rstrip("\n").split("\t")[1:]
        idx = {s: i for i, s in enumerate(header)}
        keep = keep & set(idx)
        dm = {}
        for line in text:
            p = line.rstrip("\n").split("\t")
            if p[0] in keep:
                dm[p[0]] = {s: float(p[1 + idx[s]]) for s in keep}
    return dm


def read_alpha(path):
    with zipfile.ZipFile(path) as zf:
        name = next(n for n in zf.namelist() if n.endswith("/data/alpha-diversity.tsv"))
        text = io.TextIOWrapper(zf.open(name), encoding="utf-8")
        next(text)
        return {l.split("\t")[0]: float(l.rstrip("\n").split("\t")[1]) for l in text if l.strip()}


def one(folder, stems):
    c = [p for p in Path(folder).glob("*.qza") if p.stem in stems]
    if len(c) != 1:
        sys.exit(f"STOP: expected one of {stems} in {folder}, found {c}")
    return c[0]


def bc(p, q_):
    tp, tq = sum(p.values()), sum(q_.values())
    keys = set(p) | set(q_)
    return 0.5 * sum(abs(p.get(k, 0) / tp - q_.get(k, 0) / tq) for k in keys)


def summary(v):
    return f"median {statistics.median(v):.4f}, 95th pct {q(v, 0.95):.4f}, max {max(v):.4f} (n={len(v)})"


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


# ---------------- frozen-estimate machinery (as scripts 18-20) ----------------
def load_mice():
    mice = {}
    with open("metadata/primary_analysis_samples.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m = mice.setdefault(r["mouse_id"], {"cage": r["cage_id"], "diet": r["diet"], "days": {}})
            if r["role"] in ("baseline", "window"):
                m["days"][float(r["day"])] = r["sample-id"]
    return mice


def cage_weighted(mice, value):
    cages = sorted({d["cage"] for d in mice.values()},
                   key=lambda c: (DIETS.index(next(d["diet"] for d in mice.values() if d["cage"] == c)), c))
    cage = {c: mean([value[m] for m in mice if mice[m]["cage"] == c]) for c in cages}
    n = {c: sum(1 for m in mice if mice[m]["cage"] == c) for c in cages}
    diet_of = {c: next(d["diet"] for d in mice.values() if d["cage"] == c) for c in cages}
    diet = {g: mean([cage[c] for c in cages if diet_of[c] == g]) for g in DIETS}
    return {"cages": cages, "cage": cage, "n": n, "diet_of": diet_of, "diet": diet, "delta": diet["PD"] - diet["STD"]}


def pd_below_std(r):
    return max(r["cage"][c] for c in r["cages"] if r["diet_of"][c] == "PD") < \
        min(r["cage"][c] for c in r["cages"] if r["diet_of"][c] == "STD")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rep = ["##### Stage 2 comparison: QIIME 2 (Stage 1, 200/150) vs nf-core/ampliseq 3.1.0",
           "Rules: docs/STAGE2_COMPARISON_SPEC_v1.md (tag stage2-spec-v1)", ""]
    meta = {}
    with open("metadata/humanised_sample_metadata.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r["sample-id"].startswith("#"):
                meta[r["sample-id"]] = r
    mice = load_mice()
    frozen_ids = {s for m in mice.values() for s in m["days"].values()}

    # ---- 1. read retention ----
    s1_stats = {}
    for run in RUNS:
        with open(Path("results/dada2") / S1 / f"{run}_stats" / "stats.tsv") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if not r["sample-id"].startswith("#"):
                    s1_stats[r["sample-id"]] = {"run": run, "input": float(r["input"]), "filtered": float(r["filtered"]),
                                                "denoised": float(r["denoised"]), "merged": float(r["merged"]),
                                                "nonchim": float(r["non-chimeric"])}
    a_stats = {}
    with open(AMP / "DADA2_stats.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            a_stats[r["sample"]] = {"input": float(r["DADA2_input"]), "filtered": float(r["filtered"]),
                                    "denoised": float(r["denoisedF"]), "merged": float(r["merged"]),
                                    "nonchim": float(r["nonchim"])}
    rep.append("1. Read retention (share of input reads at each stage; median per run)")
    rep.append(f"   {'run':9s} {'n':>4s}  {'stage':9s} {'QIIME 2':>8s} {'ampliseq':>9s}")
    dropped = sorted(set(s1_stats) - set(a_stats))
    extra = sorted(set(a_stats) - set(s1_stats))
    input_diff = [s for s in s1_stats if s in a_stats and s1_stats[s]["input"] != a_stats[s]["input"]]
    for run in RUNS:
        ss = [s for s in s1_stats if s1_stats[s]["run"] == run and s in a_stats]
        for st in ("filtered", "denoised", "merged", "nonchim"):
            a = statistics.median(s1_stats[s][st] / s1_stats[s]["input"] for s in ss)
            b = statistics.median(a_stats[s][st] / a_stats[s]["input"] for s in ss)
            rep.append(f"   {run:9s} {len(ss):4d}  {st:9s} {a:8.3f} {b:9.3f}")
    ratio = {s: a_stats[s]["nonchim"] / s1_stats[s]["nonchim"] for s in s1_stats if s in a_stats}
    rep.append(f"   final reads ampliseq / QIIME 2 per sample: median {statistics.median(ratio.values()):.3f}, "
               f"range {min(ratio.values()):.3f}-{max(ratio.values()):.3f}")
    rep.append(f"   samples: QIIME 2 {len(s1_stats)}, ampliseq {len(a_stats)}; dropped by ampliseq {dropped or 'none'}; "
               f"extra {extra or 'none'}; different input counts: {len(input_diff)}")

    # tables by sequence
    s1_tab, s1_seqs = stage1(S1)
    a_tab = table_tsv(AMP / "ASV_table.tsv")
    a_seqs = fasta(AMP / "ASV_seqs.fasta")
    S = by_sequence(s1_tab, s1_seqs)
    A = by_sequence(a_tab, a_seqs)
    a_depth = {s: sum(v.values()) for s, v in A.items()}
    miss_frozen = sorted(s for s in frozen_ids if a_depth.get(s, 0) < 6400)
    rep.append(f"   frozen 102 samples with >= 6,400 reads in ampliseq: {len(frozen_ids) - len(miss_frozen)} of {len(frozen_ids)}"
               + (f"; MISSING/BELOW DEPTH: {miss_frozen}" if miss_frozen else ""))
    rep.append("")

    # ---- 2. exact ASV agreement ----
    sset, aset = set(s1_seqs.values()), set(a_seqs.values())
    shared = sset & aset
    rc_only = sum(1 for x in aset - sset if rc(x) in sset)
    s_reads = Counter(); a_reads = Counter()
    for d in S.values():
        s_reads.update(d)
    for d in A.values():
        a_reads.update(d)
    s_tot, a_tot = sum(s_reads.values()), sum(a_reads.values())
    rep.append("2. Exact ASV agreement (full sequence)")
    for label, seqs, reads in [("QIIME 2", sset, s_reads), ("ampliseq", aset, a_reads)]:
        lc = Counter(len(x) for x in seqs)
        lr = Counter()
        for x, n in reads.items():
            lr[len(x)] += n
        tot = sum(reads.values())
        common = ", ".join(f"{L} bp: {lc[L]} ASVs / {lr[L] / tot:.4f} reads" for L in sorted(lc, key=lambda L: -lr[L])[:3])
        short = sorted(L for L in lc if L < 250)
        rep.append(f"   {label}: {len(seqs)} ASVs; lengths {min(lc)}-{max(lc)}; main lengths {common}; "
                   f"ASVs < 250 bp: {sum(lc[L] for L in short)} carrying {sum(lr[L] for L in short) / tot:.5f} of reads"
                   + (f" (lengths {short})" if short else ""))
    rep.append(f"   shared ASVs {len(shared)}; QIIME 2 only {len(sset - shared)}; ampliseq only {len(aset - shared)}; "
               f"ampliseq-only ASVs whose reverse complement is a QIIME 2 ASV: {rc_only}")
    rep.append(f"   reads in shared ASVs: QIIME 2 {sum(s_reads[x] for x in shared) / s_tot:.4f}; "
               f"ampliseq {sum(a_reads[x] for x in shared) / a_tot:.4f}")
    common_samples = sorted(set(S) & set(A))
    ps_s = [sum(v for x, v in S[s].items() if x in shared) / sum(S[s].values()) for s in common_samples]
    ps_a = [sum(v for x, v in A[s].items() if x in shared) / sum(A[s].values()) for s in common_samples]
    rep.append(f"   per sample, share of reads in shared ASVs: QIIME 2 median {statistics.median(ps_s):.4f} (min {min(ps_s):.4f}); "
               f"ampliseq median {statistics.median(ps_a):.4f} (min {min(ps_a):.4f})")
    for label, reads, other, tot in [("QIIME 2 only", s_reads, aset, s_tot), ("ampliseq only", a_reads, sset, a_tot)]:
        top = [(x, n) for x, n in reads.most_common() if x not in other][:3]
        rep.append(f"   largest {label}: " + "; ".join(f"{len(x)} bp, {n / tot:.5f} of reads" for x, n in top))
    rep.append("")

    # ---- 3. same-sample Bray-Curtis ----
    s1b_tab, s1b_seqs = stage1(S1B)
    SB = by_sequence(s1b_tab, s1b_seqs)
    d_wf = {s: bc(S[s], A[s]) for s in common_samples}
    d_ref = {s: bc(S[s], SB[s]) for s in sorted(set(S) & set(SB))}
    rep.append("3. Same-sample Bray-Curtis (unrarefied relative abundances, union of ASVs by sequence)")
    rep.append(f"   QIIME 2 vs ampliseq:              {summary(list(d_wf.values()))}")
    rep.append(f"   reference, QIIME 2 200/150 vs 218/166 (recomputed): {summary(list(d_ref.values()))}")
    for run in RUNS:
        d0 = [d_wf[s] for s in common_samples if meta[s]["run"] == run and meta[s]["source"] == "faeces" and float(meta[s]["day"]) == 0]
        al = [d_wf[s] for s in common_samples if meta[s]["run"] == run]
        rep.append(f"   {run}: all samples median {statistics.median(al):.4f}, max {max(al):.4f}"
                   + (f"; day-0 faecal samples n={len(d0)}, median {statistics.median(d0):.4f}, max {max(d0):.4f}" if d0 else "; no day-0 samples"))
    worst = sorted(d_wf.items(), key=lambda x: -x[1])[:5]
    rep.append("   largest: " + "; ".join(f"{s} ({meta[s]['run']}, day {meta[s]['day']}) {v:.4f}" for s, v in worst))
    rep.append("")

    # ---- 4. family composition and Muribaculaceae traces ----
    s_fam = families(Path("results/taxonomy") / S1 / "export" / "taxonomy.tsv")
    a_fam = families(Path("results/ampliseq_q2/taxonomy_export/taxonomy.tsv"))
    seq_fam_s = {s1_seqs[i]: f for i, f in s_fam.items() if i in s1_seqs}
    seq_fam_a = {a_seqs[i]: f for i, f in a_fam.items() if i in a_seqs}
    same_tax = sum(1 for x in shared if seq_fam_s.get(x, ("", ""))[1] == seq_fam_a.get(x, ("", ""))[1])

    def fam_profile(T, fm):
        out = {}
        for s, d in T.items():
            c = Counter()
            for x, n in d.items():
                c[fm.get(x, ("(no taxonomy)", ""))[0]] += n
            out[s] = c
        return out
    FS, FA = fam_profile(S, seq_fam_s), fam_profile(A, seq_fam_a)
    tot_s, tot_a = Counter(), Counter()
    for c in FS.values():
        tot_s.update(c)
    for c in FA.values():
        tot_a.update(c)
    rep.append("4. Family composition (same SILVA 144 classifier and settings)")
    rep.append(f"   shared ASVs with identical taxon string: {same_tax} of {len(shared)}")
    rep.append(f"   {'family':32s} {'QIIME 2':>8s} {'ampliseq':>9s}")
    for f, _ in tot_s.most_common(10):
        rep.append(f"   {f[:32]:32s} {tot_s[f] / sum(tot_s.values()):8.4f} {tot_a[f] / sum(tot_a.values()):9.4f}")
    fam_bc = [bc(FS[s], FA[s]) for s in common_samples]
    rep.append(f"   per-sample family-level Bray-Curtis between workflows: {summary(fam_bc)}")
    pd_win = [(m, mice[m]["cage"], s) for m in mice if mice[m]["diet"] == "PD" for d, s in mice[m]["days"].items() if d in WIN]
    trace_rows = []
    rep.append(f"   Muribaculaceae in the {len(pd_win)} eligible PD window samples:")
    for label, T, fm, other in [("QIIME 2", S, seq_fam_s, aset), ("ampliseq", A, seq_fam_a, sset)]:
        det_s, det_c, reads, asvs = 0, set(), 0, Counter()
        for m, c, s in pd_win:
            mur = {x: n for x, n in T.get(s, {}).items() if fm.get(x, ("",))[0] == "Muribaculaceae"}
            n = sum(mur.values())
            if n:
                det_s += 1; det_c.add(c)
            reads += n
            asvs.update(mur)
            for x, k in mur.items():
                trace_rows.append((label, s, m, c, x, len(x), int(k), k / sum(T[s].values()),
                                   "yes" if x in other else "no"))
        rep.append(f"     {label}: detected in {det_s} of {len(pd_win)} samples, cages {sorted(det_c) or 'none'}; "
                   f"{int(reads)} reads; {len(asvs)} distinct ASVs"
                   + (f" (top: " + "; ".join(f"{n:.0f} reads, {len(x)} bp, in other workflow: {'yes' if x in other else 'no'}"
                                             for x, n in asvs.most_common(3)) + ")" if asvs else ""))
    rep.append("")

    # ---- 5. frozen biological analysis on ampliseq ----
    rep.append("5. Frozen biological analysis (D27, D28 settings) on the ampliseq outputs")
    ok_mice = {m: d for m, d in mice.items() if all(s not in miss_frozen for s in d["days"].values())}
    if len(ok_mice) < len(mice):
        rep.append(f"   INCOMPLETE: {len(mice) - len(ok_mice)} mice lack a required sample; the analysis below uses the "
                   f"common subset of {len(ok_mice)} mice, cages {sorted({d['cage'] for d in ok_mice.values()})}")
    ids = {s for d in ok_mice.values() for s in d["days"].values()}
    res = {}
    for label, folder in [("QIIME 2", Path("data/qiime/diversity") / S1), ("ampliseq", AQ)]:
        dm = read_dm(one(folder / "distance_matrices", BC_NAMES), ids)
        mv = {m: mean([dm[d["days"][float(x)]][d["days"][0.0]] for x in WIN]) for m, d in ok_mice.items()}
        res[label] = cage_weighted(ok_mice, mv)
    rs, ra = res["QIIME 2"], res["ampliseq"]
    rep.append(f"   QIIME 2 recomputed here: delta {rs['delta']:+.4f} (frozen value {FROZEN_DELTA:+.4f})")
    rep.append(f"   ampliseq:                delta {ra['delta']:+.4f}  (STD {ra['diet']['STD']:.4f}, PD {ra['diet']['PD']:.4f})")
    rep.append("   cage values  " + "  ".join(f"{ra['diet_of'][c]} {c}: QIIME 2 {rs['cage'][c]:.4f} / ampliseq {ra['cage'][c]:.4f} (n={ra['n'][c]})" for c in ra["cages"]))
    same_dir = (ra["delta"] < 0) == (FROZEN_DELTA < 0)
    within = abs(ra["delta"] - FROZEN_DELTA) <= TOL
    rep.append(f"   HEADLINE RULE: same direction {'yes' if same_dir else 'NO'}; |{ra['delta']:+.4f} - ({FROZEN_DELTA:+.4f})| = "
               f"{abs(ra['delta'] - FROZEN_DELTA):.4f} {'<=' if within else '>'} {TOL} -> "
               + ("RESULT HOLDS under ampliseq" if same_dir and within else "RESULT DOES NOT HOLD under the rule")
               + ("" if len(ok_mice) == len(mice) else " (common subset; see INCOMPLETE above)"))
    rep.append(f"   pattern check: both PD cages below all three STD cages: QIIME 2 {'yes' if pd_below_std(rs) else 'no'}; "
               f"ampliseq {'yes' if pd_below_std(ra) else 'no'}")
    rep.append("   window levels, days 8-12 (cage-equal; script 20 A):")
    win_ids = {m: [d["days"][float(x)] for x in WIN] for m, d in ok_mice.items()}
    for metric, fmt in [("observed_features", ".1f"), ("faith_pd", ".2f")]:
        line = []
        for label, folder in [("QIIME 2", Path("data/qiime/diversity") / S1), ("ampliseq", AQ)]:
            al = read_alpha(folder / "alpha_diversities" / f"{metric}.qza")
            r = cage_weighted(ok_mice, {m: mean([al[s] for s in win_ids[m]]) for m in ok_mice})
            line.append(f"{label} STD {r['diet']['STD']:{fmt}} PD {r['diet']['PD']:{fmt}} (PD-STD {r['delta']:+{fmt}})")
        rep.append(f"     {metric:18s} " + "; ".join(line))
    for f in FAMS:
        line = []
        for label, T, fm in [("QIIME 2", S, seq_fam_s), ("ampliseq", A, seq_fam_a)]:
            sh = {s: sum(n for x, n in T[s].items() if fm.get(x, ("",))[0] == f) / sum(T[s].values()) for s in ids}
            r = cage_weighted(ok_mice, {m: mean([sh[s] for s in win_ids[m]]) for m in ok_mice})
            line.append(f"{label} STD {r['diet']['STD']:.3f} PD {r['diet']['PD']:.3f}")
        rep.append(f"     {f:18s} " + "; ".join(line))
    rep += ["", "Interpretation limit (D29): agreement supports reproducibility across the two workflow",
            "implementations; the day-0 sequencing-run difference still limits attribution to diet."]
    (OUT / "stage2_comparison_report.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))

    with open(OUT / "per_sample.tsv", "w") as fh:
        fh.write("sample\trun\tday\tsource\tqiime2_nonchim\tampliseq_nonchim\tratio\tbray_curtis_workflows\t"
                 "bray_curtis_reference_arms\tshared_read_share_qiime2\tshared_read_share_ampliseq\n")
        for s, a, b in zip(common_samples, ps_s, ps_a):
            fh.write(f"{s}\t{meta[s]['run']}\t{meta[s]['day']}\t{meta[s]['source']}\t{s1_stats[s]['nonchim']:.0f}\t"
                     f"{a_stats[s]['nonchim']:.0f}\t{ratio[s]:.4f}\t{d_wf[s]:.5f}\t{d_ref.get(s, float('nan')):.5f}\t{a:.5f}\t{b:.5f}\n")
    with open(OUT / "muribaculaceae_pd_window.tsv", "w") as fh:
        fh.write("workflow\tsample\tmouse\tcage\tsequence\tlength\treads\trelative_abundance\tin_other_workflow\n")
        for r in trace_rows:
            fh.write("\t".join(str(x) if not isinstance(x, float) else f"{x:.6f}" for x in r) + "\n")
    with open(OUT / "cage_values.tsv", "w") as fh:
        fh.write("cage_id\tdiet\tn_mice\tqiime2\tampliseq\n")
        for c in ra["cages"]:
            fh.write(f"{c}\t{ra['diet_of'][c]}\t{ra['n'][c]}\t{rs['cage'][c]:.6f}\t{ra['cage'][c]:.6f}\n")
    print(f"\nWrote {OUT}/stage2_comparison_report.txt, per_sample.tsv, muribaculaceae_pd_window.tsv, cage_values.tsv")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        sys.exit(__doc__)
    main()

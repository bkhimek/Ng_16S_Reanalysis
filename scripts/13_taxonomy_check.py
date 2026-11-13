#!/usr/bin/env python3
"""13_taxonomy_check.py - quality check of the taxonomy assignment (D24).

Per arm, read-weighted over all samples together (no diet or time contrast):
  1. Share of reads Unassigned, or assigned to Bacteria / Archaea / Eukaryota.
  2. Share of reads named at each rank (phylum ... genus). The SILVA 144
     classifier uses rank propagation: a reference sequence without a genus
     carries its family name down into the genus slot. A "genus" equal to the
     family name, or labelled uncultured/unidentified, is therefore counted
     separately as NOT resolved to genus.
  3. Non-target reads: mitochondria and chloroplast (matched at any rank),
     which are host/plant organelle 16S, not gut bacteria.
  4. The 12 most abundant families overall (sanity check that the profile
     looks like mouse gut, e.g. Muribaculaceae = the paper's S24-7), and the
     median classification confidence.
  5. With two arms: reads whose ASV sequence gets the same taxon in both arms.

Usage (repository root):
  python3 scripts/13_taxonomy_check.py trunc200_150 [trunc218_166]
"""
import csv
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

RANKS = ["d", "k", "p", "c", "o", "f", "g"]   # SILVA 144 adds a kingdom rank (k__)
NAMES = {"d": "domain", "k": "kingdom", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus"}
# labels that do not resolve a rank: placeholders, and SILVA 144 "--other" groups
VAGUE = ("uncultured", "unidentified", "unknown", "metagenome", "incertae", "--other")


def parse(taxon):
    out = {}
    for part in taxon.split(";"):
        part = part.strip()
        if len(part) > 3 and part[1:3] == "__":
            name = part[3:].strip()
            if name:
                out[part[0]] = name
    return out


def load_table(arm):
    with open(Path("results/dada2") / arm / "table_export" / "feature-table.tsv") as fh:
        rows = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("# Constructed")]
    return {r[0]: sum(float(x) for x in r[1:]) for r in rows[1:]}


def load_fasta(arm):
    seqs, name = {}, None
    for line in open(next((Path("results/dada2") / arm / "repseqs_export").glob("*.fasta"))):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = ""
        elif name:
            seqs[name] += line
    return seqs


def load_tax(arm):
    with open(Path("results/taxonomy") / arm / "export" / "taxonomy.tsv") as fh:
        return {r["Feature ID"]: (r["Taxon"], float(r["Confidence"])) for r in csv.DictReader(fh, delimiter="\t")}


def check(arm):
    reads, tax = load_table(arm), load_tax(arm)
    total = sum(reads.values())
    missing = set(reads) - set(tax)
    rep = [f"##### Taxonomy {arm}: {len(tax):,} ASVs classified; reads {int(total):,}"
           + (f"; WARNING {len(missing)} table ASVs without taxonomy" if missing else "")]
    dom, rank_named, genus_strict, fam = Counter(), Counter(), 0.0, Counter()
    mito = chloro = 0.0
    conf = []
    for a, n in reads.items():
        t, c = tax.get(a, ("Unassigned", 0.0))
        conf.append(c)
        p = parse(t)
        dom[p.get("d", "Unassigned")] += n
        for r in RANKS[1:]:
            if r in p:
                rank_named[r] += n
        g, f = p.get("g", ""), p.get("f", "")
        if g and g != f and not any(v in g.lower() for v in VAGUE):
            genus_strict += n
        tl = t.lower()   # match at any rank: label placement may differ between SILVA releases
        if "mitochondri" in tl:
            mito += n
        if "chloroplast" in tl:
            chloro += n
        fam[f or f"(no family: {p.get('o', p.get('p', p.get('d', 'Unassigned')))})"] += n
    rep.append("1. Domain: " + ", ".join(f"{k} {v/total:.4f}" for k, v in dom.most_common()))
    rep.append("2. Reads named at rank: " + ", ".join(f"{NAMES[r]} {rank_named[r]/total:.4f}" for r in RANKS[1:]))
    rep.append(f"   genus, strict (not a propagated family name, not uncultured/unidentified/Incertae_Sedis/--other): {genus_strict/total:.4f}")
    rep.append(f"3. Mitochondria {mito/total:.5f}; chloroplast {chloro/total:.5f}")
    phy = Counter()
    for a, n in reads.items():
        phy[parse(tax.get(a, ("Unassigned", 0.0))[0]).get("p", "(no phylum)")] += n
    rep.append("   phyla (SILVA 144 names): " + ", ".join(f"{k} {v/total:.4f}" for k, v in phy.most_common(8)))
    rep.append("4. Most abundant families (share of reads):")
    for k, v in fam.most_common(12):
        rep.append(f"     {v/total:.4f}  {k}")
    rep.append(f"   median confidence (per ASV): {statistics.median(conf):.3f}")
    return rep, tax


def main(arms):
    taxa = {}
    for arm in arms:
        rep, taxa[arm] = check(arm)
        Path("results/taxonomy", arm).mkdir(parents=True, exist_ok=True)
        (Path("results/taxonomy") / arm / "taxonomy_check_report.txt").write_text("\n".join(rep) + "\n")
        print("\n".join(rep) + "\n")
    if len(arms) == 2:
        a, b = arms
        sa, sb = load_fasta(a), load_fasta(b)
        ra = load_table(a)
        by_seq_b = {sb[i]: taxa[b][i][0] for i in taxa[b] if i in sb}
        shared = same = 0.0
        for i, n in ra.items():
            s = sa.get(i)
            if s in by_seq_b:
                shared += n
                if taxa[a].get(i, ("",))[0] == by_seq_b[s]:
                    same += n
        print(f"5. {a} reads in ASVs also present in {b}: {shared/sum(ra.values()):.4f}; "
              f"of these, identical taxon string: {same/shared:.4f}")


if __name__ == "__main__":
    if not 1 <= len(sys.argv[1:]) <= 2:
        sys.exit(__doc__)
    main(sys.argv[1:])

#!/usr/bin/env python3
"""17_freeze_check.py - sample list and artifact record for the frozen analysis (D28).

Reads ONLY sample metadata, sample IDs and file identities; no diversity value
is read or printed, so running it does not reveal any diet or time result.

1. Primary sample set (metadata/primary_analysis_samples.tsv):
   faecal samples, minus failed libraries (D19, metadata/excluded_samples.tsv)
   and the D23 pair; a mouse is eligible if it has exactly one day-0 faecal
   sample and exactly one sample on each of days 8, 9, 10, 11, 12.
   Day 13 samples of eligible mice are listed with role 'sens_day13' (window
   8-13 sensitivity analysis).
2. Presence check: every listed sample must be in each arm's distance matrices
   (header line only), alpha vectors (ID column only) and DADA2 table.
3. Artifact record (docs/frozen_artifacts.tsv): path, UUID, MD5 and whether the
   rarefaction seed appears in the artifact's provenance.

Usage (repository root):
  python3 scripts/17_freeze_check.py            (both arms)
Exit code 1 if any check fails.
"""
import csv
import hashlib
import io
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

ARMS = ["trunc200_150", "trunc218_166"]
D23 = {"H3058", "H3072"}           # same mouse-day label, uncertain identity (D23)
WINDOW = [8, 9, 10, 11, 12]
SEED = "20261113"
META = Path("metadata/humanised_sample_metadata.tsv")
EXCL = Path("metadata/excluded_samples.tsv")
OUT_SAMPLES = Path("metadata/primary_analysis_samples.tsv")
OUT_ART = Path("docs/frozen_artifacts.tsv")


def read_meta():
    with open(META) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if not r["sample-id"].startswith("#")]
    with open(EXCL) as fh:
        excl = {r["sample-id"] for r in csv.DictReader(fh, delimiter="\t")}
    return rows, excl | D23


def eligible(rows, excluded):
    by_mouse = defaultdict(lambda: defaultdict(list))
    info = {}
    for r in rows:
        info[r["mouse_id"]] = (r["cage_id"], r["diet"])
        if r["source"] != "faeces" or r["sample-id"] in excluded:
            continue
        by_mouse[r["mouse_id"]][float(r["day"])].append(r)
    out, problems = [], []
    for m in sorted(info, key=lambda x: int(x.split("m")[-1])):
        days = by_mouse[m]
        need = [0] + WINDOW
        if any(len(days.get(float(d), [])) > 1 for d in need + [13]):
            problems.append(f"{m}: more than one sample on a needed day")
            continue
        if all(len(days.get(float(d), [])) == 1 for d in need):
            for d in need:
                out.append((days[float(d)][0], "baseline" if d == 0 else "window"))
            if len(days.get(13.0, [])) == 1:
                out.append((days[13.0][0], "sens_day13"))
    return out, info, problems


def qza_ids_tsv(path, member_suffix, header_line):
    """Sample IDs from a .qza without reading values: header line (distance
    matrix) or first column (alpha vector)."""
    with zipfile.ZipFile(path) as zf:
        name = next(n for n in zf.namelist() if n.endswith(member_suffix))
        with zf.open(name) as fh:
            text = io.TextIOWrapper(fh, encoding="utf-8")
            if header_line:
                return set(text.readline().rstrip("\n").split("\t")[1:])
            ids = set()
            next(text)
            for line in text:
                ids.add(line.split("\t", 1)[0])
            return ids


def identity(path):
    md5 = hashlib.md5(path.read_bytes()).hexdigest()
    with zipfile.ZipFile(path) as zf:
        uuid = zf.namelist()[0].split("/")[0]
        seed = any(SEED in zf.read(n).decode("utf-8", "ignore")
                   for n in zf.namelist() if "/provenance/" in n and n.endswith(".yaml"))
    return uuid, md5, seed


def main():
    ok = True
    rows, excluded = read_meta()
    samples, info, problems = eligible(rows, excluded)
    for p in problems:
        print("PROBLEM", p); ok = False
    OUT_SAMPLES.parent.mkdir(exist_ok=True)
    with open(OUT_SAMPLES, "w") as fh:
        fh.write("sample-id\tmouse_id\tcage_id\tdiet\tday\trun\trole\n")
        for r, role in samples:
            fh.write("\t".join([r["sample-id"], r["mouse_id"], r["cage_id"], r["diet"], r["day"], r["run"], role]) + "\n")
    prim = [r for r, role in samples if role != "sens_day13"]
    mice = sorted({r["mouse_id"] for r in prim}, key=lambda x: int(x.split("m")[-1]))
    print(f"Primary set: {len(mice)} mice, {len(prim)} samples (day 0 + days 8-12); "
          f"day-13 samples for sensitivity: {sum(role == 'sens_day13' for _, role in samples)}")
    per_cage = Counter(info[m] for m in mice)
    for (cage, diet), n in sorted(per_cage.items(), key=lambda x: (x[0][1], x[0][0])):
        print(f"  {diet} {cage}: {n} mice")
    not_elig = [m for m in info if m not in mice]
    print(f"  not eligible: {', '.join(sorted(not_elig, key=lambda x: int(x.split('m')[-1]))) or 'none'}")
    runs = Counter((r["diet"], "baseline" if float(r["day"]) == 0 else "window", r["run"]) for r in prim)
    print("  runs: " + "; ".join(f"{d} {role} {run}: {n}" for (d, role, run), n in sorted(runs.items())))
    print(f"Wrote {OUT_SAMPLES}")

    ids_needed = {r["sample-id"] for r, _ in samples}
    art = []
    for arm in ARMS:
        Q = Path("data/qiime/diversity") / arm
        inputs = [Path("data/qiime/dada2") / arm / "humanised_table.qza",
                  Path("data/qiime/phylogeny") / arm / "rooted_tree.qza"]
        dms = sorted((Q / "distance_matrices").glob("*.qza"))
        alphas = sorted((Q / "alpha_diversities").glob("*.qza"))
        if not dms or not alphas:
            print(f"FAIL {arm}: diversity outputs not found under {Q}"); ok = False
            continue
        arm_ok = True
        for p in dms + alphas:
            got = qza_ids_tsv(p, "/data/distance-matrix.tsv", True) if p in dms \
                else qza_ids_tsv(p, "/data/alpha-diversity.tsv", False)
            miss = ids_needed - got
            if miss:
                print(f"FAIL {p}: {len(miss)} listed samples missing: {', '.join(sorted(miss)[:10])}"); ok = arm_ok = False
        tab = Path("results/dada2") / arm / "table_export" / "feature-table.tsv"
        with open(tab) as fh:
            for line in fh:
                if line.startswith("#OTU ID"):
                    miss = ids_needed - set(line.rstrip("\n").split("\t")[1:]); break
        if miss:
            print(f"FAIL {tab}: {len(miss)} listed samples missing"); ok = arm_ok = False
        if arm_ok:
            print(f"{arm}: all {len(ids_needed)} listed samples present in {len(dms)} distance matrices, "
                  f"{len(alphas)} alpha vectors and the DADA2 table")
        for p in inputs + dms + alphas:
            uuid, md5, seed = identity(p)
            art.append((arm, str(p), uuid, md5, ("yes" if seed else "NO") if p not in inputs else "n/a"))
            if p not in inputs and not seed:
                print(f"WARN {p}: seed {SEED} not found in provenance text (check with qiime tools peek / provenance)")
    OUT_ART.parent.mkdir(exist_ok=True)
    with open(OUT_ART, "w") as fh:
        fh.write("arm\tpath\tuuid\tmd5\tseed_20261113_in_provenance\n")
        for a in art:
            fh.write("\t".join(a) + "\n")
    print(f"Wrote {OUT_ART} ({len(art)} artifacts)")
    print("ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        sys.exit(__doc__)
    main()

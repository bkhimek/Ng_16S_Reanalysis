# Frozen analysis specification v1: humanised experiment (D28)

**Frozen:** 13 November 2026. Git tag: `analysis-frozen-v1`.
**State at freeze:** diversity outputs exist but no diversity value has been computed into a diet or time comparison or viewed by diet or time. Sample eligibility was set from the sample list only (script 17).
**Rule:** any analysis not defined here is exploratory and labelled as such. Any change to this file after the tag is a dated amendment with its reason, and results affected by it are labelled.

## 1. Data and processing

| Item | Specification | Source |
|---|---|---|
| Experiment | Humanised only; conventional on hold | D3, D4 |
| Primary processing arm | DADA2 200/150 | D20a |
| Sensitivity arm | DADA2 218/166 | D21 |
| Rarefaction | q2-boots 2026.10.0 core-metrics: depth 6,400, n = 100, without replacement, random seed 20261113; alpha averaged by mean; distance matrices by medoid (the single most representative of the 100 matrices) | D27 |
| Artifact identity | Exact files (path, UUID, MD5) listed in `docs/frozen_artifacts.tsv`, made by `scripts/17_freeze_check.py` | this spec |

The 100 rarefactions address subsampling variability. Biological replication comes from the cages, not from the rarefactions.

## 2. Samples

Listed in `metadata/primary_analysis_samples.tsv`.

- Faecal samples only. Cecal samples are excluded.
- Excluded: failed libraries H3106 and H2956 (D19), and the H_m19 day-5 pair H3058/H3072 (D23).
- **Baseline:** the single faecal day-0 sample of each mouse (D5). The day-0.5 sample is not a baseline.
- **Eligible mouse:** a valid day-0 sample and one sample on each of days 8, 9, 10, 11 and 12 (all five window days required).
- At freeze: 17 mice and 102 samples.

  | Diet | Cage | Mice |
  |---|---|---|
  | STD | H_c1 | 4 |
  | STD | H_c2 | 4 |
  | STD | H_c3 | 2 |
  | PD | H_c4 | 4 |
  | PD | H_c5 | 3 |

  Not eligible (no window samples): H_m14, H_m16, H_m17. The 1-, 3- and 5-day rules give the same 17 mice, so the rule was fixed at no cost in sample size.
- Day-13 samples of the same 17 mice (17 samples) are used only in sensitivity analysis S4.

## 3. Primary outcome and estimate

1. For each eligible mouse and each day *d* from 8 to 12: Bray–Curtis distance between the day-*d* sample and the same mouse's day-0 sample, taken from the **200/150 rarefied medoid Bray–Curtis matrix**.
2. Mouse value: mean of its five distances.
3. Cage value: mean of its mice.
4. Diet value: unweighted mean of its cage values (**each cage weighs equally**; this estimates the difference between average cages).
5. **Primary estimate:** Δ = mean of PD cage values − mean of STD cage values. A positive Δ means PD samples stay further from their own baseline, which is read as less recovery.

**Reported:** Δ; all five cage values with the number of mice in each; and the mouse values plotted within cages.

**No p-value is reported.** STD distances all cross sequencing runs (run 1 baseline to run 2 window), while PD distances all stay within run 2. Cage labels are therefore not exchangeable, and relabelling cannot recreate comparable baselines. Five cages also allow only 10 allocations.

**Supplementary (descriptive):** Δ under all 10 allocations of the five cages into three "STD" and two "PD" labels, to show where the observed contrast sits among hypothetical relabellings. This is not a test.

## 4. Secondary measures

Descriptive, using the same 17 mice and the same mouse → cage → diet averaging, with the PD − STD difference of cage means.

| Measure | Definition |
|---|---|
| Observed ASVs and Faith's PD | Rarefied mean value (200/150); per mouse, mean over days 8–12 minus day 0. Full trajectories over all days are shown as illustration. |
| Weighted UniFrac | Distance from day 0 as in §3, from the 200/150 rarefied medoid weighted UniFrac matrix. |
| Muribaculaceae, Bacteroidaceae, Akkermansiaceae | Share of reads per sample, unrarefied 200/150 table collapsed to family with the SILVA 144 taxonomy (D24); per mouse, mean over days 8–12 minus day 0. |

## 5. Sensitivity analyses

Each sensitivity analysis computes Δ exactly as in §3 on the same mice, days and cages, and is reported beside the primary Δ.

| ID | Analysis | Specification |
|---|---|---|
| S1 | Truncation (D21) | 218/166 rarefied medoid Bray–Curtis. Report Δ(200/150), Δ(218/166), their difference, and the cage values from both arms. A change in direction or a large change in magnitude is reported as such. Observed ASVs and Faith's PD are also compared between arms (D21). |
| S2 | Unrarefied Bray–Curtis | 200/150 table; each sample's counts divided by its total reads (relative abundances), so library size does not enter the distance directly. |
| S3 | Aitchison | 200/150 unrarefied counts of the 102 primary samples. Keep ASVs with ≥ 1 read in ≥ 10% of these samples (≥ 11 of 102). The filter is applied once to all 102 together, without regard to diet. Add a pseudocount of 1 to every count, apply a centred log-ratio transform per sample, and compute the Euclidean distance. |
| S4 | Window days 8–13 | Primary set plus day 13; mouse value is the mean of six distances. |

If S1 and the primary result agree closely, the cage-level pattern can be described as stable to the tested truncation change. That does not remove the run confounding.

## 6. Interpretation limits (stated with every result)

- **Baseline/run confounding:** STD baselines were sequenced in run 1 and PD baselines in run 2. Δ can partly reflect run, so attributing it to diet is limited.
- Two to three cages per diet; mice within a cage are not independent.
- No blanks, no mock community, and no confirmed technical replicates.

Conclusions describe the cage-level pattern and its consistency across cages and sensitivity analyses; they are not a test of a diet effect.

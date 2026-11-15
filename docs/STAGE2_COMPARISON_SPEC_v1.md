# Stage 2 comparison specification v1: QIIME 2 vs nf-core/ampliseq (D30)

**Status:** fixed before the full ampliseq run, and before any ampliseq ASV table is examined. Git tag: `stage2-spec-v1`. Changes after the tag are dated amendments with reasons.

**Question:** Does the frozen cage-level result hold when the same reads are processed through a matched nf-core/ampliseq workflow?

**Design:**

- Both workflows use DADA2. Inputs, DADA2 settings (`docs/STAGE2_PARAMETER_MATCH.md`), taxonomy and downstream analysis are held constant.
- Agreement would support reproducibility across the two workflow implementations.
- It would not remove the day-0 sequencing-run difference, which still limits attribution to diet (D29).

## Inputs and versions

- Identical checksummed per-sample FASTQs and run assignments for the 451 samples Stage 1 denoised (script 21; `metadata/per_sample_fastq_sha256.tsv`).
- nf-core/ampliseq 3.1.0 (`-r 3.1.0`); `-profile docker`; `config/ampliseq/params_matched.yaml` and `matched.config`.
- Nextflow version and the image digest of every container used are recorded.

## Comparisons

Stage 1 means the 200/150 arm throughout.

| # | Comparison | Specification |
|---|---|---|
| 1 | **Read retention** | Read-pair counts at equivalent stages (input, filtered, denoised, merged, non-chimeric), same input denominator, per sample and per run. List any sample dropped by ampliseq. Flag any missing day-0 or days 8–12 sample among the frozen 17 mice and 102 samples. |
| 2 | **Exact ASV agreement** | Match full sequences in the same orientation, after checking that both workflows give the same amplicon boundaries (length distribution). Report the number of shared and workflow-specific ASVs, and the share of retained reads in shared ASVs in each workflow, overall and per sample (median, minimum). This shows whether unmatched ASVs are mostly rare or carry substantial reads in particular samples. |
| 3 | **Same-sample Bray–Curtis** | Between the two workflows, for each of the 451 samples: unrarefied relative abundances over the union of ASVs, with zeros for ASVs absent from one table. Report the median, 95th percentile and maximum, overall and with day-0 samples shown by run. Stage 1 reference points are the median 0.024 and maximum 0.044 between truncation arms, computed the same way on the same sample set. If the normalisation or sample set differs, the reference is recomputed before comparison. |
| 4 | **Family composition** | ampliseq ASVs classified with the signed SILVA 144 V4 classifier (D24), with the settings of script 12 (confidence, rank handling), family aggregation and the same read denominator. Read-weighted family shares are compared overall and per sample. **Muribaculaceae traces:** exact ASV sequences, raw counts and relative abundances in each eligible PD window sample, and how many of the 35 PD window samples and how many PD cages contain them, in both workflows. |
| 5 | **Frozen biological analysis** | Rerun on the ampliseq table with the frozen downstream settings (D27, D28): tree as script 14; q2-boots depth 6,400, n = 100, seed 20261113, mean / medoid; the same 17 mice and cage weighting. Compare the primary contrast, all five cage values, and the late-window richness and family levels (script 20 A). If ampliseq loses a required sample, the original comparison is flagged incomplete, and any common-subset analysis is labelled with its revised sample and cage counts. |

## Headline decision rule

The frozen result holds under ampliseq if both conditions are met:

1. **Same direction:** ampliseq Δ is below 0, like the frozen QIIME 2 estimate Δ = −0.1454.
2. **Within tolerance:** |Δ(ampliseq) − (−0.1454)| ≤ **0.02** Bray–Curtis units.

The tolerance is a project-level reporting convention for this benchmark. It is about 14% of the Stage 1 contrast and roughly four times the 0.0046 shift between truncation settings.

**Separate pattern check:** whether both PD cages still fall below all three STD cages.

**Reporting:** comparisons 1–4 are descriptive. A result outside the rule is reported as such, with the comparison step (1–4) that explains it where possible. Thresholds are not changed after the results are seen.

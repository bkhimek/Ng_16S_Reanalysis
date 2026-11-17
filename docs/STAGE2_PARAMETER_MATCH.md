# Stage 2 parameter match: QIIME 2 (Stage 1) vs nf-core/ampliseq 3.1.0 (D30)

**Sources read on 15 November 2026:**

- **Stage 1:** q2-dada2 2026.10.0, `q2_dada2/assets/run_dada.R` (tag `2026.10.0`), called by `scripts/08_dada2_humanised.sh`.
- **Stage 2:** nf-core/ampliseq 3.1.0, commit `b9bd579c0480ee23ecf24655f51c0c67777cd6b3` (28 Oct 2026). Files read: `conf/modules.config`, `nextflow_schema.json` and `modules/local/dada2_*`.

> **Amendment, 17 November 2026.** The QIIME 2 2026.10 environment used for Stage 1 pins q2-dada2 **2026.10.1** (environment file in `config/`), not 2026.10.0. In 2026.10.1 the R script `run_dada.R` was replaced by a Python module (`q2_dada2/_run_dada.py`, tag `2026.10.1`) that calls the same DADA2 functions. That module was read on 17 November. Every call matched in the table below is the same in 2026.10.1:
> - `filterAndTrim`: truncLen, trimLeft, maxEE and truncQ as set; `rm_phix = TRUE`; maxN at its default 0;
> - `learnErrors`: `nreads = 1e6`, other arguments at DADA2 defaults;
> - `dada`: `pool = FALSE` for independent inference;
> - `mergePairs`: minOverlap, maxMismatch, trimOverhang and returnRejects as set;
> - `removeBimeraDenovo`: consensus, `minFoldParentOverAbundance = 1`, `allowOneOff = FALSE`.
>
> One detail differs from the original note on error learning. In 2026.10.1 the sample order passed to `learnErrors` follows the sample manifest, not an alphabetical file listing. This affects only residual difference 2 below. The matching, the Stage 2 run and its results are unchanged.

**Overrides:** matching settings go in `config/ampliseq/params_matched.yaml` (parameters) and `config/ampliseq/matched.config` (module arguments).

## DADA2 and surrounding steps

| Step | Stage 1 (q2-dada2 2026.10.0) | ampliseq 3.1.0 default | Stage 2 setting | Match |
|---|---|---|---|---|
| Input | Per-sample reads inside `<run>_demux.qza` | Sample sheet of FASTQ files | Same reads exported with checksums (script 21); counts checked against Stage 1 DADA2 input | Identical input |
| Primer removal | None (primers not in reads, D14) | cutadapt | `skip_cutadapt: true` | ✓ |
| Runs | DADA2 per run, tables merged | Per run when a `run` column is given (error learning, denoising, chimera removal) | `run` column | ✓ |
| `filterAndTrim` | truncLen 200/150; trimLeft 0; maxEE 2/2; truncQ 2; rm.phix TRUE; maxN 0 (default) | maxEE `max_ee`; truncQ 2; rm.phix TRUE; maxN 0; minLen 50; maxLen Inf; minQ 0 | `trunclenf 200`, `trunclenr 150`, `max_ee 2`, `max_ee_r 2`, `truncq 2` | ✓ (minLen 50 has no effect: reads are cut to 200/150, and shorter reads are discarded by truncLen) |
| `learnErrors` | `nreads = 1e6`; others at defaults (nbases 1e8, randomize FALSE, MAX_CONSIST 10, OMEGA_C 0); files in `list.files` (alphabetical) order | nbases 1e8, nreads NULL, **randomize TRUE**, seed 100 | Override (`matched.config`): `nreads = 1e6, randomize = FALSE` | ✓ by override. Both pipelines read files alphabetically; the logs state which samples were used, and these are compared |
| `dada` | Independent, sample by sample; DADA2 defaults | **pool = TRUE** (`sample_inference` default changed to pooled in 3.0.0); explicit list of DADA2 defaults | `sample_inference: independent` (pool = FALSE) | ✓ |
| `mergePairs` | minOverlap 12; maxMismatch 0; trimOverhang FALSE; returnRejects FALSE; match/mismatch/gap at defaults (1 / −64 / −64) | Strategy `merge`: minOverlap 12; maxMismatch 0; match 1; mismatch −64; gap −64; trimOverhang FALSE; **returnRejects TRUE** | `mergepairs_strategy: merge` | ✓, **to verify**: rejected pairs are expected to be dropped by `makeSequenceTable`. Smoke test checks that no empty or concatenated sequences appear and that ASV lengths match Stage 1 |
| `removeBimeraDenovo` | consensus; minFoldParentOverAbundance **1**; allowOneOff FALSE; others at defaults | consensus; minFoldParentOverAbundance **2**; minSampleFraction 0.9; ignoreNNegatives 1; minParentAbundance 8; allowOneOff FALSE; minOneOffParentDistance 4; maxShift 16 | Override (`matched.config`): minFoldParentOverAbundance 1, others unchanged | ✓ by override |
| After DADA2 | No filtering | exclude_taxa mitochondria,chloroplast; min_frequency 1; min_samples 1 | `exclude_taxa: none` | ✓ |
| DADA2 version | `bioconductor-dada2` in the rachis-qiime2-2026.10 environment: **to record** | 1.38.0 (module environment) | — | **To record**: a version difference is reported as a known difference |
| Taxonomy used for comparison | SILVA 144 V4 classifier (D24), script 12 settings | DADA2 `assignTaxonomy`, SBDI-GTDB R11-RS232-1 | ampliseq ASVs classified in QIIME 2 with the D24 classifier and identical settings; SBDI-GTDB result kept and reported separately | ✓ |
| Downstream | Project scripts 14–20 | QIIME 2 steps in ampliseq | `skip_qiime: true`, `skip_barrnap: true`; project scripts rerun on the ampliseq table | ✓ |

## Known residual differences

These are reported, not removed:

1. **DADA2 and R versions**, if the two environments differ.
2. **Exact error-learning subset**, if the sample file order differs between the pipelines. This is checked in the logs.
3. **Implementation details** outside the listed calls, such as dereplication and file handling.

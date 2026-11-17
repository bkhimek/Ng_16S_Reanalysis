# Data and samples

## Source

- **Stanford Digital Repository** object [`rm696bw4911`](https://purl.stanford.edu/rm696bw4911), cited in the data availability statement of Ng et al. 2019. Licence CC BY-NC 3.0; 41 files.
- **Library:** 16S rRNA V4 amplicon (515F/806R), Illumina MiSeq paired-end, 6 runs. The humanised experiment uses 3 runs (9 files: forward, index and reverse per run).
- **Download:** script 04 downloads the 9 files. The repository provides no checksums, so integrity rests on byte size, `gzip -t` and equal read counts across each run's three files. Project SHA-256 values are in `metadata/raw_sha256_humanised.txt`.

| Run | Read pairs | Read length | Assigned to humanised samples |
|---|---|---|---|
| hum_run1 | 10,449,045 | 2 × 301 bp | 929,573 |
| hum_run2 | 8,345,334 | 2 × 301 bp | 6,724,830 |
| hum_run3 | 11,178,474 | 2 × 301 bp | 367,220 |
| **Total** | **29,972,853** | | **8,021,623** |

In runs 1 and 3, most unassigned reads carry valid barcodes of samples outside the humanised sheet (other experiments in the same runs).

## Read layout

- Primers are not in the reads, so no trimming is needed. The amplicon is 253 bp, and 2 × 301 bp reads run through it.
- **Barcodes:** 12 bp index read. Barcodes are reused across runs and unique within each run. The sheet's barcodes are reverse complements of Golay codewords.
- **Demultiplexing:** `qiime demux emp-paired` with Golay correction and both reverse-complement options (D17).

## Samples

**Totals:** 453 humanised samples from 20 mice in 5 cages. The full table is `metadata/humanised_sample_metadata.tsv`, built by script 05 from the depositors' sheet. IDs follow D18: samples `H<ID>`, mice `H_m<n>`, cages `H_c<n>`.

**Diet groups:** STD has 10 mice in 3 cages; PD (the paper's MD, MAC-deficient arm) has 10 mice in 2 cages.

**Sequencing runs:**
- STD day-0 faecal samples are in run 1.
- PD day-0 samples and all days 8–13 faecal samples of both diets are in run 2.

**Exclusions:**
- H3106 and H2956: failed libraries (D19). Listed in `metadata/excluded_samples.tsv`.
- H3058 and H3072: same mouse-day label, uncertain identity (D23). Kept in tables, excluded from trajectories.
- Cecal samples are not used in the recovery analysis.

**Depth:** after DADA2 (200/150) every retained sample has ≥ 6,527 reads. All 451 samples pass the 6,400-read rarefaction depth.

## Frozen analysis set

- **Mice and samples:** 17 mice (10 of 10 STD, 7 of 10 PD) and 102 samples (day 0 plus days 8–12). H_m14, H_m16 and H_m17 have no window samples. Day-13 samples of the same mice (17) are used only in a sensitivity analysis.
- **Files:** the list is in `metadata/primary_analysis_samples.tsv` (script 17), and the exact input artifacts with their MD5s are in [frozen_artifacts.tsv](frozen_artifacts.tsv).

| Diet | Cage | Eligible mice |
|---|---|---|
| STD | H_c1 | 4 |
| STD | H_c2 | 4 |
| STD | H_c3 | 2 |
| PD | H_c4 | 4 |
| PD | H_c5 | 3 |

## Stage 2 input

The 451 per-sample FASTQ pairs were exported from the Stage 1 demultiplexed artifacts (script 21). Each sample's read count equals its Stage 1 DADA2 input, and SHA-256 values are in `metadata/per_sample_fastq_sha256.tsv`. The sample sheet for nf-core/ampliseq is `config/ampliseq/samplesheet.tsv`. It contains absolute paths from the original machine, and script 21 regenerates it.

## Known data limits

- **No controls:** no extraction blanks or mock communities were documented in the materials reviewed.
- **Unknown primer version:** the exact primer version (original or modified EMP) is not stated. The classifier was trained on the modified EMP primers (515F Parada, 806R Apprill).
- **Conventional run 2:** 91 of 364 listed faecal samples have no reads in the deposited file (D15). The humanised experiment is unaffected.

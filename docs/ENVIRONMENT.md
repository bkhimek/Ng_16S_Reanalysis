# Software environment

Everything ran on one Windows machine under WSL2 (Ubuntu 24.04): 6 CPUs, 27 GB RAM.

## Stage 1: QIIME 2 (conda)

- **QIIME 2 2026.10** (framework package `rachis`). The official environment file is `config/rachis-qiime2-2026.10-linux-64-conda.yml`. It pins every package, including:
  - q2-dada2 2026.10.1;
  - bioconductor-dada2 1.38.0;
  - R 4.5.2;
  - Python 3.12.16.
- **Build** (strict channel priority fails because the file pins `sortmerna=2.0`, which only bioconda provides):

  ```
  CONDA_CHANNEL_PRIORITY=flexible conda env create --name rachis-qiime2-2026.10 \
      --file config/rachis-qiime2-2026.10-linux-64-conda.yml
  conda activate rachis-qiime2-2026.10
  ```

- **Plugins used:**
  - q2-demux: `emp-paired`;
  - q2-dada2: `denoise-paired`;
  - q2-feature-classifier: `classify-sklearn`;
  - q2-phylogeny: `align-to-tree-mafft-fasttree`;
  - q2-diversity: `alpha-rarefaction`;
  - q2-boots 2026.10.0: `core-metrics`.

## Taxonomy reference: SILVA 144

Script 11 downloads and verifies the classifier, recording it in `metadata/database_manifest.tsv` (D24).

| Item | Value |
|---|---|
| File | `SILVA_144_SSURef_NR99_uniform_classifier_V4-515f-806r.qza` (SILVA directory `current-release/QIIME2/2026.7/SSU/V4-515f-806r/uniform`) |
| Size, MD5 | 79,326,205 bytes; `683a71dcd81ce90d724add22d2ad4d0f` |
| DOI | 10.82364/144/SSU/Ref-NR99/QIIME2/2026.7/V4-515f-806r/uniform |
| Artifact UUID | 76237f8f-87e0-4d1e-93b9-970bb7603cf0 |
| Primers (README) | GTGYCAGCMGCCGCGGTAA / GGACTACNVGGGTWTCTAAT (515F Parada, 806R Apprill) |
| Training | QIIME 2 2026.7.0; prokaryotes; rank propagation; genus level |
| Signature | GPG, SILVA Team key, fingerprint `01877A6A3933ED067EA5B8C3F3A26EDEDB5A4596` (public key in `config/silva/`); verified with `qiime tools signature-verify --name SILVA-Pipeline-Official` |
| Settings | Confidence 0.7, read orientation auto |

## Stage 2: nf-core/ampliseq (Nextflow, Docker)

| Item | Value |
|---|---|
| Pipeline | nf-core/ampliseq 3.1.0, commit `b9bd579c04` (run with `-r 3.1.0`) |
| Nextflow | 26.10.5 (pinned with `NXF_VER` in script 22), on OpenJDK 25.0.2 |
| Container engine | Docker 29.1.3 (`-profile docker`) |
| Settings | `config/ampliseq/params_matched.yaml` (parameters) and `config/ampliseq/matched.config` (module arguments, resource limits 6 CPUs / 24 GB); see [STAGE2_PARAMETER_MATCH.md](STAGE2_PARAMETER_MATCH.md) |
| DADA2 in containers | 1.38.0, with R 4.5.2: the same as Stage 1 |
| Run | `sick_goldstine`, 15 November 2026, 1 h 7 min |

Containers used by the run, taken from `nextflow log <run> -f container`, with their image digests:

| Image | Digest (sha256) |
|---|---|
| community.wave.seqera.io/library/bash_coreutils:25e9437236fbb54f | 4fabf92a1e0d7265d86df5f12e60b2d08dc7b1bcefc5d6f4387a80fb14052e4f |
| community.wave.seqera.io/library/bash_coreutils_findutils_gawk_pruned:63fbf75c749fd019 | 798fce497d77471e7ef9fc64da3fab48c37cc094a58315985990ab2426a804d0 |
| community.wave.seqera.io/library/bioconductor-dada2_r-base_r-digest_tbb:38acac09bac46f36 | 44d1e1927575a4e623759b68820caba070262704a5ac406781957810cf5a8314 |
| community.wave.seqera.io/library/duckdb-cli:1.5.5--9c6d18d9f687a45d | 893bed1639afaa5cf30e0a8d735cc460af1d4755fefc1604b10682350c35cd63 |
| community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc | ae074e961979ea85c72c6c67a8ed432fa3368c736b8fe37799f7223901ac84ab |
| community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-dtplyr_pruned:a6608bc81b0e6546 | 0122fcf542e6237d958035e5e9a758d3b27d4e82f938cc101d46cdf253708b14 |
| quay.io/biocontainers/bioconductor-phyloseq:1.50.0--r44hdfd78af_0 | c9120b7e9faa87cea8c7cccfb13de0eb7d67bfe14f30d40c167e82d006dcbfab |
| quay.io/biocontainers/bioconductor-treesummarizedexperiment:2.10.0--r43hdfd78af_0 | c22b31d28e96f502b87530ebcd8c9f1256d9c7d23ba08a3b2e083ad2ad22a407 |
| quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0 | e194048df39c3145d9b4e0a14f4da20b59d59250465b6f2a9cb698445fd45900 |
| quay.io/biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1 | 44acc2e05f3fd6590635df6725809ba20988e95ae311e4f6d3bf8f20fa5d8e72 |

`results/ampliseq/run_versions.txt`, written by script 22, lists every image on the machine, not only those the run used. The table above is the authoritative record.

## Practical notes

- **Pipes in scripts:** never pipe into `head` (or `grep -m`) in a script with `set -o pipefail`. When `head` closes the pipe early, the script stops silently with exit code 141 (SIGPIPE).
- **Logs:** follow them with `tail -F` (capital F), which waits for the file to appear.
- **Container downloads:** an interrupted download (for example "tls: bad record MAC") is fixed by rerunning; Nextflow's `-resume` keeps finished steps.

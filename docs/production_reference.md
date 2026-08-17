# Production-code reference

This document fixes the provenance of the fourteen scripts supplied by the
author as the code used for the magniFIND-seq production analysis. The original
files contain user- and cluster-specific absolute paths and are therefore not
redistributed. Their SHA-256 values allow a private source archive to be checked
without treating a later reconstruction as the original code.

| File | SHA-256 |
| --- | --- |
| `aggregate_umi_counts.py` | `462e09998b0241aa4a877ea44c8501dd0de4575d6c9f10c0684d7436711a2d3c` |
| `analyze_pipeline.py` | `7b64f55efdfe648d492c75fbea26bdcdcb4c8232009cfbda730270416259decf` |
| `extract_umis.py` | `6d3fb1d385a65767e3a8d6d84fbcf76a3759f54349794a13b347e0058eae9dd0` |
| `map_count_one_sample.sh` | `e8928a9c1b0699dc6c6c3e65f1bce6ed5ad8439962af0671538ad811d97e6b42` |
| `map_fusion_one_sample.sh` | `4694d16963916c37e6bc1b073c40388fe2832540716e786b24283e9e6afff2ea` |
| `preprocess_one_sample.sh` | `e5448b37bb7b54c52fd03713d2479f2bc5ac8b6f9f4a16ee1aa207a4cf78f20d` |
| `preprocessing.py` | `9846765997db92c61c16ffa8225ce1817dc0ef005432b18999218df657053a63` |
| `run_aggregate.sh` | `31096f734a9bbc0f73595aaaf7038b3cabbfe9c0362833eecd1e5df88633a624` |
| `run_analysis.sh` | `402b0b389b8b2be2536c4d4e8f16ef05890d63d10bcda4be43ad1bb295e29917` |
| `run_preprocessing.sh` | `0d6dab2bc00fa4483330437c448d2f57c86b41154ebfa27eb46afb4adcd8b6e8` |
| `star_count.sh` | `55662a97ee0a4b37d1357db8a24e19d1920c2c16124ab948d581040197cd68fc` |
| `submit_fusion_mapping_array.sh` | `67d0fe65b7af52499a09badc9dc1dc4e1cf0a999876da1af5cc62bbae2f3993b` |
| `tag_bam_with_cb_umi.py` | `bf2dd593eaec86c14dff7ceec5bf44de97ba0c8dc26670c127ec4bdb6d1c4223` |
| `umi_methods.py` | `eba2a7cd989fc3f001d45c0eb738f235ec9029edf4b0d87a137e26172b259448` |

## Production-equivalent defaults

- TSO: exact `TTGCGCAATG`; UMI: the next eight bases. Any non-empty
  post-TSO UMI substring was classified as UMI-positive.
- Both UMI and non-UMI FASTQ pairs were retained and trimmed with
  `--illumina --paired --fastqc`.
- Gene STAR attributes: `NH HI AS nM XS`; uniquely mapped alignments were
  selected with MAPQ 255.
- featureCounts: `-F GTF -t exon -g gene_id -R BAM -p -B`; the production
  command did not include `--countReadPairs`.
- UMI-tools used tag extraction, per-cell/per-gene counting, XT/XS/CB/UB tags,
  and did not explicitly add `--method` or `--paired`.
- Fusion STAR mapping used `FUSION_TWOPASS_MODE=None`, no MAPQ filter, and the
  chimeric parameters recorded in `docs/methods.md`.
- Retained production logs identify STAR 2.7.2a and featureCounts 2.0.6.

The supplied scripts do not include the BCR::ABL1 junction-scanning program,
the PBMC/K562-r embedding notebook, or the bulk/single-cell signature script.
Those parts are tied separately to retained notebooks, the archived bulk-DE
table, and figure source-data checks; they must not be described as originating
from these fourteen files.

## Behavioral parity check

A 50,000-pair subset preserving original Illumina headers was processed with
the production UMI implementation and the repository implementation. Both
emitted 48,781 UMI-positive pairs and 1,219 non-UMI pairs. The decompressed
contents were identical for all four outputs (UMI R1/R2 and non-UMI R1/R2).
The repository also contains a fixture for the production short-UMI,
empty-post-UMI, and no-TSO behavior.

## Intentional correctness hardening

The public implementation preserves scientific behavior while replacing
cluster-specific orchestration with explicit manifests and configuration
files. It also validates paired FASTQ structure, fails on missing CB/UB tags,
uses content fingerprints for restart decisions, writes deterministic gzip
files, and refuses incomplete sample aggregation. These checks do not alter a
valid production input but turn silent truncation, sample discovery, or stale
output reuse into explicit errors.

Two production edge cases are retained transparently. An all-UMI sample has an
empty non-UMI split; the archived shell wrapper attempted to run Trim Galore on
that empty pair, whereas this repository emits a valid empty trimmed pair so
the sample can continue. The archived read-tracking script counted both raw R1
and R2 despite labeling that metric as reads; raw QC therefore writes both
`readtracking_production_legacy.tsv` and the R1-only
`readtracking_corrected_r1.tsv` companion.

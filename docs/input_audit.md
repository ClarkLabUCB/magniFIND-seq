# Local study-input audit

This report records read-only validation performed on source data supplied
locally during repository preparation. The source archive remains ignored by
Git and is not redistributed by this repository.

## Bulk K562 sensitive/resistant archive

Audit date: 2026-08-05. Local archive name: `SPGZKS_fastq.zip`.

- archive size: 4,954,034,377 bytes;
- archive SHA-256:
  `945ee4d88bd5d0f736c3f9079f48026e03a2e6240650dcaa37fb6fc0a4702869`;
- ZIP test: all six stored members passed;
- inner gzip test: all six streams passed CRC/truncation validation; and
- full FASTQ scan: every record had an `@` header, `+` separator, complete
  four-line structure, and equal sequence/quality lengths.

| Library | Condition | Single-end reads |
| --- | --- | ---: |
| `SPGZKS_1_K562-r_1.fastq.gz` | K562-R | 20,633,640 |
| `SPGZKS_2_K562-r_2.fastq.gz` | K562-R | 20,073,705 |
| `SPGZKS_3_K562-r_3.fastq.gz` | K562-R | 17,973,229 |
| `SPGZKS_4_K562-s_1.fastq.gz` | K562-S | 20,340,978 |
| `SPGZKS_5_K562-s_2.fastq.gz` | K562-S | 19,437,999 |
| `SPGZKS_6_K562-s_3.fastq.gz` | K562-S | 19,214,493 |

The corresponding two-column condition file is tracked as
`03_downstream_analysis/config/study_bulk_samples.tsv`. No reduced FASTQ
fixture was created: a downsample cannot validate `filterByExpr`, dispersion,
log2FC, or FDR rankings from the complete three-versus-three experiment.

The retained production report records fastp 0.24.0, STAR 2.7.11, samtools
1.22.1, UMIcollapse 1.1.0, featureCounts 2.1.1, and edgeR 4.0.16. featureCounts
used strand-specific fractional assignment to exon and three-prime-UTR features
grouped by `gene_id`. Exact reference FASTA/GTF identities and the full command
line still require archival checksums. The archived production differential-
expression table is accepted directly by the concordance workflow so that a
different edgeR release is not silently substituted for manuscript ranking.

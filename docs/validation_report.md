# Validation report

Validation date: 2026-08-06.

## Passed locally

- After consolidation into the three numbered analysis areas, `make lint`,
  `make test-unit`, `make test-preprocessing`, `make test-fusion`,
  `make test-downstream`, and `make test-study-fastq` passed from the new
  entry points. The guided notebook also passed structural and Python-cell
  validation.
- Shell syntax, Python byte compilation, R parsing, YAML/CFF parsing, local
  Markdown links, line endings, trailing whitespace, and executable-bit audit.
- Manifest validation, malformed-input rejection, deterministic subset
  creation, byte-deterministic gzip count aggregation, cross-sample UMI-table
  rejection, multi-sample count aggregation, and non-integer-count rejection.
- Seurat object construction after `min.cells` feature filtering, including
  exact alignment of retained assay feature metadata.
- Targeted fusion fixtures including both STAR report orders, strand rejection,
  empty negatives, missing/incomplete states, QNAME/UB/junction counts, and the
  manuscript one-read binary candidate definition.
- Public K562 fixtures including SAM flag `0x900`, reverse-read restoration,
  exact 96-base FASTQ, GTF gene-body overlaps, repeated-QNAME deduplication,
  four fusion categories, called-cell CB-UMI denominators, and percentages.
- The GSM2406675 processed barcode list was downloaded from GEO, matched its
  pinned SHA-256, passed gzip validation, and contained 5,768 valid barcodes.
- Raw-QC fixture including inclusive threshold boundaries, dual-filter status,
  R-compatible type-7 65th/70th percentiles, and the archived raw R1+R2 count
  alongside its corrected R1-only companion.
- edgeR/Seurat concordance fixture in R 4.3.3, edgeR 3.42.4, Matrix 1.6-5,
  Seurat 5.0.3, ggplot2 3.5.2, and ggrepel 0.9.6. This covered FindMarkers
  effect sizes, full-signature enforcement, replicate-count enforcement,
  50-per-direction selection logic, PBMC score-only retention,
  figures, output replacement, and restart fingerprinting.
- Joint PBMC/K562-r and independently recomputed K562-r-only UMAP fixture,
  including explicit assay selection, integer-count validation, missing-K562
  candidate rejection, and byte-identical coordinate tables across repeat runs.
- Both study-derived K562-r examples: 100,000 paired records passed full FASTQ
  structure/pairing checks; exact TSO+UMI counts matched provenance; regeneration
  from the full source library remained byte-identical.
- The six full bulk libraries passed ZIP, gzip CRC, and every-record FASTQ
  structure/sequence-quality validation. Counts and archive SHA-256 are in
  `docs/input_audit.md`.
- The explicit Linux lock pins the production STAR 2.7.2a and is checked by the
  Linux GitHub Actions integration run.
- A fresh osx-64 environment created from `environment.yml` completed the full
  `make test` suite locally. The four-pair integration fixture ran STAR 2.7.2a,
  featureCounts 2.0.6, UMI-tools, Seurat, and the targeted fusion workflow from
  FASTQ inputs to exact expected outputs; a second run correctly skipped every
  completed stage.
- Production UMI extraction was compared against the archived implementation
  on the first 50,000 pairs of an author-provided K562-r library. UMI-positive
  R1/R2 (48,781 pairs) and non-UMI R1/R2 (1,219 pairs) were decompressed-byte
  identical in all four outputs. A separate edge fixture pins short UMI,
  empty-post-UMI, and no-TSO behavior.

## Environment-limited checks

The portable four-pair FASTQ integration test passes locally and remains
configured in GitHub Actions on `ubuntu-latest`. Linux CI is still the release
gate because the explicit lock file targets Linux and the production computing
environment was Linux/HPC.

A full GRCh38 STAR index was not built locally because the machine has 16 GB of
RAM. Consequently the 50,000-pair study examples were validated and regenerated
but not aligned locally against the full human reference.

## Data-limited checks

Final manuscript values cannot yet be regenerated because the study GEO/SRA
accession, full archived processed single-cell objects, original tag-preserving
SRR5082088 BAM, and reference checksums are not all available. The production
bulk method report and differential-expression table are locally available, and
the repository accepts the latter without recomputing its ranks. The processed
5,768-entry called-cell barcode list is available and
checksum-pinned, but it does not replace the original tagged BAM. These are
explicit unresolved inputs rather than silent test skips; see
`docs/information_to_complete.md`.

# Examples

Examples are user-facing walkthroughs. They complement, rather than replace,
the automated regression tests under `tests/`.

## Minimal FASTQ example

[`minimal_fastq/`](minimal_fastq/) reuses the repository's four-pair synthetic
FASTQ fixture and miniature reference. It builds a STAR index, runs
gene-expression quantification and fusion detection, verifies the expected
counts and fusion call, and preserves all outputs for inspection.

No study FASTQ or human sequence is used by this minimal example. It is
intentionally small enough to run on a workstation or in continuous
integration.

## K562-r single-bead example

[`k562r_single_bead/`](k562r_single_bead/) contains de-identified,
deterministically selected subsets of two author-provided K562-r single-bead
libraries. It provides realistic TSO/UMI structure and an optional complete
mapping workflow using a separately downloaded GENCODE v44 GRCh38 reference.
There is no BCR::ABL1-negative control, so its fusion output is exploratory and
is not a sensitivity or specificity test.

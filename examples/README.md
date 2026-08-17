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

## Extending the example with study-like data

A biologically representative example remains separate from the synthetic
regression fixture. Additional inputs should be de-identified,
publication-approved paired FASTQs generated with the same library structure,
preferably including:

- at least one expected BCR::ABL1-positive and one negative sample;
- intact read-1 TSO plus 8-base UMI structure and the corresponding read 2;
- a sample manifest containing non-identifying sample names; and
- the exact reference assembly and GTF used for the study.

Do not commit production FASTQs, BAMs, sample identifiers, or protected human
metadata to this repository. Reference external files through a manifest or a
GEO/SRA accession when one becomes available.

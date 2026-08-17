# K562-R junction scanner fixture

This compact fixture tests parsing and summarization of precomputed STAR
`Chimeric.out.junction` records. It contains coordinate/CIGAR fields and UMI
labels for two candidate-positive cases plus an empty negative case. It contains
no nucleotide sequences, FASTQ files, BAM files, or instrument identifiers.

The junction patterns were reduced from K562-R analysis output. Query names
were replaced with deterministic `fixture_read_<n>` identifiers before
inclusion. Sample names are generic fixture labels. These records are used only
to exercise parsing, report-order normalization, strand/interval filtering, UMI
summarization, and exact expected-output comparison.

This fixture does not test STAR mapping. The independent, fully synthetic FASTQ
fixture under `tests/fastq_integration/` tests mapping through both workflows.

# K562-r single-bead FASTQ example

This example uses de-identified subsets of two author-provided K562-r
single-bead RNA-sequencing libraries. It demonstrates the library structure and
provides realistic inputs for gene-expression quantification and exploratory
BCR::ABL1 fusion screening.

Both samples are K562-r single beads. There is no BCR::ABL1-negative biological
control in this fixture. A `not_detected` result means only that no qualifying
fusion-transcript read was observed in the downsampled reads; it must not be
interpreted as a BCR::ABL1-negative cell or as an estimate of assay sensitivity.

## Included data

Each sample contains 50,000 paired 100-bp reads selected across the complete
source FASTQ by the smallest seeded BLAKE2b hashes of paired read identifiers.
Original instrument headers were replaced with deterministic identifiers.
Sequences and qualities were otherwise retained unchanged.

| Sample | Source pairs | Included pairs | Included complete TSO+UMI pairs |
| --- | ---: | ---: | ---: |
| `k562r_single_bead_01` | 1,898,980 | 50,000 | 48,282 |
| `k562r_single_bead_02` | 2,193,830 | 50,000 | 48,216 |

Content hashes and selection parameters are recorded in `data/*.provenance.tsv`.
The subset can be regenerated from an authorized source copy with
[`scripts/create_deidentified_fastq_subset.py`](../../scripts/create_deidentified_fastq_subset.py).

## Input-only validation

The lightweight study-fixture test requires no genome reference:

```bash
make test-study-fastq
```

It verifies compressed-file integrity, paired records, de-identified names,
record counts, SHA-256 hashes, and exact TSO/UMI-positive counts.

## Reference preparation

The runnable mapping example uses the GENCODE release 44 GRCh38 primary
assembly FASTA and matching comprehensive GTF. These source files and the STAR
index are intentionally not committed because they require substantial disk
space. On a Linux host with at least 32 GB RAM and the pinned environment:

```bash
bash examples/k562r_single_bead/prepare_reference.sh
```

The script downloads both files from GENCODE, verifies their published MD5
checksums, and creates a STAR 2.7.2a index under
`reference/GRCh38_GENCODEv44/`. This is a documented example reference; the
exact reference release used for the final manuscript analysis should be
confirmed before claiming exact paper reproduction.

## Complete example

After reference preparation:

```bash
bash examples/k562r_single_bead/run_example.sh
```

Outputs are retained under `examples/k562r_single_bead/output/`. Full-genome
mapping is not part of CI because reference construction is resource-intensive.
The synthetic miniature example remains the portable end-to-end CI test.

## Distribution scope

The included reads are derived from a K562 cell-line experiment. 
Public release of derived study reads should remain consistent with
the manuscript's final data-sharing plan. The repository's MIT license covers
the software; the authors should confirm the intended reuse terms for these
study-derived FASTQ subsets before the public release.

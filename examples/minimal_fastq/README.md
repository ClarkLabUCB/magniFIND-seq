# Minimal FASTQ walkthrough

This runnable example shows the expected user workflow and leaves its outputs
in place for inspection. It uses the same fully synthetic input described in
[`tests/fastq_integration/README.md`](../../tests/fastq_integration/README.md):
four paired reads, two mock genes, three expected gene-level molecules, and one
mock BCR::ABL1-like chimeric alignment.

The fixture verifies software behavior only. It is not biological evidence and
cannot reproduce the manuscript's quantitative results.

## Requirements

Use Linux with the pinned root environment activated:

```bash
conda env create -f environment.yml
conda activate magnifind-seq-analysis
```

The example requires `STAR`, `trim_galore`, `samtools`, `featureCounts`,
`umi_tools`, `Rscript`, `python3`, and the R packages Matrix and Seurat. The
current Bioconda macOS STAR build can be unreliable under Apple Silicon/Rosetta;
Linux is the supported platform for this complete example.

## Run

From the repository root:

```bash
make example
```

or directly:

```bash
bash examples/minimal_fastq/run_example.sh
```

The first run builds a small STAR index and then runs both workflows. Later
runs reuse the index and valid stage fingerprints.

## Expected result

The final verification reports:

```text
Gene-level molecules
  ABL1_mock: 1
  BCR_mock: 2
Fusion screen
  sample: k562_r_fastq_fixture
  unique UMI: 1
  call: candidate_ge_min_umi
```

Inspect these main outputs:

```text
output/
├── gene_expression/
│   ├── gene_by_sample_counts.tsv.gz
│   ├── matrix.mtx
│   ├── features.tsv
│   ├── barcodes.tsv
│   ├── seurat/minimal_fastq_example.seurat.rds
│   └── run_metadata/
├── fusion_detection/
│   ├── fusion_mapped/
│   ├── fusion_scan/
│   └── run_metadata/
└── reference/
    └── star_index/
```

The `output/` directory is ignored by Git. To start from scratch, remove or
move this example's `output/` directory and rerun the command.

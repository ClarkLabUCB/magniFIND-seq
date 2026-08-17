# 01 — Preprocessing

This area converts paired magniFIND-seq FASTQs into per-sample gene counts, a
combined sparse UMI matrix, and a Seurat object. It also contains the
experiment-specific raw read/UMI QC summary used for the K562-in-EL4 analysis.

## Inputs and outputs

- Inputs: paired FASTQs, a three-column `sample/r1/r2` manifest, STAR index,
  matching GTF, and optional sample metadata.
- Gene outputs: UMI and non-UMI extracted/trimmed FASTQs, alignments, UMI tables, sparse matrix,
  feature/barcode tables, and Seurat RDS.
- QC outputs: per-sample assigned-read and deduplicated-UMI table, retained and
  removed sample lists, the production-legacy read-tracking table, and a
  corrected R1-only companion.

## Run

```bash
cp config/gene_expression.example.sh config/gene_expression.sh
cp config/samples.example.tsv config/samples.tsv
# Edit the copies.
bash code/run_gene_expression.sh config/gene_expression.sh

# Optional manuscript raw-QC summary after gene quantification:
cp config/raw_qc.example.sh config/raw_qc.sh
bash code/run_raw_qc.sh config/raw_qc.sh
```

The raw-QC cutoffs are experiment-specific; do not apply them to a different
dataset without a documented rationale. Detailed algorithms and parameters are
specified in [`../docs/methods.md`](../docs/methods.md).

## Test

From the repository root:

```bash
make test-preprocessing
```

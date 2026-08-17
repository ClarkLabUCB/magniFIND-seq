# 02 — Fusion analysis

This area contains two related workflows: targeted STAR remapping and
BCR::ABL1 evidence screening for magniFIND-seq libraries, and intergenic
chimeric-transcript quantification from the public K562 dataset SRR5082088.

## Inputs and outputs

- Targeted inputs: UMI-trimmed paired FASTQs from preprocessing and a matching
  STAR index.
- Targeted outputs: unfiltered BAMs, `Chimeric.out.junction` files, CB/UB-tagged
  BAMs, and per-sample BCR::ABL1 evidence/candidate tables.
- Public K562 inputs: a tag-preserving BAM, pinned called-cell barcodes, GTF,
  and compatible STAR index.
- Public output: intergenic, BCR, ABL1, and BCR::ABL1 read/UMI/cell summaries.

## Run

```bash
cp config/targeted.example.sh config/targeted.sh
cp config/targeted_samples.example.tsv config/targeted_samples.tsv
# Edit both copies.
bash code/run_targeted_mapping.sh config/targeted.sh
bash code/run_targeted_scan.sh config/targeted.sh

# Optional public K562 analysis:
cp config/public_k562.example.sh config/public_k562.sh
bash code/prepare_called_barcodes.sh
bash code/run_public_k562.sh config/public_k562.sh
```

Configured coordinates must match the reference assembly. Candidate calls are
screening evidence and require orthogonal validation. See
[`../docs/methods.md`](../docs/methods.md).

## Test

From the repository root:

```bash
make test-fusion
```

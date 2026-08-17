# magniFIND-seq analysis pipeline

Reproducible analysis code accompanying the magniFIND-seq manuscript,
*Single-cell DNA cytometry with magnetic- and fluorescence-activated bead
sorting*. The public-facing repository is intentionally organized into three
analysis areas:

| Area | What it does | Main result |
| --- | --- | --- |
| [`01_preprocessing/`](01_preprocessing/) | TSO/UMI extraction, gene quantification, Seurat construction, raw QC | UMI matrix, Seurat object, QC table |
| [`02_fusion_analysis/`](02_fusion_analysis/) | Targeted BCR::ABL1 screening and public K562 chimera analysis | Fusion evidence tables and unfiltered mapping files |
| [`03_downstream_analysis/`](03_downstream_analysis/) | PBMC/K562-r embeddings and bulk/single-cell concordance | UMAPs, 100-gene signature, concordance and resistance-score results |

Shared examples, tests, documentation, and provenance helpers remain at the
repository root so the scientific workflows stay easy to find without hiding
the validation machinery.

## Installation

Linux is the supported production platform. Create the complete environment:

```bash
conda env create -f environment.yml
conda activate magnifind-seq-analysis
```

For an exact Linux package recreation, use the explicit lock file:

```bash
conda create --name magnifind-seq-analysis --file environment-linux-64.lock
```

The same genome assembly and annotation must be used for the FASTA, GTF, STAR
index, and configured fusion coordinates.

## Run the three analysis areas

Each area has one short README, example configurations, executable code, and
its own fixture tests. Configuration paths are resolved relative to the config
file; FASTQ paths in a manifest are resolved relative to that manifest.

```bash
# 1. Preprocessing
cp 01_preprocessing/config/gene_expression.example.sh \
  01_preprocessing/config/gene_expression.sh
cp 01_preprocessing/config/samples.example.tsv \
  01_preprocessing/config/samples.tsv
# Edit both files, then:
bash 01_preprocessing/code/run_gene_expression.sh \
  01_preprocessing/config/gene_expression.sh

# 2. Targeted fusion mapping and BCR::ABL1 scan
cp 02_fusion_analysis/config/targeted.example.sh \
  02_fusion_analysis/config/targeted.sh
cp 02_fusion_analysis/config/targeted_samples.example.tsv \
  02_fusion_analysis/config/targeted_samples.tsv
# Edit both files, then:
bash 02_fusion_analysis/code/run_targeted_mapping.sh \
  02_fusion_analysis/config/targeted.sh
bash 02_fusion_analysis/code/run_targeted_scan.sh \
  02_fusion_analysis/config/targeted.sh

# 3. Downstream bulk/single-cell concordance
cp 03_downstream_analysis/config/concordance.example.sh \
  03_downstream_analysis/config/concordance.sh
bash 03_downstream_analysis/code/run_concordance.sh \
  03_downstream_analysis/config/concordance.sh
```

To execute any configured combination in sequence, copy
`03_downstream_analysis/config/all.example.sh`, fill in its config paths, and
run `make run-all CONFIG=03_downstream_analysis/config/all.sh`.

Configuration files are sourced by Bash and must come from a trusted source.

## Try it and validate it

The synthetic example exercises preprocessing and targeted fusion analysis
without manuscript data:

```bash
make example
```

The study-derived test contains de-identified 50,000-read-pair subsets from two
K562-r single-bead libraries and checks input integrity plus exact TSO/UMI
statistics without requiring a genome reference:

```bash
make test-study-fastq
```

Run validation by scientific area, or run everything:

```bash
make test-preprocessing
make test-fusion
make test-downstream
make test
```

Automated regression and integration checks are implemented under
[`tests/`](tests/) and the three analysis-area test directories. For a readable
guided validation interface, open
[`validation/validation_walkthrough.ipynb`](validation/validation_walkthrough.ipynb).
It invokes the same Make targets and reports pass/fail status and logs; it does
not contain a second implementation of the pipeline.

## Reproducibility boundaries

Completed workflows retain normalized manifests, configuration snapshots,
effective parameters, tool versions, SHA-256 input/source fingerprints, and
source-control state. Outputs are assembled in temporary paths before they
replace a prior result. Fusion calls are computational screening results, not
clinical diagnoses. Bulk/single-cell correlations and the resistance module
score are descriptive and may be confounded by study design.

The included fixtures verify code operation but cannot reproduce all manuscript
figures without the complete production inputs. GEO/SRA identifiers and other
facts unavailable before submission are recorded only in
[`docs/information_to_complete.md`](docs/information_to_complete.md), not guessed
in code or documentation.

Production-equivalent defaults are derived from the author-supplied scripts and
retained run logs. Their immutable checksums and the recovered commands are in
[`docs/production_reference.md`](docs/production_reference.md). A configurable
option is not described as production-equivalent unless that reference supports
it.

See [`docs/methods.md`](docs/methods.md) for the analysis specification,
[`docs/validation_report.md`](docs/validation_report.md) for the current test
scope, [`docs/production_reference.md`](docs/production_reference.md) for the
source-code audit, and [`docs/repository_layout.md`](docs/repository_layout.md)
for the full file map.

## Citation and license

Citation metadata is in [`CITATION.cff`](CITATION.cff). Software is distributed
under the [`LICENSE`](LICENSE).

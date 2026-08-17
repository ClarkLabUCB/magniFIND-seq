# Repository layout

```text
magnifind-seq_analysis_pipeline/
├── 01_preprocessing/
│   ├── README.md
│   ├── code/
│   │   ├── run_gene_expression.sh
│   │   ├── run_raw_qc.sh
│   │   ├── gene_expression/
│   │   └── raw_qc/
│   ├── config/
│   └── tests/
├── 02_fusion_analysis/
│   ├── README.md
│   ├── code/
│   │   ├── run_targeted_mapping.sh
│   │   ├── run_targeted_scan.sh
│   │   ├── run_public_k562.sh
│   │   ├── targeted/
│   │   └── public_k562/
│   ├── config/
│   └── tests/
├── 03_downstream_analysis/
│   ├── README.md
│   ├── code/
│   │   ├── run_embedding.sh
│   │   ├── run_concordance.sh
│   │   ├── run_all.sh
│   │   ├── embedding/
│   │   └── concordance/
│   ├── config/
│   └── tests/
├── examples/
│   ├── minimal_fastq/
│   └── k562r_single_bead/
├── validation/
│   └── validation_walkthrough.ipynb
├── tests/
│   ├── README.md
│   ├── run_lightweight_tests.sh
│   ├── run_fastq_integration_test.sh
│   ├── run_study_fastq_fixture_test.sh
│   └── fastq_integration/
├── scripts/
├── docs/
├── .github/workflows/ci.yml
├── environment.yml
├── environment-linux-64.lock
├── Makefile
└── README.md
```

The three numbered directories are the scientific entry points. Each contains
its own explanation, configuration templates, executable code, and focused
fixtures. Root-level `tests/` checks behavior spanning more than one area;
`validation/` is the optional human-facing dashboard for those same Make
targets; `scripts/` contains shared provenance and validation helpers. One root
environment covers all workflows. Generated outputs and production inputs are
ignored by Git, apart from explicitly documented test fixtures.

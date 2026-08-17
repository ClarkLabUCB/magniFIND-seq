PYTHON ?= python3

.PHONY: help lint notebook-check example run-all test test-unit test-preprocessing test-fusion test-downstream test-fastq test-study-fastq

help:
	@printf '%s\n' \
	  'make test              Run the complete validation suite' \
	  'make test-preprocessing Test gene-expression and raw-QC fixtures' \
	  'make test-fusion       Test targeted and public K562 fusion analyses' \
	  'make test-downstream   Test embedding and bulk/single-cell concordance' \
	  'make test-fastq        Run the full synthetic FASTQ integration test' \
	  'make test-study-fastq  Validate study-derived FASTQ subsets' \
	  'make example           Run the persistent minimal FASTQ example' \
	  'make run-all CONFIG=path  Run configured study workflows' \
	  'make notebook-check    Validate the guided Jupyter notebook' \
	  'make lint              Syntax-check all source files'

notebook-check:
	$(PYTHON) scripts/validate_notebook.py validation/validation_walkthrough.ipynb

lint: notebook-check
	@find . -type f -name '*.sh' -print0 | xargs -0 -n1 bash -n
	@PYTHONPYCACHEPREFIX=/tmp/magnifind-seq-pycache $(PYTHON) -m compileall -q \
	  scripts 01_preprocessing/code 02_fusion_analysis/code \
	  tests examples/minimal_fastq
	@find . -type f -name '*.R' -print0 | xargs -0 -n1 Rscript -e \
	  'invisible(parse(file=commandArgs(trailingOnly=TRUE)[1]))'

test-unit:
	bash tests/run_lightweight_tests.sh

test-preprocessing:
	bash 01_preprocessing/tests/run_seurat_creation_test.sh
	bash 01_preprocessing/tests/run_raw_qc_test.sh

test-fusion:
	bash 02_fusion_analysis/tests/targeted/run_fusion_detection_fixture_test.sh
	bash 02_fusion_analysis/tests/public_k562/run_test.sh

test-downstream:
	bash 03_downstream_analysis/tests/concordance/run_concordance_fixture_test.sh
	bash 03_downstream_analysis/tests/embedding/run_embedding_fixture_test.sh

test-fastq:
	bash tests/run_fastq_integration_test.sh

test-study-fastq:
	bash tests/run_study_fastq_fixture_test.sh

example:
	bash examples/minimal_fastq/run_example.sh

run-all:
	@test -n "$(CONFIG)" || { echo 'Usage: make run-all CONFIG=03_downstream_analysis/config/all.sh' >&2; exit 2; }
	bash 03_downstream_analysis/code/run_all.sh "$(CONFIG)"

test: lint test-unit test-preprocessing test-fusion test-fastq test-study-fastq test-downstream

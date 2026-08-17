# Automated tests

This directory contains repository-level regression and integration tests.
They are non-interactive, return a nonzero exit status on failure, and are run
by both `make test` and GitHub Actions.

- `run_lightweight_tests.sh` checks manifests, malformed-input rejection,
  production UMI edge behavior, recovered command defaults, and deterministic
  count aggregation.
- `run_fastq_integration_test.sh` builds a miniature STAR index and tests the
  complete FASTQ-to-Seurat and targeted-fusion workflows against exact expected
  outputs.
- `run_study_fastq_fixture_test.sh` verifies the integrity, pairing,
  de-identification, and TSO/UMI structure of the two study-derived K562-r
  subsets.
- `fastq_integration/` stores the synthetic integration fixture and expected
  outputs.

Focused tests that belong to one scientific area remain next to that area in
`01_preprocessing/tests/`, `02_fusion_analysis/tests/`, and
`03_downstream_analysis/tests/`.

The notebook under [`../validation/`](../validation/) is only a guided front
end for these Make targets; it does not contain a separate test implementation.

# Concordance fixture

The fixture is fully synthetic and tests software behavior only. It contains
three bulk replicates per condition, four single cells per condition, and two
additional PBMC-labeled score-only cells. Two
protein-coding genes are designed to increase in each direction. Strong
mitochondrial- and ribosomal-prefix controls verify the configured optional
exclusion regular expression. An intentionally impossible DEG fold-change
threshold verifies that manuscript-mode signature selection is independent of
the separate DEG table.

The test checks assay feature-metadata matching, expected selected genes,
comparison direction, concordance metrics, resistance-score direction, output
figures and provenance, and fingerprint-based restart behavior.

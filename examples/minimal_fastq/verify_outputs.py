#!/usr/bin/env python3
"""Verify and summarize outputs from the runnable minimal FASTQ example."""

import csv
import gzip
from pathlib import Path


EXAMPLE_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXAMPLE_DIR.parent.parent
EXPECTED_DIR = REPO_ROOT / "tests" / "fastq_integration" / "expected"
GENE_OUT = EXAMPLE_DIR / "output" / "gene_expression"
FUSION_OUT = EXAMPLE_DIR / "output" / "fusion_detection"
SAMPLE = "k562_r_fastq_fixture"


def read_expected_counts():
    with (EXPECTED_DIR / "gene_counts.tsv").open(encoding="utf-8") as handle:
        return {
            row["gene"]: int(row["count"])
            for row in csv.DictReader(handle, delimiter="\t")
        }


def read_actual_counts():
    path = GENE_OUT / "gene_by_sample_counts.tsv.gz"
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or SAMPLE not in reader.fieldnames:
            raise SystemExit(f"Missing sample column {SAMPLE!r} in {path}")
        gene_column = reader.fieldnames[0]
        return {row[gene_column]: int(row[SAMPLE]) for row in reader}


def read_fusion_summary():
    path = (
        FUSION_OUT
        / "fusion_scan"
        / "minimal_fastq_example_BCR_ABL1_expected_orientation_summary_all_samples.tsv"
    )
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise SystemExit(f"Expected one fusion-summary row in {path}, found {len(rows)}")
    return rows[0], path


def main():
    expected_counts = read_expected_counts()
    actual_counts = read_actual_counts()
    if actual_counts != expected_counts:
        raise SystemExit(
            f"Gene counts differ from the fixture expectation: "
            f"{actual_counts} != {expected_counts}"
        )

    expected_features = (EXPECTED_DIR / "features.tsv").read_text(encoding="utf-8")
    actual_features = (GENE_OUT / "features.tsv").read_text(encoding="utf-8")
    if actual_features != expected_features:
        raise SystemExit("features.tsv differs from the fixture expectation")
    if (GENE_OUT / "barcodes.tsv").read_text(encoding="utf-8") != f"{SAMPLE}\n":
        raise SystemExit("barcodes.tsv differs from the fixture expectation")

    fusion, fusion_path = read_fusion_summary()
    expected_call = "candidate_ge_min_umi"
    if fusion["sample"] != SAMPLE or fusion["expected_BCR_ABL1_candidate_call"] != expected_call:
        raise SystemExit(f"Unexpected fusion summary: {fusion}")
    if int(fusion["unique_UB"]) != 1:
        raise SystemExit(f"Expected one fusion-supporting UMI, found {fusion['unique_UB']}")

    seurat_path = GENE_OUT / "seurat" / "minimal_fastq_example.seurat.rds"
    if not seurat_path.is_file() or seurat_path.stat().st_size == 0:
        raise SystemExit(f"Missing Seurat output: {seurat_path}")

    print("[PASS] Minimal FASTQ example matched its exact expected result")
    print("Gene-level molecules")
    for gene, count in sorted(actual_counts.items()):
        print(f"  {gene}: {count}")
    print("Fusion screen")
    print(f"  sample: {fusion['sample']}")
    print(f"  unique UMI: {fusion['unique_UB']}")
    print(f"  call: {fusion['expected_BCR_ABL1_candidate_call']}")
    print("Main outputs")
    print(f"  counts: {GENE_OUT / 'gene_by_sample_counts.tsv.gz'}")
    print(f"  Seurat: {seurat_path}")
    print(f"  fusion summary: {fusion_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import csv
import gzip
import os

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


def parse_gtf_attributes(attr_field):
    attrs = {}
    for item in attr_field.strip().strip(";").split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def load_gtf_gene_info(gtf_path):
    gene_info = {}
    with open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue
            gene_name = attrs.get("gene_name", gene_id)
            gene_type = attrs.get("gene_biotype", attrs.get("gene_type", "NA"))
            annotation = (gene_name, gene_type)
            if gene_id not in gene_info or fields[2] == "gene":
                gene_info[gene_id] = annotation
    return gene_info


def load_manifest_samples(manifest_path):
    with open(manifest_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["sample", "r1", "r2"]:
            raise ValueError(
                "Manifest header must contain exactly: sample, r1, r2"
            )
        samples = [row["sample"].strip() for row in reader if row.get("sample", "").strip()]
    if not samples:
        raise ValueError(f"No samples found in manifest: {manifest_path}")
    if len(samples) != len(set(samples)):
        raise ValueError(f"Duplicate samples found in manifest: {manifest_path}")
    return samples


def load_counts(counts_root, expected_samples):
    count_tables = {}
    for sample in expected_samples:
        path = os.path.join(counts_root, sample, f"{sample}_counts.tsv.gz")
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Missing UMI count table for {sample}: {path}")
        with gzip.open(path, "rt") as handle:
            df = pd.read_csv(handle, sep="\t", usecols=["gene", "cell", "count"])
        invalid_cell = df["cell"].isna() | df["cell"].astype(str).str.strip().ne(sample)
        if invalid_cell.any():
            raise ValueError(
                f"UMI count table contains a cell ID other than {sample!r}: {path}"
            )
        invalid_gene = df["gene"].isna() | df["gene"].astype(str).str.strip().eq("")
        if invalid_gene.any():
            raise ValueError(f"Missing gene identifier in UMI count table: {path}")
        numeric_counts = pd.to_numeric(df["count"], errors="coerce")
        if (
            numeric_counts.isna().any()
            or not np.isfinite(numeric_counts).all()
            or (numeric_counts < 0).any()
            or (numeric_counts != np.floor(numeric_counts)).any()
        ):
            raise ValueError(
                f"UMI counts must be finite, non-negative integers: {path}"
            )
        df["count"] = numeric_counts.astype(np.int64)
        df = df.groupby("gene", sort=False)["count"].sum().to_frame(sample)
        count_tables[sample] = df
    return count_tables


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate umi_tools per-sample count tables into Seurat-compatible files."
    )
    parser.add_argument("--counts-root", required=True, help="Directory containing *_counts.tsv.gz")
    parser.add_argument("--manifest", required=True, help="Normalized sample manifest")
    parser.add_argument("--gtf", required=True, help="GTF used for featureCounts")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    expected_samples = load_manifest_samples(args.manifest)
    count_tables = load_counts(args.counts_root, expected_samples)
    if not count_tables:
        raise FileNotFoundError(f"No *_counts.tsv.gz files found under {args.counts_root}")

    merged = pd.concat(count_tables.values(), axis=1).fillna(0).astype(int)
    merged = merged.reindex(columns=expected_samples)
    merged.sort_index(inplace=True)
    if merged.shape[0] == 0 or int(merged.to_numpy().sum()) == 0:
        raise ValueError("No positive gene-level UMI counts were available to aggregate")

    gene_info = load_gtf_gene_info(args.gtf)
    features_path = os.path.join(args.out_dir, "features.tsv")
    barcodes_path = os.path.join(args.out_dir, "barcodes.tsv")
    matrix_path = os.path.join(args.out_dir, "matrix.mtx")

    mmwrite(matrix_path, sparse.csr_matrix(merged.values))

    with open(features_path, "wt") as handle:
        for gene_id in merged.index:
            gene_name, gene_type = gene_info.get(gene_id, (gene_id, "NA"))
            handle.write(f"{gene_id}\t{gene_name}\t{gene_type}\n")

    with open(barcodes_path, "wt") as handle:
        for sample in merged.columns:
            handle.write(f"{sample}\n")

    merged.to_csv(
        os.path.join(args.out_dir, "gene_by_sample_counts.tsv.gz"),
        sep="\t",
        index_label="gene_id",
        compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
    )
    print(f"[INFO] Wrote {matrix_path}")
    print(f"[INFO] Wrote {features_path}")
    print(f"[INFO] Wrote {barcodes_path}")


if __name__ == "__main__":
    main()

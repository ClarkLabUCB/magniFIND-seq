#!/usr/bin/env python3
"""Summarize manuscript raw read/UMI QC from gene-pipeline outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import os
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Collect featureCounts-assigned alignments and deduplicated UMI totals, "
            "then apply the manuscript's inclusive raw-QC thresholds."
        )
    )
    parser.add_argument("--pipeline-output", required=True)
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--assigned-read-threshold", type=int, default=478877)
    parser.add_argument("--deduplicated-umi-threshold", type=int, default=12149)
    parser.add_argument("--assigned-read-percentile", type=float, default=0.65)
    parser.add_argument("--deduplicated-umi-percentile", type=float, default=0.70)
    return parser.parse_args()


def read_samples(path: Path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "sample" not in reader.fieldnames:
            raise ValueError(f"Manifest must contain a sample column: {path}")
        samples = [row["sample"].strip() for row in reader if row["sample"].strip()]
    if not samples:
        raise ValueError(f"Manifest contains no samples: {path}")
    if len(samples) != len(set(samples)):
        raise ValueError(f"Manifest contains duplicate samples: {path}")
    return samples


def assigned_count(path: Path):
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = {row[0]: row[1] for row in reader if len(row) >= 2}
    if "Assigned" not in rows:
        raise ValueError(f"featureCounts summary has no Assigned row: {path}")
    try:
        value = int(rows["Assigned"])
    except ValueError as exc:
        raise ValueError(f"Invalid Assigned count in {path}: {rows['Assigned']!r}") from exc
    if value < 0:
        raise ValueError(f"Assigned count must be non-negative: {path}")
    return value


def umi_total(path: Path, expected_cell: str):
    total = 0
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"cell", "count"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"UMI table must contain cell and count columns: {path}")
        for line_number, row in enumerate(reader, start=2):
            if row.get("cell", "").strip() != expected_cell:
                raise ValueError(
                    f"Unexpected cell ID at {path}:{line_number}; "
                    f"expected {expected_cell!r}, found {row.get('cell', '')!r}"
                )
            raw = row.get("count", "")
            try:
                value = int(raw)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid UMI count at {path}:{line_number}: {raw!r}"
                ) from exc
            if value < 0:
                raise ValueError(f"Negative UMI count at {path}:{line_number}")
            total += value
    return total


def quantile_type7(values, probability):
    """R-compatible quantile(..., type=7), including the one-value case."""
    if not values:
        raise ValueError("Cannot calculate a percentile from an empty sequence")
    ordered = sorted(values)
    if len(ordered) == 1:
        return float(ordered[0])
    index = (len(ordered) - 1) * probability
    lower = math.floor(index)
    fraction = index - lower
    return ordered[lower] + fraction * (ordered[lower + 1] - ordered[lower])


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        with tmp.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def main():
    args = parse_args()
    if args.assigned_read_threshold < 0 or args.deduplicated_umi_threshold < 0:
        raise ValueError("QC thresholds must be non-negative integers")
    for name, value in (
        ("assigned-read-percentile", args.assigned_read_percentile),
        ("deduplicated-umi-percentile", args.deduplicated_umi_percentile),
    ):
        if not 0 <= value <= 1:
            raise ValueError(f"--{name} must be between 0 and 1")

    pipeline = Path(args.pipeline_output).expanduser().resolve()
    manifest = (
        Path(args.manifest).expanduser().resolve()
        if args.manifest
        else pipeline / "run_metadata" / "samples.normalized.tsv"
    )
    samples = read_samples(manifest)
    rows = []
    for sample in samples:
        sample_dir = pipeline / "mapped" / sample
        summary = sample_dir / f"{sample}_gene_assigned.txt.summary"
        counts = sample_dir / f"{sample}_counts.tsv.gz"
        if not summary.is_file():
            raise FileNotFoundError(f"Missing featureCounts summary: {summary}")
        if not counts.is_file():
            raise FileNotFoundError(f"Missing UMI count table: {counts}")
        assigned = assigned_count(summary)
        umis = umi_total(counts, sample)
        rows.append(
            {
                "sample": sample,
                "featurecounts_assigned": assigned,
                "total_deduplicated_umis": umis,
                "assigned_read_pass": int(assigned <= args.assigned_read_threshold),
                "deduplicated_umi_pass": int(umis <= args.deduplicated_umi_threshold),
                "raw_qc_keep": int(
                    assigned <= args.assigned_read_threshold
                    and umis <= args.deduplicated_umi_threshold
                ),
                "raw_qc_status": (
                    "keep"
                    if assigned <= args.assigned_read_threshold
                    and umis <= args.deduplicated_umi_threshold
                    else "raw_qc_removal"
                ),
                "featurecounts_summary_sha256": sha256(summary),
                "umi_count_table_sha256": sha256(counts),
            }
        )

    empirical_assigned = quantile_type7(
        [row["featurecounts_assigned"] for row in rows],
        args.assigned_read_percentile,
    )
    empirical_umi = quantile_type7(
        [row["total_deduplicated_umis"] for row in rows],
        args.deduplicated_umi_percentile,
    )
    out_dir = Path(args.out_dir).expanduser().resolve()
    atomic_tsv(
        out_dir / "raw_qc_by_sample.tsv",
        rows,
        [
            "sample",
            "featurecounts_assigned",
            "total_deduplicated_umis",
            "assigned_read_pass",
            "deduplicated_umi_pass",
            "raw_qc_keep",
            "raw_qc_status",
            "featurecounts_summary_sha256",
            "umi_count_table_sha256",
        ],
    )
    parameter_rows = [
        {"parameter": "manifest", "value": str(manifest)},
        {"parameter": "sample_count", "value": len(samples)},
        {"parameter": "assigned_read_percentile", "value": args.assigned_read_percentile},
        {"parameter": "empirical_assigned_read_percentile_value", "value": f"{empirical_assigned:.10g}"},
        {"parameter": "applied_assigned_read_threshold", "value": args.assigned_read_threshold},
        {"parameter": "deduplicated_umi_percentile", "value": args.deduplicated_umi_percentile},
        {"parameter": "empirical_deduplicated_umi_percentile_value", "value": f"{empirical_umi:.10g}"},
        {"parameter": "applied_deduplicated_umi_threshold", "value": args.deduplicated_umi_threshold},
        {"parameter": "threshold_comparison", "value": "inclusive_le"},
        {"parameter": "retained_sample_count", "value": sum(row["raw_qc_keep"] for row in rows)},
        {"parameter": "removed_sample_count", "value": sum(1 - row["raw_qc_keep"] for row in rows)},
    ]
    atomic_tsv(out_dir / "raw_qc_parameters.tsv", parameter_rows, ["parameter", "value"])
    print(f"[DONE] summarized {len(samples)} samples in {out_dir}")


if __name__ == "__main__":
    main()

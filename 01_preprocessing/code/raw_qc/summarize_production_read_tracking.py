#!/usr/bin/env python3
"""Reproduce the production read-tracking table and emit a corrected companion."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pipeline-output", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args()


def read_samples(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "sample" not in reader.fieldnames:
            raise ValueError(f"Manifest must contain a sample column: {path}")
        samples = [row["sample"].strip() for row in reader if row["sample"].strip()]
    if not samples or len(samples) != len(set(samples)):
        raise ValueError(f"Manifest has no samples or contains duplicates: {path}")
    return samples


def fastq_records(path):
    lines = 0
    with gzip.open(path, "rt") as handle:
        for _ in handle:
            lines += 1
    if lines % 4:
        raise ValueError(f"FASTQ does not contain complete four-line records: {path}")
    return lines // 4


def umi_total(path):
    total = 0
    if not path.is_file():
        return total
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"gene", "cell", "count"}.issubset(reader.fieldnames):
            raise ValueError(f"Unexpected UMI count-table header: {path}")
        for row in reader:
            total += int(row["count"])
    return total


def bam_metrics(path):
    if not path.is_file():
        return {}, 0, 0
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("pysam is required when an assigned BAM is present") from exc

    status_counts = defaultdict(int)
    unique_umis = set()
    tagged_reads = 0
    accepted = {
        "Assigned",
        "Unassigned_NoFeatures",
        "Unassigned_Ambiguity",
        "Unassigned_MultiMapping",
    }
    with pysam.AlignmentFile(path, "rb") as bam:
        for read in bam:
            if not read.is_unmapped and read.reference_id is not None:
                reference = bam.get_reference_name(read.reference_id)
                species = "mouse" if reference.startswith("mouse_") else "human"
                status = read.get_tag("XS") if read.has_tag("XS") else "Unassigned_NoXS"
                if status not in accepted:
                    status = "Unassigned_Other"
                status_counts[f"reads_{species}_{status}"] += 1
            if all(read.has_tag(tag) for tag in ("CB", "UB", "XT", "XS")) and \
                    read.get_tag("XS") == "Assigned":
                unique_umis.add((read.get_tag("CB"), read.get_tag("XT"), read.get_tag("UB")))
                tagged_reads += 1
    return status_counts, tagged_reads, len(unique_umis)


def append(rows, sample, step, value):
    rows.append({"sample": sample, "step": step, "value": value})


def write_rows(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("sample", "step", "value"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    pipeline = Path(args.pipeline_output).expanduser().resolve()
    manifest = Path(args.manifest).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    legacy_rows = []
    corrected_rows = []

    for sample in read_samples(manifest):
        extracted = pipeline / "umi_extracted" / sample
        trimmed = pipeline / "umi_trimmed" / sample
        raw_paths = {
            "raw_umi": (
                extracted / f"{sample}_R1_umi_out.fastq.gz",
                extracted / f"{sample}_R2_umi_out.fastq.gz",
            ),
            "raw_noumi": (
                extracted / f"{sample}_R1_out.fastq.gz",
                extracted / f"{sample}_R2_out.fastq.gz",
            ),
        }
        trimmed_paths = {
            "trimmed_umi": trimmed / f"{sample}_R1_umi_out_val_1.fq.gz",
            "trimmed_noumi": trimmed / f"{sample}_R1_out_val_1.fq.gz",
        }
        for step, pair in raw_paths.items():
            counts = [fastq_records(path) if path.is_file() else 0 for path in pair]
            append(legacy_rows, sample, step, sum(counts))
            append(corrected_rows, sample, step, counts[0])
        for step, path in trimmed_paths.items():
            value = fastq_records(path) if path.is_file() else 0
            append(legacy_rows, sample, step, value)
            append(corrected_rows, sample, step, value)

        mapped = pipeline / "mapped" / sample
        status, tagged, unique = bam_metrics(mapped / f"{sample}_assigned_sorted.bam")
        shared = []
        for step in sorted(status):
            append(shared, sample, step, status[step])
        human = status.get("reads_human_Assigned", 0)
        mouse = status.get("reads_mouse_Assigned", 0)
        if human + mouse:
            append(shared, sample, "Percent_Human_Assigned", f"{human / (human + mouse) * 100:.2f}%")
        append(shared, sample, "Tagged_Reads", tagged)
        append(shared, sample, "UMIs_tagged", unique)
        append(shared, sample, "UMIs_deduplicated", umi_total(mapped / f"{sample}_counts.tsv.gz"))
        legacy_rows.extend(shared)
        corrected_rows.extend(shared)

    write_rows(out_dir / "readtracking_production_legacy.tsv", legacy_rows)
    write_rows(out_dir / "readtracking_corrected_r1.tsv", corrected_rows)
    print(f"[DONE] production and corrected read tracking: {out_dir}")


if __name__ == "__main__":
    main()

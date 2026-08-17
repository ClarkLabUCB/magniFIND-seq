#!/usr/bin/env python3
"""Validate and normalize a paired-end FASTQ sample manifest."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
from pathlib import Path
from typing import TextIO


SAMPLE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate sample IDs and paired-end FASTQ files, then write a manifest "
            "with absolute paths. Relative FASTQ paths are resolved from the input "
            "manifest directory."
        )
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--check-fastq",
        action="store_true",
        help="Read every FASTQ record and verify pairing, format, and record counts.",
    )
    return parser.parse_args()


def open_fastq(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def normalized_read_id(header: str) -> str:
    token = header[1:].split()[0]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def read_record(handle: TextIO, path: Path, record_number: int):
    lines = [handle.readline() for _ in range(4)]
    if not any(lines):
        return None
    if any(line == "" for line in lines):
        raise ValueError(f"Truncated FASTQ record {record_number} in {path}")

    header, sequence, plus, quality = [line.rstrip("\r\n") for line in lines]
    if not header.startswith("@"):
        raise ValueError(f"Invalid FASTQ header at record {record_number} in {path}")
    if not plus.startswith("+"):
        raise ValueError(f"Invalid FASTQ separator at record {record_number} in {path}")
    if len(sequence) != len(quality):
        raise ValueError(
            f"Sequence/quality length mismatch at record {record_number} in {path}"
        )
    return header


def validate_fastq_pair(r1: Path, r2: Path, sample: str) -> int:
    count = 0
    with open_fastq(r1) as r1_handle, open_fastq(r2) as r2_handle:
        while True:
            record_number = count + 1
            h1 = read_record(r1_handle, r1, record_number)
            h2 = read_record(r2_handle, r2, record_number)
            if h1 is None and h2 is None:
                break
            if h1 is None or h2 is None:
                raise ValueError(f"FASTQ record counts differ for sample {sample}")
            if normalized_read_id(h1) != normalized_read_id(h2):
                raise ValueError(
                    f"FASTQ read IDs differ at record {record_number} for sample {sample}: "
                    f"{h1} != {h2}"
                )
            count += 1
    if count == 0:
        raise ValueError(f"FASTQ pair contains no records for sample {sample}")
    return count


def load_manifest(path: Path, check_fastq: bool):
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Manifest not found: {path}")

    rows = []
    seen = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["sample", "r1", "r2"]:
            raise ValueError(
                "Manifest header must contain exactly three tab-separated columns: "
                "sample, r1, r2"
            )
        for line_number, row in enumerate(reader, start=2):
            sample = (row.get("sample") or "").strip()
            r1_value = (row.get("r1") or "").strip()
            r2_value = (row.get("r2") or "").strip()
            if not sample and not r1_value and not r2_value:
                continue
            if not SAMPLE_RE.fullmatch(sample):
                raise ValueError(
                    f"Unsafe sample ID at line {line_number}: {sample!r}. Use letters, "
                    "numbers, periods, underscores, and hyphens only."
                )
            if sample in seen:
                raise ValueError(f"Duplicate sample ID at line {line_number}: {sample}")
            if not r1_value or not r2_value:
                raise ValueError(f"Missing FASTQ path at line {line_number}: {sample}")

            r1 = Path(r1_value).expanduser()
            r2 = Path(r2_value).expanduser()
            if not r1.is_absolute():
                r1 = path.parent / r1
            if not r2.is_absolute():
                r2 = path.parent / r2
            r1 = r1.resolve()
            r2 = r2.resolve()
            if not r1.is_file() or r1.stat().st_size == 0:
                raise FileNotFoundError(f"Missing or empty R1 FASTQ for {sample}: {r1}")
            if not r2.is_file() or r2.stat().st_size == 0:
                raise FileNotFoundError(f"Missing or empty R2 FASTQ for {sample}: {r2}")
            if r1 == r2:
                raise ValueError(f"R1 and R2 refer to the same file for sample {sample}")

            records = validate_fastq_pair(r1, r2, sample) if check_fastq else None
            rows.append((sample, r1, r2, records))
            seen.add(sample)

    if not rows:
        raise ValueError(f"No samples found in manifest: {path}")
    return rows


def main() -> None:
    args = parse_args()
    rows = load_manifest(args.manifest, args.check_fastq)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(f".{args.output.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample", "r1", "r2"])
            for sample, r1, r2, _ in rows:
                writer.writerow([sample, r1, r2])
        os.replace(temporary, args.output)
    finally:
        if temporary.exists():
            temporary.unlink()

    print(f"[OK] Validated {len(rows)} sample(s): {args.manifest}")
    if args.check_fastq:
        total = sum(records or 0 for *_, records in rows)
        print(f"[OK] Validated {total} paired FASTQ record(s)")
    print(f"[OK] Wrote normalized manifest: {args.output}")


if __name__ == "__main__":
    main()

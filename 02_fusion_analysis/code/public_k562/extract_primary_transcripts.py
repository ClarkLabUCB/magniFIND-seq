#!/usr/bin/env python3
"""Recover primary BAM/SAM transcript records as tagged single-end FASTQ."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import subprocess
from pathlib import Path
from urllib.parse import quote


FILTER_FLAGS = 0x900  # secondary (0x100) plus supplementary (0x800)
COMPLEMENT = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Export one 96-bp transcript read per primary BAM/SAM record"
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--input-format", choices=("bam", "sam"), default="bam")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--read-length", type=int, default=96)
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--metadata", required=True)
    return parser.parse_args()


def encode(value):
    return quote(value, safe="-._~")


def tag_values(optional_fields):
    tags = {}
    for field in optional_fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def sam_lines(path: Path, input_format: str, samtools: str):
    if input_format == "sam":
        with path.open() as handle:
            for line in handle:
                yield line
        return
    process = subprocess.Popen(
        [samtools, "view", "-F", str(FILTER_FLAGS), str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None
    for line in process.stdout:
        yield line
    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code:
        raise RuntimeError(f"samtools view failed ({return_code}): {stderr.strip()}")


def main():
    args = parse_args()
    if args.read_length < 1:
        raise ValueError("--read-length must be positive")
    input_path = Path(args.input).expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    fastq_path = Path(args.fastq).expanduser().resolve()
    metadata_path = Path(args.metadata).expanduser().resolve()
    if len({input_path, fastq_path, metadata_path}) != 3:
        raise ValueError("Input BAM/SAM, FASTQ output, and metadata output must differ")
    fastq_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    fastq_tmp = fastq_path.with_name(f".{fastq_path.name}.tmp.{os.getpid()}")
    metadata_tmp = metadata_path.with_name(f".{metadata_path.name}.tmp.{os.getpid()}")
    record_count = 0
    try:
        with open(fastq_tmp, "wb") as raw_fastq, gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_fastq, mtime=0
        ) as gz_fastq, metadata_tmp.open("w", newline="") as meta_handle:
            writer = csv.DictWriter(
                meta_handle,
                fieldnames=["star_read_name", "original_qname", "CB", "UB", "CR", "UR"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for line_number, raw in enumerate(
                sam_lines(input_path, args.input_format, args.samtools), start=1
            ):
                if not raw.strip() or raw.startswith("@"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 11:
                    raise ValueError(f"Malformed SAM record at input line {line_number}")
                try:
                    flag = int(fields[1])
                except ValueError as exc:
                    raise ValueError(f"Invalid SAM flag at input line {line_number}") from exc
                if flag & FILTER_FLAGS:
                    continue
                qname, sequence, quality = fields[0], fields[9], fields[10]
                if sequence == "*" or quality == "*":
                    raise ValueError(f"Missing sequence or quality for {qname!r}")
                if len(sequence) != len(quality):
                    raise ValueError(f"Sequence/quality length mismatch for {qname!r}")
                if flag & 0x10:
                    sequence = sequence.translate(COMPLEMENT)[::-1]
                    quality = quality[::-1]
                if len(sequence) < args.read_length:
                    raise ValueError(
                        f"Read {qname!r} is shorter than {args.read_length} bases"
                    )
                sequence = sequence[: args.read_length]
                quality = quality[: args.read_length]
                tags = tag_values(fields[11:])
                retained = {name: tags.get(name, "") for name in ("CB", "UB", "CR", "UR")}
                star_name = "Q=" + encode(qname) + "".join(
                    f"|{name}={encode(retained[name])}" for name in ("CB", "UB", "CR", "UR")
                )
                gz_fastq.write(
                    f"@{star_name}\n{sequence}\n+\n{quality}\n".encode("ascii")
                )
                writer.writerow(
                    {
                        "star_read_name": star_name,
                        "original_qname": qname,
                        **retained,
                    }
                )
                record_count += 1
        if record_count == 0:
            raise ValueError("No primary transcript records were exported")
        os.replace(fastq_tmp, fastq_path)
        os.replace(metadata_tmp, metadata_path)
    finally:
        for tmp in (fastq_tmp, metadata_tmp):
            if tmp.exists():
                tmp.unlink()
    print(f"[DONE] exported {record_count} primary records")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Create a deterministic, paired, de-identified FASTQ subset."""

import argparse
import gzip
import hashlib
import heapq
import io
import os
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--out-r1", required=True, type=Path)
    parser.add_argument("--out-r2", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--read-pairs", type=int, default=50000)
    parser.add_argument("--seed", default="magnifind-seq-study-example-v1")
    parser.add_argument("--tso", default="TTGCGCAATG")
    parser.add_argument("--umi-length", type=int, default=8)
    args = parser.parse_args()
    if args.read_pairs < 1:
        parser.error("--read-pairs must be positive")
    if args.umi_length < 1:
        parser.error("--umi-length must be positive")
    if not args.sample or any(char not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-" for char in args.sample):
        parser.error("--sample may contain only letters, numbers, periods, underscores, and hyphens")
    if not args.tso:
        parser.error("--tso must not be empty")
    return args


def open_fastq(path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii", newline="")
    return path.open("rt", encoding="ascii", newline="")


def canonical_read_id(header):
    if not header.startswith("@"):
        raise ValueError("FASTQ header does not start with '@'")
    token = header[1:].split()[0]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def read_record(handle, path):
    lines = [handle.readline() for _ in range(4)]
    if not lines[0]:
        if any(lines[1:]):
            raise ValueError(f"Truncated FASTQ record in {path}")
        return None
    if any(line == "" for line in lines[1:]):
        raise ValueError(f"Truncated FASTQ record in {path}")
    record = tuple(line.rstrip("\r\n") for line in lines)
    if not record[2].startswith("+"):
        raise ValueError(f"Malformed FASTQ separator in {path}")
    if len(record[1]) != len(record[3]):
        raise ValueError(f"Sequence/quality length mismatch in {path}")
    return record


def selection_score(seed, read_id):
    digest = hashlib.blake2b(
        f"{seed}\0{read_id}".encode("utf-8"), digest_size=16
    ).digest()
    return int.from_bytes(digest, "big")


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def select_pairs(args):
    selected = []
    total_pairs = 0
    complete_tso_umi = 0
    with open_fastq(args.r1) as r1_handle, open_fastq(args.r2) as r2_handle:
        while True:
            r1_record = read_record(r1_handle, args.r1)
            r2_record = read_record(r2_handle, args.r2)
            if r1_record is None and r2_record is None:
                break
            if r1_record is None or r2_record is None:
                raise ValueError("R1 and R2 contain different numbers of records")

            r1_id = canonical_read_id(r1_record[0])
            r2_id = canonical_read_id(r2_record[0])
            if r1_id != r2_id:
                raise ValueError(
                    f"Paired read identifiers differ at pair {total_pairs + 1}"
                )

            total_pairs += 1
            tso_position = r1_record[1].find(args.tso)
            if (
                tso_position >= 0
                and len(r1_record[1])
                >= tso_position + len(args.tso) + args.umi_length
            ):
                complete_tso_umi += 1

            score = selection_score(args.seed, r1_id)
            item = (-score, -total_pairs, total_pairs, r1_record, r2_record)
            if len(selected) < args.read_pairs:
                heapq.heappush(selected, item)
            elif item > selected[0]:
                heapq.heapreplace(selected, item)

    if total_pairs < args.read_pairs:
        raise ValueError(
            f"Requested {args.read_pairs} pairs, but inputs contain {total_pairs}"
        )
    return sorted(selected, key=lambda item: item[2]), total_pairs, complete_tso_umi


def deterministic_gzip_text(path):
    raw_handle = path.open("wb")
    gzip_handle = gzip.GzipFile(
        filename="", mode="wb", fileobj=raw_handle, compresslevel=6, mtime=0
    )
    text_handle = io.TextIOWrapper(gzip_handle, encoding="ascii", newline="\n")
    return raw_handle, gzip_handle, text_handle


def write_subset(selected, args, out_r1, out_r2):
    r1_raw, r1_gzip, r1_text = deterministic_gzip_text(out_r1)
    r2_raw, r2_gzip, r2_text = deterministic_gzip_text(out_r2)
    try:
        for output_index, (_, _, _, r1_record, r2_record) in enumerate(selected, 1):
            root = f"{args.sample}_read_{output_index:06d}"
            r1_text.write(f"@{root}/1\n{r1_record[1]}\n+\n{r1_record[3]}\n")
            r2_text.write(f"@{root}/2\n{r2_record[1]}\n+\n{r2_record[3]}\n")
    finally:
        r1_text.close()
        r2_text.close()
        r1_gzip.close()
        r2_gzip.close()
        r1_raw.close()
        r2_raw.close()


def main():
    args = parse_args()
    for source in (args.r1, args.r2):
        if not source.is_file() or source.stat().st_size == 0:
            raise SystemExit(f"Missing or empty input FASTQ: {source}")

    outputs = (args.out_r1, args.out_r2, args.summary)
    for path in outputs:
        if path.exists():
            raise SystemExit(f"Refusing to overwrite existing output: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
    temporary = tuple(
        path.with_name(f".{path.name}.tmp.{os.getpid()}") for path in outputs
    )
    try:
        selected, total_pairs, complete_tso_umi = select_pairs(args)
        write_subset(selected, args, temporary[0], temporary[1])
    except (OSError, ValueError) as error:
        for path in temporary:
            if path.exists():
                path.unlink()
        raise SystemExit(f"[ERROR] {error}") from error

    selected_complete = sum(
        record[3][1].find(args.tso) >= 0
        and len(record[3][1])
        >= record[3][1].find(args.tso) + len(args.tso) + args.umi_length
        for record in selected
    )
    metrics = [
        ("sample", args.sample),
        ("selection_method", "smallest_blake2b_read_id_hashes"),
        ("selection_seed", args.seed),
        ("source_r1_sha256", sha256_file(args.r1)),
        ("source_r2_sha256", sha256_file(args.r2)),
        ("source_total_read_pairs", total_pairs),
        ("source_complete_tso_umi_pairs", complete_tso_umi),
        ("selected_read_pairs", len(selected)),
        ("selected_complete_tso_umi_pairs", selected_complete),
        ("tso_sequence", args.tso),
        ("umi_length", args.umi_length),
        ("output_r1_sha256", sha256_file(temporary[0])),
        ("output_r2_sha256", sha256_file(temporary[1])),
    ]
    try:
        with temporary[2].open("w", encoding="utf-8", newline="") as handle:
            handle.write("metric\tvalue\n")
            for metric, value in metrics:
                handle.write(f"{metric}\t{value}\n")
        # The summary is promoted last and acts as the completion record.
        os.replace(temporary[0], args.out_r1)
        os.replace(temporary[1], args.out_r2)
        os.replace(temporary[2], args.summary)
    finally:
        for path in temporary:
            if path.exists():
                path.unlink()

    print(
        f"[OK] Wrote {len(selected)} de-identified pairs for {args.sample}; "
        f"{selected_complete} contain complete TSO+UMI structure"
    )


if __name__ == "__main__":
    main()

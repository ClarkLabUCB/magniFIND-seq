#!/usr/bin/env python3
"""Write SHA-256 identities for files/directories referenced by a TSV manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path


def file_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--path-column", default="path")
    parser.add_argument("--label-column", default="dataset")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    manifest = args.manifest.expanduser().resolve()
    rows = []
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {args.path_column, args.label_column}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Manifest lacks columns {sorted(required)}: {manifest}")
        for row_number, row in enumerate(reader, start=2):
            raw = row[args.path_column].strip()
            label = row[args.label_column].strip()
            path = Path(raw).expanduser()
            if not path.is_absolute():
                path = manifest.parent / path
            path = path.resolve()
            if not path.exists():
                raise FileNotFoundError(f"Missing input at manifest line {row_number}: {path}")
            files = [path] if path.is_file() else sorted(
                p for p in path.rglob("*") if p.is_file()
            )
            if not files:
                raise ValueError(f"Input directory contains no files: {path}")
            for file_path in files:
                rows.append(
                    {
                        "manifest_row": row_number,
                        "label": label,
                        "input_path": str(path),
                        "relative_file": "." if path.is_file() else str(file_path.relative_to(path)),
                        "bytes": file_path.stat().st_size,
                        "sha256": file_hash(file_path),
                    }
                )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(f".{args.output.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "manifest_row", "label", "input_path", "relative_file",
                    "bytes", "sha256",
                ],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, args.output)
    finally:
        if temporary.exists():
            temporary.unlink()
    print(f"[DONE] hashed {len(rows)} manifest input file(s)")


if __name__ == "__main__":
    main()

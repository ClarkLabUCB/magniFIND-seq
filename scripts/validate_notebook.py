#!/usr/bin/env python3
"""Validate committed Jupyter notebooks without requiring Jupyter packages."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


PRIVATE_MARKERS = ("/Users/", "/Volumes/SW_Data/", "BEGIN OPENSSH PRIVATE KEY")


def validate(path: Path) -> None:
    with path.open(encoding="utf-8") as handle:
        notebook = json.load(handle)

    if notebook.get("nbformat") != 4:
        raise ValueError(f"{path}: expected notebook format 4")
    cells = notebook.get("cells")
    if not isinstance(cells, list) or not cells:
        raise ValueError(f"{path}: notebook must contain cells")

    code_cells = 0
    for index, cell in enumerate(cells, start=1):
        cell_type = cell.get("cell_type")
        if cell_type not in {"code", "markdown", "raw"}:
            raise ValueError(f"{path}: cell {index} has invalid type {cell_type!r}")
        source = cell.get("source", "")
        if isinstance(source, list):
            source = "".join(source)
        if not isinstance(source, str):
            raise ValueError(f"{path}: cell {index} source is not text")
        for marker in PRIVATE_MARKERS:
            if marker in source:
                raise ValueError(f"{path}: cell {index} contains private marker {marker!r}")
        if cell_type == "code":
            code_cells += 1
            compile(source, f"{path}:cell-{index}", "exec")
            if cell.get("execution_count") is not None or cell.get("outputs") not in ([], None):
                raise ValueError(
                    f"{path}: cell {index} contains committed execution state or output"
                )

    if code_cells == 0:
        raise ValueError(f"{path}: notebook contains no executable cells")
    kernelspec = notebook.get("metadata", {}).get("kernelspec", {})
    if kernelspec.get("name") != "python3":
        raise ValueError(f"{path}: kernelspec must be python3")
    print(f"[OK] {path}: {len(cells)} cells, {code_cells} executable")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("notebooks", nargs="+", type=Path)
    args = parser.parse_args()
    for path in args.notebooks:
        validate(path)


if __name__ == "__main__":
    main()

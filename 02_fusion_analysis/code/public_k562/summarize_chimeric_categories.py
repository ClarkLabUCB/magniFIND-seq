#!/usr/bin/env python3
"""Annotate STAR junction ends and quantify manuscript chimeric categories."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
from collections import defaultdict
from pathlib import Path


BIN_SIZE = 1_000_000
CATEGORIES = ("whole_intergene", "BCR_any", "ABL1_any", "BCR_ABL1")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--junctions", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--read-metadata", required=True)
    parser.add_argument("--called-barcodes", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args()


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def attributes(text):
    result = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        result[key] = value.strip().strip('"')
    return result


def load_gene_bins(path):
    bins = defaultdict(list)
    genes = 0
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed GTF record at {path}:{line_number}")
            if fields[2] != "gene":
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f"Invalid GTF coordinates at {path}:{line_number}") from exc
            attrs = attributes(fields[8])
            gene_id = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", gene_id)
            if not gene_id or not gene_name or start < 1 or end < start:
                raise ValueError(f"Invalid gene feature at {path}:{line_number}")
            record = (start, end, gene_id, gene_name, fields[6])
            for bucket in range(start // BIN_SIZE, end // BIN_SIZE + 1):
                bins[(fields[0], bucket)].append(record)
            genes += 1
    if not genes:
        raise ValueError(f"GTF contains no gene features: {path}")
    return bins


def overlap_genes(bins, chrom, position):
    matches = {
        (gene_id, gene_name, strand)
        for start, end, gene_id, gene_name, strand in bins.get(
            (chrom, position // BIN_SIZE), []
        )
        if start <= position <= end
    }
    return sorted(matches)


def load_metadata(path):
    by_star_name = {}
    called_pairs_all = set()
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"star_read_name", "original_qname", "CB", "UB", "CR", "UR"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Read metadata lacks required columns: {path}")
        for row in reader:
            key = row["star_read_name"]
            if not key:
                raise ValueError(f"Empty star_read_name in {path}")
            identity = {name: row[name] for name in required}
            if key in by_star_name and by_star_name[key] != identity:
                raise ValueError(f"Conflicting metadata for STAR read name: {key}")
            by_star_name[key] = identity
            if row["CB"] and row["UB"]:
                called_pairs_all.add((row["CB"], row["UB"]))
    if not by_star_name:
        raise ValueError(f"No read metadata records: {path}")
    return by_star_name, called_pairs_all


def called_barcodes(path):
    values = []
    with open_text(path) as handle:
        for raw in handle:
            value = raw.strip().split("\t")[0]
            if value and value.lower() not in {"barcode", "cell_barcode"}:
                values.append(value)
    if not values or len(values) != len(set(values)):
        raise ValueError("Called-barcode list must be non-empty and unique")
    return set(values)


def atomic_tsv(path, rows, fields):
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
    gtf_bins = load_gene_bins(Path(args.gtf))
    metadata, all_cb_umi = load_metadata(Path(args.read_metadata))
    called = called_barcodes(Path(args.called_barcodes))
    total_called_pairs = {(cb, ub) for cb, ub in all_cb_umi if cb in called}
    annotated = []
    support = {category: {} for category in CATEGORIES}
    malformed = 0
    with open_text(Path(args.junctions)) as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 14:
                malformed += 1
                continue
            try:
                pos1, pos2 = int(fields[1]), int(fields[4])
            except ValueError:
                malformed += 1
                continue
            star_name = fields[9]
            if star_name not in metadata:
                raise ValueError(
                    f"Junction line {line_number} has no matching read metadata: {star_name}"
                )
            meta = metadata[star_name]
            genes1 = overlap_genes(gtf_bins, fields[0], pos1)
            genes2 = overlap_genes(gtf_bins, fields[3], pos2)
            for gene1 in genes1:
                for gene2 in genes2:
                    if gene1[0] == gene2[0]:
                        continue
                    names = {gene1[1], gene2[1]}
                    categories = {"whole_intergene"}
                    if "BCR" in names:
                        categories.add("BCR_any")
                    if "ABL1" in names:
                        categories.add("ABL1_any")
                    if names == {"BCR", "ABL1"}:
                        categories.add("BCR_ABL1")
                    annotated.append(
                        {
                            "source_line": line_number,
                            "original_qname": meta["original_qname"],
                            "CB": meta["CB"],
                            "UB": meta["UB"],
                            "chrom1": fields[0],
                            "pos1": pos1,
                            "strand1": fields[2],
                            "gene1_id": gene1[0],
                            "gene1_name": gene1[1],
                            "gene1_strand": gene1[2],
                            "chrom2": fields[3],
                            "pos2": pos2,
                            "strand2": fields[5],
                            "gene2_id": gene2[0],
                            "gene2_name": gene2[1],
                            "gene2_strand": gene2[2],
                            "categories": ",".join(sorted(categories)),
                        }
                    )
                    for category in categories:
                        support[category][meta["original_qname"]] = meta
    if malformed:
        raise ValueError(f"Encountered {malformed} malformed STAR junction records")

    support_rows = []
    summary_rows = []
    for category in CATEGORIES:
        records = support[category]
        for qname in sorted(records):
            meta = records[qname]
            support_rows.append(
                {
                    "category": category,
                    "original_qname": qname,
                    "CB": meta["CB"],
                    "UB": meta["UB"],
                    "called_cell": int(meta["CB"] in called),
                }
            )
        called_records = [meta for meta in records.values() if meta["CB"] in called]
        pairs = {(meta["CB"], meta["UB"]) for meta in called_records if meta["UB"]}
        cells = {meta["CB"] for meta in called_records}
        summary_rows.append(
            {
                "category": category,
                "physical_reads_all": len(records),
                "physical_reads_called_cells": len(called_records),
                "fusion_cb_umi_pairs": len(pairs),
                "total_called_cb_umi_pairs": len(total_called_pairs),
                "fusion_umi_percent": (
                    f"{100 * len(pairs) / len(total_called_pairs):.12g}"
                    if total_called_pairs
                    else "NA"
                ),
                "positive_called_cells": len(cells),
                "total_called_cells": len(called),
                "positive_cell_percent": f"{100 * len(cells) / len(called):.12g}",
            }
        )

    out = Path(args.out_dir).expanduser().resolve()
    atomic_tsv(
        out / "annotated_intergenic_junctions.tsv",
        annotated,
        [
            "source_line", "original_qname", "CB", "UB", "chrom1", "pos1",
            "strand1", "gene1_id", "gene1_name", "gene1_strand", "chrom2",
            "pos2", "strand2", "gene2_id", "gene2_name", "gene2_strand",
            "categories",
        ],
    )
    atomic_tsv(
        out / "category_read_support.tsv",
        support_rows,
        ["category", "original_qname", "CB", "UB", "called_cell"],
    )
    atomic_tsv(
        out / "category_summary.tsv",
        summary_rows,
        [
            "category", "physical_reads_all", "physical_reads_called_cells",
            "fusion_cb_umi_pairs", "total_called_cb_umi_pairs", "fusion_umi_percent",
            "positive_called_cells", "total_called_cells", "positive_cell_percent",
        ],
    )
    print(f"[DONE] wrote public K562 chimeric summaries to {out}")


if __name__ == "__main__":
    main()

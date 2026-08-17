#!/usr/bin/env python3
"""Screen STAR chimeric junctions for strand-aware BCR::ABL1 evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import subprocess
from pathlib import Path


UB_RE = re.compile(r"(?:^|\|)UB:([^|]+)")
PROJECT_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Scan STAR Chimeric.out.junction files for BCR::ABL1 candidate "
            "junctions in either STAR report order while enforcing gene strands."
        )
    )
    parser.add_argument("--config", required=True, help="Trusted shell-style config file")
    parser.add_argument("--out-dir", default=None, help="Override OUT_DIR from config")
    parser.add_argument("--project-name", default=None, help="Override PROJECT_NAME")
    return parser.parse_args()


def parse_shell_config(path):
    """Source a trusted Bash config and capture its fully expanded environment."""
    config_path = Path(path).expanduser().resolve()
    command = 'set -a; source "$1"; env -0'
    result = subprocess.run(
        ["bash", "-c", command, "bash", str(config_path)],
        cwd=config_path.parent,
        check=True,
        stdout=subprocess.PIPE,
    )
    values = {}
    for entry in result.stdout.split(b"\0"):
        if b"=" not in entry:
            continue
        key, value = entry.split(b"=", 1)
        values[key.decode()] = value.decode()
    return values, config_path.parent


def as_int(values, key, default):
    try:
        return int(values.get(key, default))
    except ValueError as exc:
        raise ValueError(f"{key} must be an integer: {values.get(key)!r}") from exc


def config_path(value, config_dir):
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = config_dir / path
    return path.resolve()


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def in_interval(chrom, pos, target_chrom, start, end):
    return chrom == target_chrom and start <= pos <= end


def parse_ub(read_name):
    match = UB_RE.search(read_name)
    return match.group(1) if match else ""


def classify_pair(chrom1, pos1, strand1, chrom2, pos2, strand2, cfg):
    bcr1 = in_interval(chrom1, pos1, cfg["BCR_CHROM"], cfg["BCR_START"], cfg["BCR_END"])
    bcr2 = in_interval(chrom2, pos2, cfg["BCR_CHROM"], cfg["BCR_START"], cfg["BCR_END"])
    abl1 = in_interval(chrom1, pos1, cfg["ABL1_CHROM"], cfg["ABL1_START"], cfg["ABL1_END"])
    abl2 = in_interval(chrom2, pos2, cfg["ABL1_CHROM"], cfg["ABL1_START"], cfg["ABL1_END"])

    if bcr1 and abl2:
        report_order = "BCR_then_ABL1"
        bcr = (chrom1, pos1, strand1)
        abl = (chrom2, pos2, strand2)
    elif abl1 and bcr2:
        report_order = "ABL1_then_BCR"
        bcr = (chrom2, pos2, strand2)
        abl = (chrom1, pos1, strand1)
    else:
        return {
            "interval_match": False,
            "strand_match": False,
            "report_order": "outside_target_intervals",
            "bcr": ("", "", ""),
            "abl": ("", "", ""),
        }

    strand_match = bcr[2] == cfg["BCR_STRAND"] and abl[2] == cfg["ABL1_STRAND"]
    return {
        "interval_match": True,
        "strand_match": strand_match,
        "report_order": report_order,
        "bcr": bcr,
        "abl": abl,
    }


def read_chimeric_file(path, sample, cfg):
    rows_chr9_chr22 = []
    rows_candidates = []
    malformed_records = 0
    total_records = 0

    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            raw = raw.rstrip("\n")
            # STAR --chimOutJunctionFormat 1 prepends comment/header records.
            # They describe the format and are not chimeric alignment records.
            if not raw or raw.lstrip().startswith("#"):
                continue
            total_records += 1
            parts = raw.split("\t")
            if len(parts) < 14:
                malformed_records += 1
                continue
            try:
                pos1 = int(parts[1])
                pos2 = int(parts[4])
            except ValueError:
                malformed_records += 1
                continue

            chrom1, strand1 = parts[0], parts[2]
            chrom2, strand2 = parts[3], parts[5]
            if {chrom1, chrom2} != {cfg["BCR_CHROM"], cfg["ABL1_CHROM"]}:
                continue

            classification = classify_pair(
                chrom1, pos1, strand1, chrom2, pos2, strand2, cfg
            )
            candidate = classification["interval_match"] and classification["strand_match"]
            bcr_chrom, bcr_pos, bcr_strand = classification["bcr"]
            abl_chrom, abl_pos, abl_strand = classification["abl"]
            read_name = parts[9]
            row = {
                "sample": sample,
                "source_line": line_number,
                "chrom1": chrom1,
                "pos1": pos1,
                "strand1": strand1,
                "chrom2": chrom2,
                "pos2": pos2,
                "strand2": strand2,
                "junction_type": parts[6],
                "repeat_left": parts[7],
                "repeat_right": parts[8],
                "read_name": read_name,
                "UB": parse_ub(read_name),
                "start1": parts[10],
                "cigar1": parts[11],
                "start2": parts[12],
                "cigar2": parts[13],
                "report_order": classification["report_order"],
                "bcr_abl_interval_match": int(classification["interval_match"]),
                "bcr_abl_strand_match": int(classification["strand_match"]),
                "bcr_abl_expected_orientation": int(candidate),
                "abl_bcr_reverse_report": int(
                    classification["report_order"] == "ABL1_then_BCR"
                ),
                "canonical_bcr_chrom": bcr_chrom,
                "canonical_bcr_pos": bcr_pos,
                "canonical_bcr_strand": bcr_strand,
                "canonical_abl1_chrom": abl_chrom,
                "canonical_abl1_pos": abl_pos,
                "canonical_abl1_strand": abl_strand,
                "raw_record": raw,
            }
            rows_chr9_chr22.append(row)
            if candidate:
                rows_candidates.append(row)

    return rows_chr9_chr22, rows_candidates, total_records, malformed_records


def write_table(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with open(tmp_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    os.replace(tmp_path, path)


def read_manifest_samples(path):
    if not path.is_file():
        raise FileNotFoundError(f"Sample manifest not found: {path}")
    samples = []
    seen = set()
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "sample" not in reader.fieldnames:
            raise ValueError(f"Manifest has no sample column: {path}")
        for row in reader:
            sample = row.get("sample", "").strip()
            if not sample:
                continue
            if not PROJECT_RE.fullmatch(sample):
                raise ValueError(f"Unsafe sample ID in manifest: {sample!r}")
            if sample in seen:
                raise ValueError(f"Duplicate sample in manifest: {sample}")
            samples.append(sample)
            seen.add(sample)
    return samples


def main():
    args = parse_args()
    values, config_dir = parse_shell_config(args.config)

    project_name = args.project_name or values.get("PROJECT_NAME", "magnifind_seq")
    if not PROJECT_RE.fullmatch(project_name):
        raise ValueError(f"Unsafe PROJECT_NAME: {project_name!r}")
    out_value = args.out_dir or values.get("OUT_DIR", "output/fusion_detection")
    out_dir = config_path(out_value, config_dir)
    mapped_dir = out_dir / "fusion_mapped"
    scan_dir = out_dir / "fusion_scan"

    cfg = {
        "BCR_CHROM": values.get("BCR_CHROM", "chr22"),
        "BCR_START": as_int(values, "BCR_START", 23179704),
        "BCR_END": as_int(values, "BCR_END", 23318037),
        "BCR_STRAND": values.get("BCR_STRAND", "+"),
        "ABL1_CHROM": values.get("ABL1_CHROM", "chr9"),
        "ABL1_START": as_int(values, "ABL1_START", 130713946),
        "ABL1_END": as_int(values, "ABL1_END", 130887675),
        "ABL1_STRAND": values.get("ABL1_STRAND", "+"),
    }
    for key in ("BCR_STRAND", "ABL1_STRAND"):
        if cfg[key] not in {"+", "-"}:
            raise ValueError(f"{key} must be '+' or '-': {cfg[key]!r}")
    for gene in ("BCR", "ABL1"):
        if cfg[f"{gene}_START"] > cfg[f"{gene}_END"]:
            raise ValueError(f"{gene}_START must be <= {gene}_END")

    min_unique_umi = as_int(values, "MIN_UNIQUE_UMI", 2)
    if min_unique_umi < 1:
        raise ValueError("MIN_UNIQUE_UMI must be at least 1")
    require_mapping_done = as_int(values, "REQUIRE_MAPPING_DONE", 1)
    fail_on_missing = as_int(values, "FAIL_ON_MISSING_JUNCTIONS", 1)
    fail_on_malformed = as_int(values, "FAIL_ON_MALFORMED_RECORDS", 1)
    for key, value in (
        ("REQUIRE_MAPPING_DONE", require_mapping_done),
        ("FAIL_ON_MISSING_JUNCTIONS", fail_on_missing),
        ("FAIL_ON_MALFORMED_RECORDS", fail_on_malformed),
    ):
        if value not in {0, 1}:
            raise ValueError(f"{key} must be 0 or 1")

    manifest_value = values.get("SAMPLES_TSV", "")
    if manifest_value:
        manifest_path = config_path(manifest_value, config_dir)
        normalized = out_dir / "run_metadata" / "fusion_samples.normalized.tsv"
        sample_names = read_manifest_samples(normalized if normalized.is_file() else manifest_path)
    else:
        sample_names = []
    if not sample_names:
        sample_names = sorted(p.name for p in mapped_dir.glob("*") if p.is_dir())
    for sample in sample_names:
        if not PROJECT_RE.fullmatch(sample):
            raise ValueError(f"Unsafe sample ID: {sample!r}")
    if not sample_names:
        raise ValueError(f"No samples found in manifest or {mapped_dir}")

    all_chr9_chr22 = []
    all_candidates = []
    summary = []
    evaluation_errors = []

    for sample in sample_names:
        sample_dir = mapped_dir / sample
        junction_path = sample_dir / f"{sample}_Chimeric.out.junction"
        done_path = sample_dir / f"{sample}_fusion_mapping.done"

        if not junction_path.exists():
            mapping_status = "missing_junction_file"
            call = "not_evaluated_missing_junction"
            chr_rows, candidate_rows = [], []
            total_records = malformed_records = 0
            evaluation_errors.append(f"{sample}: missing {junction_path}")
        elif require_mapping_done and (
            not done_path.is_file() or done_path.stat().st_size == 0
        ):
            mapping_status = "incomplete_mapping"
            call = "not_evaluated_incomplete_mapping"
            chr_rows, candidate_rows = [], []
            total_records = malformed_records = 0
            evaluation_errors.append(f"{sample}: missing completion marker {done_path}")
        else:
            chr_rows, candidate_rows, total_records, malformed_records = read_chimeric_file(
                junction_path, sample, cfg
            )
            if malformed_records:
                mapping_status = "malformed_junction_records"
                call = "not_evaluated_malformed_junctions"
                evaluation_errors.append(
                    f"{sample}: {malformed_records} malformed junction record(s)"
                )
            else:
                mapping_status = "complete"
                unique_ub = {row["UB"] for row in candidate_rows if row["UB"]}
                if not candidate_rows:
                    call = "not_detected"
                elif len(unique_ub) >= min_unique_umi:
                    call = "candidate_ge_min_umi"
                else:
                    call = "candidate_low_umi"

        all_chr9_chr22.extend(chr_rows)
        all_candidates.extend(candidate_rows)
        read_names = {row["read_name"] for row in candidate_rows}
        unique_ub = {row["UB"] for row in candidate_rows if row["UB"]}
        junctions = {
            (
                row["canonical_bcr_chrom"],
                row["canonical_bcr_pos"],
                row["canonical_bcr_strand"],
                row["canonical_abl1_chrom"],
                row["canonical_abl1_pos"],
                row["canonical_abl1_strand"],
            )
            for row in candidate_rows
        }
        summary.append(
            {
                "sample": sample,
                "mapping_status": mapping_status,
                "chimeric_records_total": total_records,
                "malformed_records": malformed_records,
                "expected_orientation_read_records": len(candidate_rows),
                "unique_read_names": len(read_names),
                "unique_UB": len(unique_ub),
                "unique_junctions": len(junctions),
                "expected_BCR_ABL1_candidate": (
                    int(bool(candidate_rows)) if mapping_status == "complete" else ""
                ),
                "expected_BCR_ABL1_candidate_call": call,
            }
        )

    common_fields = [
        "sample",
        "source_line",
        "chrom1",
        "pos1",
        "strand1",
        "chrom2",
        "pos2",
        "strand2",
        "junction_type",
        "repeat_left",
        "repeat_right",
        "read_name",
        "UB",
        "start1",
        "cigar1",
        "start2",
        "cigar2",
        "report_order",
        "bcr_abl_interval_match",
        "bcr_abl_strand_match",
        "bcr_abl_expected_orientation",
        "abl_bcr_reverse_report",
        "canonical_bcr_chrom",
        "canonical_bcr_pos",
        "canonical_bcr_strand",
        "canonical_abl1_chrom",
        "canonical_abl1_pos",
        "canonical_abl1_strand",
        "raw_record",
    ]
    write_table(
        scan_dir / f"{project_name}_chr9_chr22_chimeric_junctions.tsv",
        all_chr9_chr22,
        common_fields,
    )
    write_table(
        scan_dir / f"{project_name}_BCR_ABL1_expected_orientation_junctions.tsv",
        all_candidates,
        common_fields,
    )
    write_table(
        scan_dir / f"{project_name}_BCR_ABL1_expected_orientation_summary_all_samples.tsv",
        summary,
        [
            "sample",
            "mapping_status",
            "chimeric_records_total",
            "malformed_records",
            "expected_orientation_read_records",
            "unique_read_names",
            "unique_UB",
            "unique_junctions",
            "expected_BCR_ABL1_candidate",
            "expected_BCR_ABL1_candidate_call",
        ],
    )
    write_table(
        scan_dir / f"{project_name}_fusion_scan_parameters.tsv",
        [
            {"parameter": "config", "value": str(Path(args.config).resolve())},
            {"parameter": "sample_manifest", "value": str(manifest_path) if manifest_value else ""},
            {"parameter": "BCR_CHROM", "value": cfg["BCR_CHROM"]},
            {"parameter": "BCR_START", "value": cfg["BCR_START"]},
            {"parameter": "BCR_END", "value": cfg["BCR_END"]},
            {"parameter": "BCR_STRAND", "value": cfg["BCR_STRAND"]},
            {"parameter": "ABL1_CHROM", "value": cfg["ABL1_CHROM"]},
            {"parameter": "ABL1_START", "value": cfg["ABL1_START"]},
            {"parameter": "ABL1_END", "value": cfg["ABL1_END"]},
            {"parameter": "ABL1_STRAND", "value": cfg["ABL1_STRAND"]},
            {"parameter": "MIN_UNIQUE_UMI", "value": min_unique_umi},
            {"parameter": "REQUIRE_MAPPING_DONE", "value": require_mapping_done},
            {"parameter": "FAIL_ON_MISSING_JUNCTIONS", "value": fail_on_missing},
            {"parameter": "FAIL_ON_MALFORMED_RECORDS", "value": fail_on_malformed},
        ],
        ["parameter", "value"],
    )

    print(f"[DONE] scanned {len(sample_names)} samples")
    print(f"[DONE] wrote {scan_dir}")
    should_fail = (
        fail_on_missing
        and any("missing" in error for error in evaluation_errors)
    ) or (fail_on_malformed and any("malformed" in error for error in evaluation_errors))
    if should_fail:
        for error in evaluation_errors:
            print(f"[ERROR] {error}", file=os.sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()

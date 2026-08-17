#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
example_dir="$repo_root/examples/k562r_single_bead"
tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

echo "[TEST] Validate de-identified K562-r single-bead FASTQs"
python3 "$repo_root/scripts/validate_sample_manifest.py" \
  --manifest "$example_dir/data/samples.tsv" \
  --output "$tmp_dir/normalized.tsv" \
  --check-fastq

python3 - "$repo_root" "$example_dir" <<'PY'
import csv
import gzip
import hashlib
import re
import sys
import types
from pathlib import Path

repo_root = Path(sys.argv[1])
example_dir = Path(sys.argv[2])


def xopen(path, mode):
    path = Path(path)
    if path.name.endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode)


stub = types.ModuleType("xopen")
stub.xopen = xopen
sys.modules["xopen"] = stub
sys.path.insert(0, str(repo_root / "01_preprocessing" / "code" / "gene_expression"))
from umi_methods import SeqExperiment


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


with (example_dir / "data" / "samples.tsv").open() as handle:
    samples = list(csv.DictReader(handle, delimiter="\t"))

for row in samples:
    sample = row["sample"]
    r1 = example_dir / "data" / row["r1"]
    r2 = example_dir / "data" / row["r2"]
    provenance_path = example_dir / "data" / f"{sample}.provenance.tsv"
    with provenance_path.open() as handle:
        provenance = {
            item["metric"]: item["value"]
            for item in csv.DictReader(handle, delimiter="\t")
        }

    assert sha256_file(r1) == provenance["output_r1_sha256"]
    assert sha256_file(r2) == provenance["output_r2_sha256"]
    assert int(provenance["selected_read_pairs"]) == 50000

    header_pattern = re.compile(rf"^@{re.escape(sample)}_read_[0-9]{{6}}/1$")
    with gzip.open(r1, "rt") as handle:
        observed_records = 0
        while True:
            header = handle.readline().rstrip("\r\n")
            if not header:
                break
            sequence = handle.readline().rstrip("\r\n")
            plus = handle.readline().rstrip("\r\n")
            quality = handle.readline().rstrip("\r\n")
            observed_records += 1
            assert header_pattern.fullmatch(header)
            assert plus == "+"
            assert len(sequence) == len(quality) == 100
    assert observed_records == 50000

    experiment = SeqExperiment(r1, r2, tso="TTGCGCAATG", umi_length=8)
    observed_umi_reads = sum(1 for _ in experiment.iter_umi_reads())
    expected_umi_reads = int(provenance["selected_complete_tso_umi_pairs"])
    assert experiment.total_reads == 50000
    assert experiment.umi_reads == observed_umi_reads == expected_umi_reads
    print(
        f"[OK] {sample}: 50000 pairs, "
        f"{observed_umi_reads} complete TSO+UMI pairs"
    )
PY

gzip -t "$example_dir"/data/fastq/*.fastq.gz
echo "[PASS] Study-derived FASTQ fixture integrity and UMI structure matched provenance"

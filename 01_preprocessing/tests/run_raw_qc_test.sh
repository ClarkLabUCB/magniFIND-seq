#!/usr/bin/env bash
set -euo pipefail

analysis_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/pipeline/run_metadata"
printf 'sample\tr1\tr2\nA\ta\tb\nB\ta\tb\nC\ta\tb\nD\ta\tb\n' > "$tmp/pipeline/run_metadata/samples.normalized.tsv"
for sample in A B C D; do mkdir -p "$tmp/pipeline/mapped/$sample"; done
mkdir -p "$tmp/pipeline/umi_extracted/A" "$tmp/pipeline/umi_trimmed/A"
for path in \
  "$tmp/pipeline/umi_extracted/A/A_R1_umi_out.fastq.gz" \
  "$tmp/pipeline/umi_extracted/A/A_R2_umi_out.fastq.gz" \
  "$tmp/pipeline/umi_trimmed/A/A_R1_umi_out_val_1.fq.gz"; do
  printf '@read1\nACGT\n+\nIIII\n' | gzip -n > "$path"
done
for pair in A:100 B:478877 C:478878 D:900000; do
  sample="${pair%%:*}"; value="${pair#*:}"
  printf 'Status\t%s\nAssigned\t%s\nUnassigned_NoFeatures\t0\n' "$sample" "$value" > \
    "$tmp/pipeline/mapped/$sample/${sample}_gene_assigned.txt.summary"
done
for pair in A:100 B:12149 C:12150 D:20000; do
  sample="${pair%%:*}"; value="${pair#*:}"
  printf 'gene\tcell\tcount\ngene1\t%s\t%s\n' "$sample" "$value" | gzip -n > \
    "$tmp/pipeline/mapped/$sample/${sample}_counts.tsv.gz"
done
cat > "$tmp/config.sh" <<EOF
PIPELINE_OUTPUT="$tmp/pipeline"
OUT_DIR="$tmp/output"
EOF
bash "$analysis_dir/code/run_raw_qc.sh" "$tmp/config.sh"
python3 - "$tmp/output" <<'PY'
import csv
import sys
from pathlib import Path

out = Path(sys.argv[1])
with (out / "raw_qc_by_sample.tsv").open(newline="") as handle:
    rows = {r["sample"]: r for r in csv.DictReader(handle, delimiter="\t")}
assert rows["A"]["raw_qc_status"] == "keep"
assert rows["B"]["raw_qc_status"] == "keep"  # inclusive boundary
assert rows["C"]["raw_qc_status"] == "raw_qc_removal"
assert rows["D"]["raw_qc_status"] == "raw_qc_removal"
assert all(len(row["featurecounts_summary_sha256"]) == 64 for row in rows.values())
assert all(len(row["umi_count_table_sha256"]) == 64 for row in rows.values())
with (out / "raw_qc_parameters.tsv").open(newline="") as handle:
    params = {r["parameter"]: r["value"] for r in csv.DictReader(handle, delimiter="\t")}
assert params["sample_count"] == "4"
assert params["retained_sample_count"] == "2"
assert params["removed_sample_count"] == "2"
assert params["empirical_assigned_read_percentile_value"] == "478877.95"
assert params["empirical_deduplicated_umi_percentile_value"] == "12935"
assert (out / "run_metadata/config.sh").is_file()
assert (out / "run_metadata/source_provenance.tsv").is_file()
for name in ("readtracking_production_legacy.tsv", "readtracking_corrected_r1.tsv"):
    assert (out / name).is_file()
def tracking(path):
    with path.open(newline="") as handle:
        return {
            (row["sample"], row["step"]): row["value"]
            for row in csv.DictReader(handle, delimiter="\t")
        }
legacy = tracking(out / "readtracking_production_legacy.tsv")
corrected = tracking(out / "readtracking_corrected_r1.tsv")
assert legacy[("A", "raw_umi")] == "2"       # archived code counted R1 + R2
assert corrected[("A", "raw_umi")] == "1"    # corrected metric counts R1
assert legacy[("A", "trimmed_umi")] == corrected[("A", "trimmed_umi")] == "1"
PY

echo "[TEST] Reject a UMI table assigned to a different sample"
printf 'gene\tcell\tcount\ngene1\tWRONG_SAMPLE\t100\n' | gzip -n > \
  "$tmp/pipeline/mapped/A/A_counts.tsv.gz"
if python3 "$analysis_dir/code/raw_qc/summarize_raw_qc.py" \
    --pipeline-output "$tmp/pipeline" --out-dir "$tmp/invalid_output" \
    2> "$tmp/wrong_cell.err"; then
  echo "[ERROR] A mismatched UMI-table cell ID was accepted" >&2
  exit 1
fi
grep -q "Unexpected cell ID" "$tmp/wrong_cell.err"
echo "[PASS] Raw-QC fixture passed"

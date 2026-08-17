#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
analysis_dir="$(cd "$script_dir/../.." && pwd)"

tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

mkdir -p "$tmp_dir/data"
cp -R "$script_dir/data/k562_r_fusion_fixture" "$tmp_dir/data/"

config_file="$tmp_dir/config.sh"
cat > "$config_file" <<EOF
FIXTURE_ROOT="$tmp_dir/data/k562_r_fusion_fixture"
PROJECT_NAME="k562_r_fusion_fixture_test"
SAMPLES_TSV="\${FIXTURE_ROOT}/samples.tsv"
OUT_DIR="\${FIXTURE_ROOT}"
BCR_CHROM="chr22"
BCR_START=23179704
BCR_END=23318037
BCR_STRAND="+"
ABL1_CHROM="chr9"
ABL1_START=130713946
ABL1_END=130887675
ABL1_STRAND="+"
MIN_UNIQUE_UMI=2
REQUIRE_MAPPING_DONE=0
FAIL_ON_MISSING_JUNCTIONS=1
FAIL_ON_MALFORMED_RECORDS=1
EOF

python3 "$analysis_dir/code/targeted/scan_bcr_abl1_chimeric_junctions.py" --config "$config_file"

scan_dir="$tmp_dir/data/k562_r_fusion_fixture/fusion_scan"
actual="$scan_dir/k562_r_fusion_fixture_test_BCR_ABL1_expected_orientation_summary_all_samples.tsv"
expected="$script_dir/expected/k562_r_fusion_fixture/BCR_ABL1_expected_orientation_summary_all_samples.tsv"
diff -u "$expected" "$actual"

python3 - "$scan_dir" <<'PY'
import csv
import sys
from pathlib import Path

scan_dir = Path(sys.argv[1])
candidate_path = scan_dir / "k562_r_fusion_fixture_test_BCR_ABL1_expected_orientation_junctions.tsv"
chromosome_path = scan_dir / "k562_r_fusion_fixture_test_chr9_chr22_chimeric_junctions.tsv"
with candidate_path.open() as handle:
    candidates = list(csv.DictReader(handle, delimiter="\t"))
with chromosome_path.open() as handle:
    chromosome_rows = list(csv.DictReader(handle, delimiter="\t"))
assert len(candidates) == 131
assert len(chromosome_rows) == 131
assert {row["sample"] for row in candidates} == {"k562_r_candidate_1", "k562_r_candidate_2"}
assert {row["report_order"] for row in candidates} == {"BCR_then_ABL1"}
assert {row["bcr_abl_strand_match"] for row in candidates} == {"1"}
assert all(row["canonical_bcr_chrom"] == "chr22" for row in candidates)
assert all(row["canonical_abl1_chrom"] == "chr9" for row in candidates)
PY

edge_root="$tmp_dir/edge_cases"
mkdir -p "$edge_root/fusion_mapped"/{forward,reverse,wrong_strand,empty}
cat > "$edge_root/samples.tsv" <<'EOF'
sample	r1	r2
forward	unused_R1.fastq.gz	unused_R2.fastq.gz
reverse	unused_R1.fastq.gz	unused_R2.fastq.gz
wrong_strand	unused_R1.fastq.gz	unused_R2.fastq.gz
empty	unused_R1.fastq.gz	unused_R2.fastq.gz
EOF
printf '# STAR chimeric junction format 1 header\n# header records must not be counted as alignments\nchr22\t200\t+\tchr9\t300\t+\t-1\t0\t0\tread_forward|CB:forward|UB:AAAAAAAA\t180\t20M\t300\t20M\n' > "$edge_root/fusion_mapped/forward/forward_Chimeric.out.junction"
printf 'chr9\t300\t+\tchr22\t200\t+\t-1\t0\t0\tread_reverse|CB:reverse|UB:CCCCCCCC\t300\t20M\t180\t20M\n' > "$edge_root/fusion_mapped/reverse/reverse_Chimeric.out.junction"
printf 'chr22\t200\t-\tchr9\t300\t+\t-1\t0\t0\tread_wrong|CB:wrong_strand|UB:GGGGGGGG\t180\t20M\t300\t20M\n' > "$edge_root/fusion_mapped/wrong_strand/wrong_strand_Chimeric.out.junction"
: > "$edge_root/fusion_mapped/empty/empty_Chimeric.out.junction"

edge_config="$tmp_dir/edge.config.sh"
cat > "$edge_config" <<EOF
ROOT="$edge_root"
PROJECT_NAME="edge_cases"
SAMPLES_TSV="\${ROOT}/samples.tsv"
OUT_DIR="\${ROOT}"
BCR_CHROM="chr22"
BCR_START=100
BCR_END=250
BCR_STRAND="+"
ABL1_CHROM="chr9"
ABL1_START=250
ABL1_END=400
ABL1_STRAND="+"
MIN_UNIQUE_UMI=2
REQUIRE_MAPPING_DONE=0
FAIL_ON_MISSING_JUNCTIONS=1
FAIL_ON_MALFORMED_RECORDS=1
EOF
python3 "$analysis_dir/code/targeted/scan_bcr_abl1_chimeric_junctions.py" --config "$edge_config"
diff -u \
  "$script_dir/expected/edge_case_summary.tsv" \
  "$edge_root/fusion_scan/edge_cases_BCR_ABL1_expected_orientation_summary_all_samples.tsv"

python3 - "$edge_root/fusion_scan/edge_cases_BCR_ABL1_expected_orientation_junctions.tsv" <<'PY'
import csv
import sys

with open(sys.argv[1]) as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
assert [row["sample"] for row in rows] == ["forward", "reverse"]
assert [row["report_order"] for row in rows] == ["BCR_then_ABL1", "ABL1_then_BCR"]
assert {row["canonical_bcr_chrom"] for row in rows} == {"chr22"}
assert {row["canonical_abl1_chrom"] for row in rows} == {"chr9"}
PY

missing_root="$tmp_dir/missing_case"
mkdir -p "$missing_root/fusion_mapped/missing_sample"
printf 'sample\tr1\tr2\nmissing_sample\tunused_R1.fastq.gz\tunused_R2.fastq.gz\n' > "$missing_root/samples.tsv"
sed "s|ROOT=\"$edge_root\"|ROOT=\"$missing_root\"|; s/PROJECT_NAME=\"edge_cases\"/PROJECT_NAME=\"missing_case\"/" "$edge_config" > "$tmp_dir/missing.config.sh"
if python3 "$analysis_dir/code/targeted/scan_bcr_abl1_chimeric_junctions.py" \
    --config "$tmp_dir/missing.config.sh" \
    > "$tmp_dir/missing.stdout" 2> "$tmp_dir/missing.stderr"; then
  echo "[ERROR] Missing junction file was accepted" >&2
  exit 1
fi
grep -q "missing.*Chimeric.out.junction" "$tmp_dir/missing.stderr"
grep -q $'missing_sample\tmissing_junction_file\t0\t0\t0\t0\t0\t0\t\tnot_evaluated_missing_junction' \
  "$missing_root/fusion_scan/missing_case_BCR_ABL1_expected_orientation_summary_all_samples.tsv"

echo "[PASS] Fusion scanner fixtures, report-order handling, strand checks, and missing-input checks passed"

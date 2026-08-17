#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

echo "[TEST] Validate the embedded paired FASTQ manifest"
python3 "$repo_root/scripts/validate_sample_manifest.py" \
  --manifest "$repo_root/tests/fastq_integration/data/samples.tsv" \
  --output "$tmp_dir/normalized.tsv" \
  --check-fastq

fixture_r1="$repo_root/tests/fastq_integration/data/fastq/k562_r_fastq_fixture_R1.fastq.gz"
fixture_r2="$repo_root/tests/fastq_integration/data/fastq/k562_r_fastq_fixture_R2.fastq.gz"
cat > "$tmp_dir/duplicate.tsv" <<EOF
sample	r1	r2
duplicate	$fixture_r1	$fixture_r2
duplicate	$fixture_r1	$fixture_r2
EOF
if python3 "$repo_root/scripts/validate_sample_manifest.py" \
    --manifest "$tmp_dir/duplicate.tsv" --output "$tmp_dir/invalid.tsv" \
    2> "$tmp_dir/duplicate.err"; then
  echo "[ERROR] Duplicate sample IDs were accepted" >&2
  exit 1
fi
grep -q "Duplicate sample ID" "$tmp_dir/duplicate.err"

echo "[TEST] Reject malformed FASTQ sequence/quality lengths"
cat > "$tmp_dir/bad_R1.fastq" <<'EOF'
@bad/1
ACGT
+
III
EOF
cat > "$tmp_dir/bad_R2.fastq" <<'EOF'
@bad/2
TGCA
+
IIII
EOF
cat > "$tmp_dir/bad_fastq.tsv" <<EOF
sample	r1	r2
bad	$tmp_dir/bad_R1.fastq	$tmp_dir/bad_R2.fastq
EOF
if python3 "$repo_root/scripts/validate_sample_manifest.py" \
    --manifest "$tmp_dir/bad_fastq.tsv" --output "$tmp_dir/bad.normalized.tsv" \
    --check-fastq 2> "$tmp_dir/bad_fastq.err"; then
  echo "[ERROR] Sequence/quality length mismatch was accepted" >&2
  exit 1
fi
grep -q "Sequence/quality length mismatch" "$tmp_dir/bad_fastq.err"

echo "[TEST] Create a byte-deterministic, de-identified paired subset"
python3 "$repo_root/scripts/create_deidentified_fastq_subset.py" \
  --r1 "$fixture_r1" --r2 "$fixture_r2" --read-pairs 2 --sample subset \
  --out-r1 "$tmp_dir/subset_a_R1.fastq.gz" --out-r2 "$tmp_dir/subset_a_R2.fastq.gz" \
  --summary "$tmp_dir/subset_a.tsv"
python3 "$repo_root/scripts/create_deidentified_fastq_subset.py" \
  --r1 "$fixture_r1" --r2 "$fixture_r2" --read-pairs 2 --sample subset \
  --out-r1 "$tmp_dir/subset_b_R1.fastq.gz" --out-r2 "$tmp_dir/subset_b_R2.fastq.gz" \
  --summary "$tmp_dir/subset_b.tsv"
cmp "$tmp_dir/subset_a_R1.fastq.gz" "$tmp_dir/subset_b_R1.fastq.gz"
cmp "$tmp_dir/subset_a_R2.fastq.gz" "$tmp_dir/subset_b_R2.fastq.gz"
grep -q '^@subset_read_000001/1$' <(gzip -cd "$tmp_dir/subset_a_R1.fastq.gz")
if python3 "$repo_root/scripts/create_deidentified_fastq_subset.py" \
    --r1 "$fixture_r1" --r2 "$tmp_dir/bad_R2.fastq" --read-pairs 2 --sample failed \
    --out-r1 "$tmp_dir/failed_R1.fastq.gz" --out-r2 "$tmp_dir/failed_R2.fastq.gz" \
    --summary "$tmp_dir/failed.tsv" 2> "$tmp_dir/failed.err"; then
  echo "[ERROR] Invalid input was accepted by subset creation" >&2
  exit 1
fi
[[ ! -e "$tmp_dir/failed_R1.fastq.gz" && ! -e "$tmp_dir/failed_R2.fastq.gz" && ! -e "$tmp_dir/failed.tsv" ]]

echo "[TEST] Accept conventional /1 and /2 paired read suffixes during UMI extraction"
cat > "$tmp_dir/suffix_R1.fastq" <<'EOF'
@pair001/1
TTGCGCAATGACGTACGTGATTACA
+
IIIIIIIIIIIIIIIIIIIIIIIII
EOF
cat > "$tmp_dir/suffix_R2.fastq" <<'EOF'
@pair001/2
GATTACAGATTACA
+
IIIIIIIIIIIIII
EOF
python3 - "$repo_root" "$tmp_dir" <<'PY'
import sys
import types
import hashlib
from pathlib import Path

repo_root = Path(sys.argv[1])
tmp_dir = Path(sys.argv[2])
stub = types.ModuleType("xopen")
stub.xopen = open
sys.modules["xopen"] = stub
sys.path.insert(0, str(repo_root / "01_preprocessing" / "code" / "gene_expression"))
from umi_methods import SeqExperiment

records = list(
    SeqExperiment(
        tmp_dir / "suffix_R1.fastq",
        tmp_dir / "suffix_R2.fastq",
    ).iter_umi_reads()
)
assert len(records) == 1
assert records[0][0] == "pair001"
assert records[0][1] == "ACGTACGT"

outputs = []
for name in ("extract_a", "extract_b"):
    out_dir = tmp_dir / name
    outputs.append(
        SeqExperiment(
            tmp_dir / "suffix_R1.fastq",
            tmp_dir / "suffix_R2.fastq",
        ).write_umi_fastqs(out_dir, "fixture")
    )
for first, second in zip(outputs[0], outputs[1]):
    assert hashlib.sha256(Path(first).read_bytes()).digest() == \
        hashlib.sha256(Path(second).read_bytes()).digest()
PY

echo "[TEST] Reproduce production UMI edge handling and non-UMI output"
cat > "$tmp_dir/production_edge_R1.fastq" <<'EOF'
@edge1 1:N:0:X
TTGCGCAATGACGTACGTGATTACA
+
IIIIIIIIIIIIIIIIIIIIIIIII
@edge2 1:N:0:X
TTGCGCAATGACGT
+
IIIIIIIIIIIIII
@edge3 1:N:0:X
TTGCGCAATGACGTACGT
+
IIIIIIIIIIIIIIIIII
@edge4 1:N:0:X
GGGGGGGGGGGG
+
IIIIIIIIIIII
EOF
cat > "$tmp_dir/production_edge_R2.fastq" <<'EOF'
@edge1 2:N:0:X
TGCATGCA
+
IIIIIIII
@edge2 2:N:0:X
TGCATGCA
+
IIIIIIII
@edge3 2:N:0:X
TGCATGCA
+
IIIIIIII
@edge4 2:N:0:X
TGCATGCA
+
IIIIIIII
EOF
python3 - "$repo_root" "$tmp_dir" <<'PY'
import gzip
import sys
import types
from pathlib import Path

repo_root = Path(sys.argv[1])
tmp_dir = Path(sys.argv[2])
stub = types.ModuleType("xopen")
stub.xopen = open
sys.modules["xopen"] = stub
sys.path.insert(0, str(repo_root / "01_preprocessing" / "code" / "gene_expression"))
from umi_methods import SeqExperiment

experiment = SeqExperiment(
    tmp_dir / "production_edge_R1.fastq",
    tmp_dir / "production_edge_R2.fastq",
)
outputs = [Path(p) for p in experiment.write_umi_fastqs(tmp_dir / "edge_out", "edge")]

def records(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    assert len(lines) % 4 == 0
    return [lines[i:i + 4] for i in range(0, len(lines), 4)]

umi = records(outputs[0])
non_umi = records(outputs[2])
assert experiment.total_reads == 4
assert experiment.tso_reads == 3
assert experiment.umi_reads == len(umi) == 3
assert experiment.non_umi_reads == len(non_umi) == 1
assert umi[1][0].endswith("|UB:ACGT")  # production retained a short UMI
assert umi[2][1] == ""                 # production retained empty post-UMI R1
assert non_umi[0][0].endswith("|UB:")
PY

echo "[TEST] Pin recovered production command defaults"
grep -Fq 'TRIM_GALORE_OPTS="${TRIM_GALORE_OPTS:---illumina --paired --fastqc}"' \
  "$repo_root/01_preprocessing/code/run_gene_expression.sh"
grep -Fq 'COUNT_READ_PAIRS="${COUNT_READ_PAIRS:-0}"' \
  "$repo_root/01_preprocessing/code/run_gene_expression.sh"
grep -Fq 'UMI_TOOLS_PAIRED="${UMI_TOOLS_PAIRED:-0}"' \
  "$repo_root/01_preprocessing/code/run_gene_expression.sh"
grep -Fq 'STAR_OUT_SAM_ATTRIBUTES="${STAR_OUT_SAM_ATTRIBUTES:-NH HI AS nM XS}"' \
  "$repo_root/01_preprocessing/code/run_gene_expression.sh"
grep -Fq 'FUSION_TWOPASS_MODE="${FUSION_TWOPASS_MODE:-None}"' \
  "$repo_root/02_fusion_analysis/code/targeted/run_star_mapping.sh"
grep -Fq 'star-2.7.2a-0.tar.bz2' "$repo_root/environment-linux-64.lock"

echo "[TEST] Preserve manifest order during matrix aggregation"
mkdir -p "$tmp_dir/counts/cell_b" "$tmp_dir/counts/cell_a"
printf 'gene\tcell\tcount\ngene_a\tcell_b\t2\n' | \
  gzip -c > "$tmp_dir/counts/cell_b/cell_b_counts.tsv.gz"
printf 'gene\tcell\tcount\ngene_b\tcell_a\t3\n' | \
  gzip -c > "$tmp_dir/counts/cell_a/cell_a_counts.tsv.gz"
cat > "$tmp_dir/samples.tsv" <<'EOF'
sample	r1	r2
cell_b	unused_R1.fastq.gz	unused_R2.fastq.gz
cell_a	unused_R1.fastq.gz	unused_R2.fastq.gz
EOF
cat > "$tmp_dir/genes.gtf" <<'EOF'
chr1	fixture	gene	1	100	.	+	.	gene_id "gene_a"; gene_name "GENE_A"; gene_type "test";
chr1	fixture	gene	201	300	.	+	.	gene_id "gene_b"; gene_name "GENE_B"; gene_type "test";
EOF
python3 "$repo_root/01_preprocessing/code/gene_expression/aggregate_umi_counts.py" \
  --counts-root "$tmp_dir/counts" \
  --manifest "$tmp_dir/samples.tsv" \
  --gtf "$tmp_dir/genes.gtf" \
  --out-dir "$tmp_dir/matrix"

first_gzip_sha="$(sha256sum "$tmp_dir/matrix/gene_by_sample_counts.tsv.gz" | awk '{print $1}')"
sleep 1
python3 "$repo_root/01_preprocessing/code/gene_expression/aggregate_umi_counts.py" \
  --counts-root "$tmp_dir/counts" \
  --manifest "$tmp_dir/samples.tsv" \
  --gtf "$tmp_dir/genes.gtf" \
  --out-dir "$tmp_dir/matrix_repeat"
second_gzip_sha="$(sha256sum "$tmp_dir/matrix_repeat/gene_by_sample_counts.tsv.gz" | awk '{print $1}')"
[[ "$first_gzip_sha" == "$second_gzip_sha" ]] || {
  echo "[ERROR] Aggregated gzip output is not byte-deterministic" >&2
  exit 1
}

python3 - "$tmp_dir/matrix" <<'PY'
import sys
import gzip
from pathlib import Path
from scipy.io import mmread

root = Path(sys.argv[1])
assert root.joinpath("barcodes.tsv").read_text() == "cell_b\ncell_a\n"
assert root.joinpath("features.tsv").read_text() == (
    "gene_a\tGENE_A\ttest\n"
    "gene_b\tGENE_B\ttest\n"
)
assert mmread(root / "matrix.mtx").toarray().tolist() == [[2, 0], [0, 3]]
with gzip.open(root / "gene_by_sample_counts.tsv.gz", "rt") as handle:
    assert handle.readline().rstrip("\n") == "gene_id\tcell_b\tcell_a"
PY

echo "[TEST] Reject non-integer UMI counts"
printf 'gene\tcell\tcount\ngene_a\tcell_b\t1.5\n' | \
  gzip -c > "$tmp_dir/counts/cell_b/cell_b_counts.tsv.gz"
if python3 "$repo_root/01_preprocessing/code/gene_expression/aggregate_umi_counts.py" \
    --counts-root "$tmp_dir/counts" \
    --manifest "$tmp_dir/samples.tsv" \
    --gtf "$tmp_dir/genes.gtf" \
    --out-dir "$tmp_dir/invalid_matrix" \
    2> "$tmp_dir/invalid_count.err"; then
  echo "[ERROR] Non-integer UMI count was accepted" >&2
  exit 1
fi
grep -q "finite, non-negative integers" "$tmp_dir/invalid_count.err"

echo "[PASS] Manifest validation and multi-sample aggregation tests passed"

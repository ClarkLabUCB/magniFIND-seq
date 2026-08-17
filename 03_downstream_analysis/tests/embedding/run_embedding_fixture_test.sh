#!/usr/bin/env bash
set -euo pipefail
test_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
analysis_dir="$(cd "$test_dir/../.." && pwd)"
tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT
Rscript "$test_dir/create_fixture_objects.R" "$tmp"
cat > "$tmp/inputs.tsv" <<EOF
dataset	path	format	group	assay	feature_column
PBMC	$tmp/pbmc.rds	seurat_rds	PBMC	RNA	gene_name
K562r	$tmp/k562r.rds	seurat_rds	K562-R	RNA	gene_name
EOF
cat > "$tmp/candidates.tsv" <<'EOF'
sample	expected_BCR_ABL1_candidate
K1	1
K2	0
K3	1
K4	0
K5	0
K6	0
EOF
cat > "$tmp/config.sh" <<EOF
INPUT_MANIFEST="$tmp/inputs.tsv"
FUSION_CANDIDATES="$tmp/candidates.tsv"
OUT_DIR="$tmp/output"
K562_R_GROUP="K562-R"
SEED=20260730
PCA_COMPONENTS=5
JOINT_DIMS=5
K562_DIMS=5
JOINT_VARIABLE_FEATURES=20
K562_VARIABLE_FEATURES=20
EOF
bash "$analysis_dir/code/run_embedding.sh" "$tmp/config.sh"
for path in \
  results/joint_pbmc_k562r.seurat.rds \
  results/k562r_only.seurat.rds \
  results/joint_umap_coordinates.tsv \
  results/k562r_only_umap_coordinates.tsv \
  figures/joint_umap_by_group.pdf \
  figures/joint_umap_fusion_candidates.pdf \
  figures/k562r_only_umap_fusion_candidates.pdf \
  run_metadata/input_file_checksums.tsv \
  run_metadata/analysis_checksums.sha256 \
  run_metadata/source_provenance.tsv \
  sessionInfo.txt; do
  [[ -s "$tmp/output/$path" ]] || { echo "[ERROR] Missing output: $path" >&2; exit 1; }
done
python3 - "$tmp/output" <<'PY'
import csv, sys
from pathlib import Path
out = Path(sys.argv[1])
with (out / "results/joint_umap_coordinates.tsv").open(newline="") as h:
    joint = list(csv.DictReader(h, delimiter="\t"))
with (out / "results/k562r_only_umap_coordinates.tsv").open(newline="") as h:
    k562 = list(csv.DictReader(h, delimiter="\t"))
assert len(joint) == 12 and {r["group"] for r in joint} == {"PBMC", "K562-R"}
assert len(k562) == 6 and sum(int(r["fusion_candidate"]) for r in k562) == 2
assert all(r["group"] == "K562-R" for r in k562)
PY

echo "[TEST] Reproduce identical UMAP coordinate tables"
sed "s|OUT_DIR=\"$tmp/output\"|OUT_DIR=\"$tmp/output_repeat\"|" \
  "$tmp/config.sh" > "$tmp/config_repeat.sh"
bash "$analysis_dir/code/run_embedding.sh" "$tmp/config_repeat.sh"
cmp "$tmp/output/results/joint_umap_coordinates.tsv" \
  "$tmp/output_repeat/results/joint_umap_coordinates.tsv"
cmp "$tmp/output/results/k562r_only_umap_coordinates.tsv" \
  "$tmp/output_repeat/results/k562r_only_umap_coordinates.tsv"

echo "[TEST] Reject a missing K562-r fusion-candidate row"
sed '$d' "$tmp/candidates.tsv" > "$tmp/candidates_missing.tsv"
mkdir "$tmp/missing_candidate_output"
if Rscript "$analysis_dir/code/embedding/run_joint_embedding.R" \
    --input-manifest "$tmp/inputs.tsv" \
    --fusion-candidates "$tmp/candidates_missing.tsv" \
    --out-dir "$tmp/missing_candidate_output" \
    --k562-r-group K562-R --seed 20260730 --pca-components 5 \
    --joint-dims 5 --k562-dims 5 \
    --joint-variable-features 20 --k562-variable-features 20 \
    2> "$tmp/missing_candidate.err"; then
  echo "[ERROR] A missing K562-r candidate row was accepted" >&2
  exit 1
fi
grep -q "missing K562-r samples: K6" "$tmp/missing_candidate.err"
echo "[PASS] Single-cell embedding fixture passed"

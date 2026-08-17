#!/usr/bin/env bash
set -euo pipefail

test_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
analysis_dir="$(cd "$test_dir/../.." && pwd)"
data_dir="$test_dir/data"
expected_dir="$test_dir/expected"

for cmd in Rscript python3; do
  command -v "$cmd" >/dev/null 2>&1 || {
    echo "[ERROR] Required command not found: $cmd" >&2
    exit 2
  }
done

Rscript -e \
  'needed <- c("edgeR", "Matrix", "Seurat", "ggplot2", "ggrepel"); stopifnot(all(vapply(needed, requireNamespace, logical(1), quietly=TRUE)))' \
  >/dev/null

tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

echo "[TEST] Create synthetic two-condition Seurat object"
Rscript "$data_dir/create_fixture_seurat.R" "$tmp_dir/single_cell.rds"

cat > "$tmp_dir/production_bulk_de.csv" <<'EOF'
EnsemblID,gene_name,gene_biotype,logFC,PValue,FDR
gene_a1,GENEA1,protein_coding,-4,0.000001,0.00001
gene_a2,GENEA2,protein_coding,-3,0.000002,0.00002
gene_b1,GENEB1,protein_coding,4,0.000001,0.00001
gene_b2,GENEB2,protein_coding,3,0.000002,0.00002
gene_stable1,GENESTABLE1,protein_coding,0,0.9,0.9
gene_stable2,GENESTABLE2,protein_coding,0,0.9,0.9
gene_mt,MT-FIXTURE,protein_coding,-5,0.0000001,0.000001
gene_rpl,RPL-FIXTURE,protein_coding,5,0.0000001,0.000001
EOF

cat > "$tmp_dir/config.sh" <<EOF
PROJECT_NAME="fixture"
BULK_MATRIX="$data_dir/bulk_expression_matrix.tsv"
BULK_SAMPLE_METADATA="$data_dir/bulk_samples.tsv"
BULK_DE_TABLE="$tmp_dir/production_bulk_de.csv"
BULK_DE_LOGFC_MULTIPLIER=-1
SEURAT_RDS="$tmp_dir/single_cell.rds"
OUT_DIR="$tmp_dir/output"
CONDITION_A="Resistant"
CONDITION_B="Sensitive"
SC_FEATURE_COLUMN="gene_name"
SC_KEEP_COLUMN="analysis_keep"
FDR_THRESHOLD=0.05
MIN_ABS_LOG2FC=1
TOP_N_PER_DIRECTION=2
SIGNATURE_REQUIRE_BULK_DEG=1
REQUIRE_FULL_SIGNATURE=1
MIN_BULK_REPLICATES_PER_CONDITION=2
SC_MIN_TOTAL_COUNTS=1
SC_CPM_PSEUDOCOUNT=1
SC_DE_TEST="wilcox"
SC_DE_LOGFC_THRESHOLD=0
SC_DE_MIN_PCT=0
EXCLUDE_GENE_REGEX='^MT-|^RPL|^RPS'
MODULE_SCORE_SEED=123
MODULE_SCORE_CTRL=2
MODULE_SCORE_NBIN=2
SKIP_COMPLETED=1
EOF

echo "[TEST] Run bulk-to-single-cell concordance workflow"
bash "$analysis_dir/code/run_concordance.sh" "$tmp_dir/config.sh"

python3 - "$tmp_dir/output" "$expected_dir/top_genes.tsv" <<'PY'
import csv
import math
import sys
from pathlib import Path

out = Path(sys.argv[1])
expected_path = Path(sys.argv[2])

required = [
    "results/bulk_edger_all_tested_genes.tsv",
    "results/bulk_edger_deg.tsv",
    "results/production_bulk_de_standardized.tsv",
    "results/single_cell_descriptive_pseudobulk.tsv",
    "results/single_cell_findmarkers.tsv",
    "results/bulk_single_cell_matched_genes.tsv",
    "results/bulk_top_genes.tsv",
    "results/concordance_metrics.tsv",
    "results/single_cell_resistance_scores.tsv",
    "results/resistance_score_group_summary.tsv",
    "figures/bulk_mds.pdf",
    "figures/bulk_bcv.pdf",
    "figures/bulk_volcano.pdf",
    "figures/bulk_signature_single_cell_dotplot.pdf",
    "figures/bulk_single_cell_concordance.pdf",
    "figures/bulk_single_cell_concordance.png",
    "figures/bulk_signature_resistance_score.pdf",
    "run_metadata/analysis_parameters.tsv",
    "run_metadata/bulk_sample_metadata.tsv",
    "run_metadata/single_cell_group_summary.tsv",
    "run_metadata/R_sessionInfo.txt",
    "run_metadata/input_checksums.txt",
    "fixture_concordance.done",
]
for relative in required:
    path = out / relative
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"Missing or empty output: {relative}")

def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

actual = {
    (row["gene"], row["signature"], row["high_condition"])
    for row in read_tsv(out / "results/bulk_top_genes.tsv")
}
expected = {
    (row["gene"], row["signature"], row["high_condition"])
    for row in read_tsv(expected_path)
}
if actual != expected:
    raise SystemExit(f"Selected genes differ: {actual} != {expected}")

metrics = {row["scope"]: row for row in read_tsv(out / "results/concordance_metrics.tsv")}
selected = metrics["selected_genes_combined"]
if int(selected["n_genes"]) != 4:
    raise SystemExit("The combined selected-gene metric must contain four genes")
if not math.isclose(float(selected["direction_concordance"]), 1.0):
    raise SystemExit("Selected genes did not retain exact directional concordance")
if float(selected["pearson_r"]) < 0.95 or float(selected["spearman_rho"]) < 0.8:
    raise SystemExit("Synthetic concordance is unexpectedly weak")

top_rows = read_tsv(out / "results/bulk_top_genes.tsv")
if {row["comparison"] for row in top_rows} != {"Resistant_vs_Sensitive"}:
    raise SystemExit("Comparison direction is incorrect")
if any(row["gene"].startswith(("MT-", "RPL", "RPS")) for row in top_rows):
    raise SystemExit("Excluded prefix was selected")

groups = {row["condition"]: row for row in read_tsv(out / "run_metadata/single_cell_group_summary.tsv")}
if {name: int(row["retained_cells"]) for name, row in groups.items()} != {
    "Resistant": 4,
    "Sensitive": 4,
}:
    raise SystemExit("Single-cell group sizes are incorrect")

score_rows = read_tsv(out / "results/single_cell_resistance_scores.tsv")
if len(score_rows) != 10 or {row["condition"] for row in score_rows} != {
    "Resistant", "Sensitive", "PBMC"
}:
    raise SystemExit("Module scoring did not retain the additional PBMC group")
score_summary = read_tsv(out / "results/resistance_score_group_summary.tsv")
if {row["condition"] for row in score_summary} != {"Resistant", "Sensitive", "PBMC"}:
    raise SystemExit("Resistance-score summary omitted the additional PBMC group")
score_by_group = {}
for row in score_rows:
    score_by_group.setdefault(row["condition"], []).append(float(row["resistance_score"]))
if sum(score_by_group["Resistant"]) / len(score_by_group["Resistant"]) <= \
        sum(score_by_group["Sensitive"]) / len(score_by_group["Sensitive"]):
    raise SystemExit("Resistance score direction is inconsistent with the fixture")
parameters = {
    row["parameter"]: row["value"]
    for row in read_tsv(out / "run_metadata/analysis_parameters.tsv")
}
if parameters["bulk_signature_source"] != "production_bulk_de_table":
    raise SystemExit("The production bulk-DE table was not used for signature ranking")
PY

echo "[TEST] Confirm completed fingerprint skips recomputation"
done_file="$tmp_dir/output/fixture_concordance.done"
mtime_before="$(python3 -c 'import os,sys; print(os.stat(sys.argv[1]).st_mtime_ns)' "$done_file")"
sleep 1
bash "$analysis_dir/code/run_concordance.sh" "$tmp_dir/config.sh"
mtime_after="$(python3 -c 'import os,sys; print(os.stat(sys.argv[1]).st_mtime_ns)' "$done_file")"
[[ "$mtime_before" == "$mtime_after" ]] || {
  echo "[ERROR] Completed output was unexpectedly recomputed" >&2
  exit 1
}

echo "[TEST] Reject a signature shorter than the configured size"
sed \
  -e 's/TOP_N_PER_DIRECTION=2/TOP_N_PER_DIRECTION=100/' \
  -e "s|OUT_DIR=\"$tmp_dir/output\"|OUT_DIR=\"$tmp_dir/short_signature\"|" \
  "$tmp_dir/config.sh" > "$tmp_dir/short_signature_config.sh"
if bash "$analysis_dir/code/run_concordance.sh" "$tmp_dir/short_signature_config.sh" \
    > "$tmp_dir/short_signature.log" 2>&1; then
  echo "[ERROR] An undersized signature was accepted" >&2
  exit 1
fi
grep -q "eligible genes were available" "$tmp_dir/short_signature.log"

echo "[TEST] Reject insufficient bulk replication"
sed \
  -e 's/MIN_BULK_REPLICATES_PER_CONDITION=2/MIN_BULK_REPLICATES_PER_CONDITION=4/' \
  -e "s|OUT_DIR=\"$tmp_dir/output\"|OUT_DIR=\"$tmp_dir/insufficient_replicates\"|" \
  "$tmp_dir/config.sh" > "$tmp_dir/insufficient_replicates_config.sh"
if bash "$analysis_dir/code/run_concordance.sh" "$tmp_dir/insufficient_replicates_config.sh" \
    > "$tmp_dir/insufficient_replicates.log" 2>&1; then
  echo "[ERROR] Insufficient bulk replication was accepted" >&2
  exit 1
fi
grep -q "must contain at least 4 replicates" "$tmp_dir/insufficient_replicates.log"

echo "[PASS] Bulk-to-single-cell concordance fixture passed"

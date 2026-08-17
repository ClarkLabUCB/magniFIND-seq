#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 <config.sh>" >&2
  exit 2
}

[[ $# -eq 1 ]] || usage

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/../.." && pwd)"
config_path="$1"
[[ "$config_path" = /* ]] || config_path="$(pwd)/$config_path"
config_path="$(cd "$(dirname "$config_path")" && pwd)/$(basename "$config_path")"
[[ -f "$config_path" ]] || { echo "[ERROR] Config not found: $config_path" >&2; exit 2; }
config_dir="$(dirname "$config_path")"

# Configuration files are trusted shell input, matching the other workflows.
# shellcheck disable=SC1090
source "$config_path"

required=(PROJECT_NAME BULK_MATRIX BULK_SAMPLE_METADATA SEURAT_RDS OUT_DIR CONDITION_A CONDITION_B)
for name in "${required[@]}"; do
  [[ -n "${!name:-}" ]] || { echo "[ERROR] Missing required setting: $name" >&2; exit 2; }
done
[[ "$PROJECT_NAME" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
  echo "[ERROR] PROJECT_NAME must contain only letters, numbers, dot, underscore, or hyphen" >&2
  exit 2
}

resolve_path() {
  local value="$1"
  if [[ "$value" = /* ]]; then
    printf '%s\n' "$value"
  else
    printf '%s/%s\n' "$config_dir" "$value"
  fi
}

BULK_MATRIX="$(resolve_path "$BULK_MATRIX")"
BULK_SAMPLE_METADATA="$(resolve_path "$BULK_SAMPLE_METADATA")"
SEURAT_RDS="$(resolve_path "$SEURAT_RDS")"
OUT_DIR="$(resolve_path "$OUT_DIR")"

R_SCRIPT="${R_SCRIPT:-$script_dir/concordance/run_bulk_single_cell_concordance.R}"
R_SCRIPT="$(resolve_path "$R_SCRIPT")"
RSCRIPT="${RSCRIPT:-Rscript}"

BULK_COUNT_SUFFIX="${BULK_COUNT_SUFFIX:-_count}"
BULK_GENE_ID_COLUMN="${BULK_GENE_ID_COLUMN:-gene_id}"
BULK_GENE_NAME_COLUMN="${BULK_GENE_NAME_COLUMN:-gene_name}"
BULK_GENE_BIOTYPE_COLUMN="${BULK_GENE_BIOTYPE_COLUMN:-gene_biotype}"
BULK_DE_TABLE="${BULK_DE_TABLE:-}"
BULK_DE_GENE_ID_COLUMN="${BULK_DE_GENE_ID_COLUMN:-EnsemblID}"
BULK_DE_GENE_NAME_COLUMN="${BULK_DE_GENE_NAME_COLUMN:-gene_name}"
BULK_DE_GENE_BIOTYPE_COLUMN="${BULK_DE_GENE_BIOTYPE_COLUMN:-gene_biotype}"
BULK_DE_LOGFC_COLUMN="${BULK_DE_LOGFC_COLUMN:-logFC}"
BULK_DE_PVALUE_COLUMN="${BULK_DE_PVALUE_COLUMN:-PValue}"
BULK_DE_FDR_COLUMN="${BULK_DE_FDR_COLUMN:-FDR}"
BULK_DE_LOGFC_MULTIPLIER="${BULK_DE_LOGFC_MULTIPLIER:--1}"
if [[ -n "$BULK_DE_TABLE" ]]; then
  BULK_DE_TABLE="$(resolve_path "$BULK_DE_TABLE")"
fi
SC_ASSAY="${SC_ASSAY:-RNA}"
SC_COUNTS_LAYER="${SC_COUNTS_LAYER:-counts}"
SC_CONDITION_COLUMN="${SC_CONDITION_COLUMN:-condition}"
SC_FEATURE_COLUMN="${SC_FEATURE_COLUMN:-}"
SC_KEEP_COLUMN="${SC_KEEP_COLUMN:-}"
FDR_THRESHOLD="${FDR_THRESHOLD:-0.05}"
MIN_ABS_LOG2FC="${MIN_ABS_LOG2FC:-1}"
TOP_N_PER_DIRECTION="${TOP_N_PER_DIRECTION:-50}"
SIGNATURE_REQUIRE_BULK_DEG="${SIGNATURE_REQUIRE_BULK_DEG:-1}"
REQUIRE_FULL_SIGNATURE="${REQUIRE_FULL_SIGNATURE:-1}"
MIN_BULK_REPLICATES_PER_CONDITION="${MIN_BULK_REPLICATES_PER_CONDITION:-2}"
SC_MIN_TOTAL_COUNTS="${SC_MIN_TOTAL_COUNTS:-1}"
SC_CPM_PSEUDOCOUNT="${SC_CPM_PSEUDOCOUNT:-1}"
SC_DE_TEST="${SC_DE_TEST:-wilcox}"
SC_DE_LOGFC_THRESHOLD="${SC_DE_LOGFC_THRESHOLD:-0}"
SC_DE_MIN_PCT="${SC_DE_MIN_PCT:-0}"
EXCLUDE_GENE_REGEX="${EXCLUDE_GENE_REGEX:-^MT-|^RP[SL]}"
MODULE_SCORE_SEED="${MODULE_SCORE_SEED:-20260729}"
MODULE_SCORE_CTRL="${MODULE_SCORE_CTRL:-50}"
MODULE_SCORE_NBIN="${MODULE_SCORE_NBIN:-24}"
SKIP_COMPLETED="${SKIP_COMPLETED:-1}"

for path in "$BULK_MATRIX" "$BULK_SAMPLE_METADATA" "$SEURAT_RDS" "$R_SCRIPT"; do
  [[ -s "$path" ]] || { echo "[ERROR] Required input is missing or empty: $path" >&2; exit 2; }
done
if [[ -n "$BULK_DE_TABLE" && ! -s "$BULK_DE_TABLE" ]]; then
  echo "[ERROR] Production bulk differential-expression table is missing: $BULK_DE_TABLE" >&2
  exit 2
fi
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "[ERROR] Rscript not found: $RSCRIPT" >&2; exit 2; }
command -v sha256sum >/dev/null 2>&1 || { echo "[ERROR] sha256sum is required" >&2; exit 2; }
[[ "$SKIP_COMPLETED" == "0" || "$SKIP_COMPLETED" == "1" ]] || {
  echo "[ERROR] SKIP_COMPLETED must be 0 or 1" >&2
  exit 2
}
[[ "$SIGNATURE_REQUIRE_BULK_DEG" == "0" || "$SIGNATURE_REQUIRE_BULK_DEG" == "1" ]] || {
  echo "[ERROR] SIGNATURE_REQUIRE_BULK_DEG must be 0 or 1" >&2
  exit 2
}
[[ "$REQUIRE_FULL_SIGNATURE" == "0" || "$REQUIRE_FULL_SIGNATURE" == "1" ]] || {
  echo "[ERROR] REQUIRE_FULL_SIGNATURE must be 0 or 1" >&2
  exit 2
}
[[ "$MIN_BULK_REPLICATES_PER_CONDITION" =~ ^[1-9][0-9]*$ ]] || {
  echo "[ERROR] MIN_BULK_REPLICATES_PER_CONDITION must be a positive integer" >&2
  exit 2
}

if [[ -z "$OUT_DIR" || "$OUT_DIR" == "/" ]]; then
  echo "[ERROR] Refusing unsafe OUT_DIR: $OUT_DIR" >&2
  exit 2
fi

"$RSCRIPT" -e \
  'needed <- c("edgeR", "Matrix", "Seurat", "ggplot2", "ggrepel"); missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) stop("Missing R packages: ", paste(missing, collapse=", "))' \
  >/dev/null

r_versions="$("$RSCRIPT" -e 'cat(paste(R.version.string, paste0("edgeR=", packageVersion("edgeR")), paste0("Matrix=", packageVersion("Matrix")), paste0("Seurat=", packageVersion("Seurat")), paste0("ggplot2=", packageVersion("ggplot2")), paste0("ggrepel=", packageVersion("ggrepel")), sep="|"))')"
fingerprint_inputs=("$BULK_MATRIX" "$BULK_SAMPLE_METADATA" "$SEURAT_RDS" "$R_SCRIPT" "$config_path")
[[ -z "$BULK_DE_TABLE" ]] || fingerprint_inputs+=("$BULK_DE_TABLE")
fingerprint="$({
  sha256sum "${fingerprint_inputs[@]}"
  printf '%s\n' \
    "$CONDITION_A" "$CONDITION_B" "$BULK_COUNT_SUFFIX" \
    "$BULK_GENE_ID_COLUMN" "$BULK_GENE_NAME_COLUMN" "$BULK_GENE_BIOTYPE_COLUMN" \
    "$BULK_DE_TABLE" "$BULK_DE_GENE_ID_COLUMN" "$BULK_DE_GENE_NAME_COLUMN" \
    "$BULK_DE_GENE_BIOTYPE_COLUMN" "$BULK_DE_LOGFC_COLUMN" \
    "$BULK_DE_PVALUE_COLUMN" "$BULK_DE_FDR_COLUMN" "$BULK_DE_LOGFC_MULTIPLIER" \
    "$SC_ASSAY" "$SC_COUNTS_LAYER" "$SC_CONDITION_COLUMN" "$SC_FEATURE_COLUMN" \
    "$SC_KEEP_COLUMN" "$FDR_THRESHOLD" "$MIN_ABS_LOG2FC" \
    "$TOP_N_PER_DIRECTION" "$SIGNATURE_REQUIRE_BULK_DEG" \
    "$REQUIRE_FULL_SIGNATURE" "$MIN_BULK_REPLICATES_PER_CONDITION" \
    "$SC_MIN_TOTAL_COUNTS" "$SC_CPM_PSEUDOCOUNT" \
    "$SC_DE_TEST" "$SC_DE_LOGFC_THRESHOLD" "$SC_DE_MIN_PCT" \
    "$EXCLUDE_GENE_REGEX" "$MODULE_SCORE_SEED" "$MODULE_SCORE_CTRL" \
    "$MODULE_SCORE_NBIN" "$r_versions"
} | sha256sum | awk '{print $1}')"

done_file="$OUT_DIR/${PROJECT_NAME}_concordance.done"
if [[ "$SKIP_COMPLETED" == "1" && -s "$done_file" ]] && \
   grep -Fxq "fingerprint=$fingerprint" "$done_file" && \
   [[ -s "$OUT_DIR/results/concordance_metrics.tsv" && \
      -s "$OUT_DIR/results/bulk_top_genes.tsv" && \
      -s "$OUT_DIR/results/single_cell_resistance_scores.tsv" && \
      -s "$OUT_DIR/results/resistance_score_group_summary.tsv" && \
      -s "$OUT_DIR/figures/bulk_signature_resistance_score.pdf" && \
      -s "$OUT_DIR/figures/bulk_single_cell_concordance.pdf" ]]; then
  echo "[SKIP] Completed concordance analysis: $OUT_DIR"
  exit 0
fi

mkdir -p "$(dirname "$OUT_DIR")"
tmp_dir="${OUT_DIR}.tmp.$$"
backup_dir="${OUT_DIR}.previous.$$"
rm -rf "$tmp_dir"
[[ ! -e "$backup_dir" ]] || {
  echo "[ERROR] Refusing to replace an existing backup: $backup_dir" >&2
  exit 1
}
trap 'rm -rf "$tmp_dir"' EXIT
mkdir -p "$tmp_dir"

args=(
  --bulk-matrix "$BULK_MATRIX"
  --bulk-sample-metadata "$BULK_SAMPLE_METADATA"
  --seurat-rds "$SEURAT_RDS"
  --out-dir "$tmp_dir"
  --project-name "$PROJECT_NAME"
  --condition-a "$CONDITION_A"
  --condition-b "$CONDITION_B"
  --bulk-count-suffix "$BULK_COUNT_SUFFIX"
  --bulk-gene-id-column "$BULK_GENE_ID_COLUMN"
  --bulk-gene-name-column "$BULK_GENE_NAME_COLUMN"
  --bulk-gene-biotype-column "$BULK_GENE_BIOTYPE_COLUMN"
  --sc-assay "$SC_ASSAY"
  --sc-counts-layer "$SC_COUNTS_LAYER"
  --sc-condition-column "$SC_CONDITION_COLUMN"
  --fdr-threshold "$FDR_THRESHOLD"
  --min-abs-log2fc "$MIN_ABS_LOG2FC"
  --top-n-per-direction "$TOP_N_PER_DIRECTION"
  --signature-require-bulk-deg "$SIGNATURE_REQUIRE_BULK_DEG"
  --require-full-signature "$REQUIRE_FULL_SIGNATURE"
  --min-bulk-replicates-per-condition "$MIN_BULK_REPLICATES_PER_CONDITION"
  --sc-min-total-counts "$SC_MIN_TOTAL_COUNTS"
  --sc-cpm-pseudocount "$SC_CPM_PSEUDOCOUNT"
  --sc-de-test "$SC_DE_TEST"
  --sc-de-logfc-threshold "$SC_DE_LOGFC_THRESHOLD"
  --sc-de-min-pct "$SC_DE_MIN_PCT"
  --exclude-gene-regex "$EXCLUDE_GENE_REGEX"
  --module-score-seed "$MODULE_SCORE_SEED"
  --module-score-ctrl "$MODULE_SCORE_CTRL"
  --module-score-nbin "$MODULE_SCORE_NBIN"
)
[[ -z "$SC_FEATURE_COLUMN" ]] || args+=(--sc-feature-column "$SC_FEATURE_COLUMN")
[[ -z "$SC_KEEP_COLUMN" ]] || args+=(--sc-keep-column "$SC_KEEP_COLUMN")
if [[ -n "$BULK_DE_TABLE" ]]; then
  args+=(
    --bulk-de-table "$BULK_DE_TABLE"
    --bulk-de-gene-id-column "$BULK_DE_GENE_ID_COLUMN"
    --bulk-de-gene-name-column "$BULK_DE_GENE_NAME_COLUMN"
    --bulk-de-gene-biotype-column "$BULK_DE_GENE_BIOTYPE_COLUMN"
    --bulk-de-logfc-column "$BULK_DE_LOGFC_COLUMN"
    --bulk-de-pvalue-column "$BULK_DE_PVALUE_COLUMN"
    --bulk-de-fdr-column "$BULK_DE_FDR_COLUMN"
    --bulk-de-logfc-multiplier "$BULK_DE_LOGFC_MULTIPLIER"
  )
fi

echo "[RUN] Bulk and single-cell concordance: $PROJECT_NAME"
"$RSCRIPT" "$R_SCRIPT" "${args[@]}"

mkdir -p "$tmp_dir/run_metadata"
cp "$config_path" "$tmp_dir/run_metadata/config.sh"
bash "$repo_root/scripts/write_source_provenance.sh" \
  "$repo_root" "$tmp_dir/run_metadata/source_provenance.tsv"
{
  echo "fingerprint=$fingerprint"
  sha256sum "${fingerprint_inputs[@]}"
  printf 'r_packages\t%s\n' "$r_versions"
} > "$tmp_dir/run_metadata/input_checksums.txt"
{
  echo "fingerprint=$fingerprint"
  echo "completed_at_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$tmp_dir/${PROJECT_NAME}_concordance.done"

if [[ -e "$OUT_DIR" ]]; then
  mv "$OUT_DIR" "$backup_dir"
fi
if ! mv "$tmp_dir" "$OUT_DIR"; then
  if [[ -e "$backup_dir" ]] && ! mv "$backup_dir" "$OUT_DIR"; then
    echo "[ERROR] Previous output remains at: $backup_dir" >&2
  fi
  echo "[ERROR] Could not promote completed output" >&2
  exit 1
fi
rm -rf "$backup_dir"
trap - EXIT

echo "[DONE] Results: $OUT_DIR"

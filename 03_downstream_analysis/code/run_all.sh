#!/usr/bin/env bash
set -euo pipefail

[[ $# -eq 1 ]] || { echo "Usage: $0 <config.sh>" >&2; exit 2; }
code_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$code_dir/../.." && pwd)"
config="$1"; [[ "$config" = /* ]] || config="$(pwd)/$config"
config="$(cd "$(dirname "$config")" && pwd)/$(basename "$config")"
config_dir="$(dirname "$config")"
# shellcheck disable=SC1090
source "$config"

resolve_optional() {
  local value="$1"
  [[ -z "$value" ]] && return 0
  [[ "$value" = /* ]] && printf '%s\n' "$value" || printf '%s/%s\n' "$config_dir" "$value"
}
run_configured() {
  local label="$1" runner="$2" configured="$3"
  [[ -n "$configured" ]] || { echo "[SKIP] $label (no config)"; return 0; }
  local resolved
  resolved="$(resolve_optional "$configured")"
  [[ -s "$resolved" ]] || { echo "[ERROR] $label config not found: $resolved" >&2; exit 2; }
  echo "[RUN] $label"
  bash "$runner" "$resolved"
}

configured_count=0
for name in GENE_EXPRESSION_CONFIG RAW_QC_CONFIG TARGETED_FUSION_CONFIG \
  PUBLIC_K562_CONFIG SINGLE_CELL_EMBEDDING_CONFIG BULK_CONCORDANCE_CONFIG; do
  [[ -z "${!name:-}" ]] || configured_count=$((configured_count + 1))
done
[[ "$configured_count" -gt 0 ]] || { echo "[ERROR] No paper workflow is configured" >&2; exit 2; }

run_configured "gene-expression quantification" \
  "$repo_root/01_preprocessing/code/run_gene_expression.sh" \
  "${GENE_EXPRESSION_CONFIG:-}"
run_configured "single-bead raw QC" \
  "$repo_root/01_preprocessing/code/run_raw_qc.sh" "${RAW_QC_CONFIG:-}"

if [[ -n "${TARGETED_FUSION_CONFIG:-}" ]]; then
  targeted="$(resolve_optional "$TARGETED_FUSION_CONFIG")"
  [[ -s "$targeted" ]] || { echo "[ERROR] Targeted fusion config not found: $targeted" >&2; exit 2; }
  echo "[RUN] targeted fusion mapping and candidate scan"
  bash "$repo_root/02_fusion_analysis/code/run_targeted_mapping.sh" "$targeted"
  bash "$repo_root/02_fusion_analysis/code/run_targeted_scan.sh" "$targeted"
else
  echo "[SKIP] targeted fusion mapping and candidate scan (no config)"
fi

run_configured "public K562 chimeric quantification" \
  "$repo_root/02_fusion_analysis/code/run_public_k562.sh" "${PUBLIC_K562_CONFIG:-}"
run_configured "PBMC/K562-r embeddings" \
  "$repo_root/03_downstream_analysis/code/run_embedding.sh" \
  "${SINGLE_CELL_EMBEDDING_CONFIG:-}"
run_configured "bulk/single-cell concordance and resistance scores" \
  "$repo_root/03_downstream_analysis/code/run_concordance.sh" \
  "${BULK_CONCORDANCE_CONFIG:-}"

echo "[DONE] Configured manuscript workflows completed"

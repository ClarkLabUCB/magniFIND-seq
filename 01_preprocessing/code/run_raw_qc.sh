#!/usr/bin/env bash
set -euo pipefail

[[ $# -eq 1 ]] || { echo "Usage: $0 <config.sh>" >&2; exit 2; }
code_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$code_dir/../.." && pwd)"
config="$1"
[[ "$config" = /* ]] || config="$(pwd)/$config"
config="$(cd "$(dirname "$config")" && pwd)/$(basename "$config")"
config_dir="$(dirname "$config")"
# shellcheck disable=SC1090
source "$config"

: "${PIPELINE_OUTPUT:?Missing PIPELINE_OUTPUT}"
: "${OUT_DIR:?Missing OUT_DIR}"
resolve_path() { [[ "$1" = /* ]] && printf '%s\n' "$1" || printf '%s/%s\n' "$config_dir" "$1"; }
PIPELINE_OUTPUT="$(resolve_path "$PIPELINE_OUTPUT")"
OUT_DIR="$(resolve_path "$OUT_DIR")"
PYTHON="${PYTHON:-python3}"
ASSIGNED_READ_THRESHOLD="${ASSIGNED_READ_THRESHOLD:-478877}"
DEDUPLICATED_UMI_THRESHOLD="${DEDUPLICATED_UMI_THRESHOLD:-12149}"
ASSIGNED_READ_PERCENTILE="${ASSIGNED_READ_PERCENTILE:-0.65}"
DEDUPLICATED_UMI_PERCENTILE="${DEDUPLICATED_UMI_PERCENTILE:-0.70}"
[[ -s "$PIPELINE_OUTPUT/run_metadata/samples.normalized.tsv" ]] || {
  echo "[ERROR] Gene-pipeline manifest not found under: $PIPELINE_OUTPUT" >&2
  exit 2
}
[[ -n "$OUT_DIR" && "$OUT_DIR" != / ]] || { echo "[ERROR] Unsafe OUT_DIR" >&2; exit 2; }

tmp="${OUT_DIR}.tmp.$$"
backup="${OUT_DIR}.previous.$$"
trap 'rm -rf "$tmp"' EXIT
[[ ! -e "$tmp" && ! -e "$backup" ]] || {
  echo "[ERROR] Temporary or backup path already exists" >&2
  exit 1
}
"$PYTHON" "$code_dir/raw_qc/summarize_raw_qc.py" \
  --pipeline-output "$PIPELINE_OUTPUT" \
  --out-dir "$tmp" \
  --assigned-read-threshold "$ASSIGNED_READ_THRESHOLD" \
  --deduplicated-umi-threshold "$DEDUPLICATED_UMI_THRESHOLD" \
  --assigned-read-percentile "$ASSIGNED_READ_PERCENTILE" \
  --deduplicated-umi-percentile "$DEDUPLICATED_UMI_PERCENTILE"
"$PYTHON" "$code_dir/raw_qc/summarize_production_read_tracking.py" \
  --pipeline-output "$PIPELINE_OUTPUT" \
  --manifest "$PIPELINE_OUTPUT/run_metadata/samples.normalized.tsv" \
  --out-dir "$tmp"
mkdir -p "$tmp/run_metadata"
cp "$config" "$tmp/run_metadata/config.sh"
bash "$repo_root/scripts/write_source_provenance.sh" \
  "$repo_root" "$tmp/run_metadata/source_provenance.tsv"
if [[ -e "$OUT_DIR" ]]; then mv "$OUT_DIR" "$backup"; fi
if ! mv "$tmp" "$OUT_DIR"; then
  [[ ! -e "$backup" ]] || mv "$backup" "$OUT_DIR"
  exit 1
fi
rm -rf "$backup"
trap - EXIT
echo "[DONE] Raw QC: $OUT_DIR"

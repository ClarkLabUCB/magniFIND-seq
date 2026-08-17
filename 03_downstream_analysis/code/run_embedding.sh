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
for name in INPUT_MANIFEST FUSION_CANDIDATES OUT_DIR K562_R_GROUP; do
  [[ -n "${!name:-}" ]] || { echo "[ERROR] Missing setting: $name" >&2; exit 2; }
done
resolve_path() { [[ "$1" = /* ]] && printf '%s\n' "$1" || printf '%s/%s\n' "$config_dir" "$1"; }
INPUT_MANIFEST="$(resolve_path "$INPUT_MANIFEST")"
FUSION_CANDIDATES="$(resolve_path "$FUSION_CANDIDATES")"
OUT_DIR="$(resolve_path "$OUT_DIR")"
RSCRIPT="${RSCRIPT:-Rscript}"; PYTHON="${PYTHON:-python3}"
SEED="${SEED:-20260730}"
PCA_COMPONENTS="${PCA_COMPONENTS:-30}"
JOINT_DIMS="${JOINT_DIMS:-20}"
K562_DIMS="${K562_DIMS:-15}"
JOINT_VARIABLE_FEATURES="${JOINT_VARIABLE_FEATURES:-3000}"
K562_VARIABLE_FEATURES="${K562_VARIABLE_FEATURES:-2000}"
for path in "$INPUT_MANIFEST" "$FUSION_CANDIDATES"; do [[ -s "$path" ]] || { echo "[ERROR] Missing input: $path" >&2; exit 2; }; done
[[ -n "$OUT_DIR" && "$OUT_DIR" != / ]] || { echo "[ERROR] Unsafe OUT_DIR" >&2; exit 2; }
tmp="${OUT_DIR}.tmp.$$"; backup="${OUT_DIR}.previous.$$"
trap 'rm -rf "$tmp"' EXIT
[[ ! -e "$tmp" && ! -e "$backup" ]] || { echo "[ERROR] Temporary or backup path exists" >&2; exit 1; }
mkdir -p "$tmp"
"$RSCRIPT" "$code_dir/embedding/run_joint_embedding.R" \
  --input-manifest "$INPUT_MANIFEST" --fusion-candidates "$FUSION_CANDIDATES" \
  --out-dir "$tmp" --k562-r-group "$K562_R_GROUP" --seed "$SEED" \
  --pca-components "$PCA_COMPONENTS" \
  --joint-dims "$JOINT_DIMS" --k562-dims "$K562_DIMS" \
  --joint-variable-features "$JOINT_VARIABLE_FEATURES" \
  --k562-variable-features "$K562_VARIABLE_FEATURES"
cp "$config" "$tmp/config.sh"
mkdir -p "$tmp/run_metadata"
cp "$INPUT_MANIFEST" "$tmp/run_metadata/inputs.tsv"
cp "$FUSION_CANDIDATES" "$tmp/run_metadata/fusion_candidates.tsv"
"$PYTHON" "$repo_root/scripts/hash_manifest_inputs.py" \
  --manifest "$INPUT_MANIFEST" \
  --output "$tmp/run_metadata/input_file_checksums.tsv"
sha256sum "$INPUT_MANIFEST" "$FUSION_CANDIDATES" \
  "$code_dir/embedding/run_joint_embedding.R" > \
  "$tmp/run_metadata/analysis_checksums.sha256"
bash "$repo_root/scripts/write_source_provenance.sh" \
  "$repo_root" "$tmp/run_metadata/source_provenance.tsv"
[[ ! -e "$backup" ]] || { echo "[ERROR] Existing backup: $backup" >&2; exit 1; }
if [[ -e "$OUT_DIR" ]]; then mv "$OUT_DIR" "$backup"; fi
if ! mv "$tmp" "$OUT_DIR"; then [[ ! -e "$backup" ]] || mv "$backup" "$OUT_DIR"; exit 1; fi
rm -rf "$backup"; trap - EXIT
echo "[DONE] Embedding: $OUT_DIR"

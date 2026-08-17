#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: bash validate_inputs.sh <config.sh>" >&2
  exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/../../.." && pwd)"
config_file="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
config_dir="$(dirname "$config_file")"

# shellcheck disable=SC1090
source "$config_file"

: "${SAMPLES_TSV:?SAMPLES_TSV is required}"
: "${OUT_DIR:?OUT_DIR is required}"

if [[ "$SAMPLES_TSV" != /* ]]; then
  SAMPLES_TSV="$config_dir/$SAMPLES_TSV"
fi
if [[ "$OUT_DIR" != /* ]]; then
  OUT_DIR="$config_dir/$OUT_DIR"
fi
PYTHON="${PYTHON:-python3}"
if [[ "$PYTHON" == */* && "$PYTHON" != /* ]]; then
  PYTHON="$config_dir/$PYTHON"
fi

mkdir -p "$OUT_DIR/run_metadata"
normalized="$OUT_DIR/run_metadata/fusion_samples.normalized.tsv"
tmp_normalized="${normalized}.tmp.$$"

"$PYTHON" "$repo_root/scripts/validate_sample_manifest.py" \
  --manifest "$SAMPLES_TSV" \
  --output "$tmp_normalized" \
  --check-fastq
mv "$tmp_normalized" "$normalized"

echo "[OK] Fusion mapping inputs are valid"
echo "[OK] Normalized manifest: $normalized"

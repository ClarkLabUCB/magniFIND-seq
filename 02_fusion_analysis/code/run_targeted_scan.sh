#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: bash 02_fusion_analysis/code/run_targeted_scan.sh <config.sh>" >&2
  exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
config_dir="$(dirname "$config_file")"

# shellcheck disable=SC1090
source "$config_file"
PYTHON="${PYTHON:-python3}"
if [[ "$PYTHON" == */* ]]; then
  [[ "$PYTHON" = /* ]] || PYTHON="$config_dir/$PYTHON"
  [[ -x "$PYTHON" ]] || { echo "[ERROR] Python executable not found: $PYTHON" >&2; exit 1; }
else
  PYTHON="$(command -v "$PYTHON")"
fi

"$PYTHON" "$script_dir/targeted/scan_bcr_abl1_chimeric_junctions.py" --config "$config_file"

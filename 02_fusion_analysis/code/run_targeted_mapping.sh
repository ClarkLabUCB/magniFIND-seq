#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: bash 02_fusion_analysis/code/run_targeted_mapping.sh <config.sh>" >&2
  exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/../.." && pwd)"
config_file="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
config_dir="$(dirname "$config_file")"

# shellcheck disable=SC1090
source "$config_file"
: "${OUT_DIR:?OUT_DIR is required}"
: "${STAR:?STAR is required}"
: "${STAR_INDEX:?STAR_INDEX is required}"
if [[ "$OUT_DIR" != /* ]]; then
  OUT_DIR="$config_dir/$OUT_DIR"
fi
if [[ "$STAR_INDEX" != /* ]]; then
  STAR_INDEX="$config_dir/$STAR_INDEX"
fi
if [[ -z "$OUT_DIR" || "$OUT_DIR" == "/" ]]; then
  echo "[ERROR] Refusing unsafe OUT_DIR: $OUT_DIR" >&2
  exit 1
fi

SAMTOOLS="${SAMTOOLS:-samtools}"
if [[ "$STAR" == */* ]]; then
  [[ "$STAR" = /* ]] || STAR="$config_dir/$STAR"
else
  STAR="$(command -v "$STAR")"
fi
if [[ "$SAMTOOLS" == */* ]]; then
  [[ "$SAMTOOLS" = /* ]] || SAMTOOLS="$config_dir/$SAMTOOLS"
else
  SAMTOOLS="$(command -v "$SAMTOOLS")"
fi

bash "$script_dir/targeted/validate_inputs.sh" "$config_file"
normalized_manifest="$OUT_DIR/run_metadata/fusion_samples.normalized.tsv"

index_checksums=()
for index_file in genomeParameters.txt Genome SA SAindex chrName.txt; do
  index_path="$STAR_INDEX/$index_file"
  [[ -s "$index_path" ]] || {
    echo "[ERROR] STAR index file not found: $index_path" >&2
    exit 1
  }
  index_checksums+=("$index_file=$(sha256sum "$index_path" | awk '{print $1}')")
done
FUSION_STAR_INDEX_SIGNATURE="$({
  printf '%s\n' "${index_checksums[@]}"
} | sha256sum | awk '{print $1}')"
export FUSION_STAR_INDEX_SIGNATURE

versions_file="$OUT_DIR/run_metadata/fusion_tool_versions.txt"
{
  printf 'project\t%s\n' "${PROJECT_NAME:-magnifind_seq}"
  printf 'config\t%s\n' "$config_file"
  printf 'manifest\t%s\n' "$normalized_manifest"
  printf 'star_index\t%s\n' "$STAR_INDEX"
  printf 'STAR\t%s\n' "$("$STAR" --version 2>&1 | head -n 1)"
  printf 'samtools\t%s\n' "$("$SAMTOOLS" --version 2>&1 | head -n 1)"
  printf 'threads\t%s\n' "${THREADS:-8}"
  printf 'twopass_mode\t%s\n' "${FUSION_TWOPASS_MODE:-None}"
  printf 'chim_segment_min\t%s\n' "${CHIM_SEGMENT_MIN:-12}"
  printf 'chim_junction_overhang_min\t%s\n' "${CHIM_JUNCTION_OVERHANG_MIN:-12}"
  printf 'chim_main_segment_mult_nmax\t%s\n' "${CHIM_MAIN_SEGMENT_MULT_NMAX:-1}"
  printf 'tag_bam_with_cb_ub\t%s\n' "${TAG_BAM_WITH_CB_UB:-1}"
  printf 'predecompress_fastq\t%s\n' "${PREDECOMPRESS_FASTQ:-0}"
  printf 'star_index_sha256\t%s\n' "$FUSION_STAR_INDEX_SIGNATURE"
} > "${versions_file}.tmp.$$"
mv "${versions_file}.tmp.$$" "$versions_file"
cp "$config_file" "$OUT_DIR/run_metadata/config.sh"
bash "$repo_root/scripts/write_source_provenance.sh" \
  "$repo_root" "$OUT_DIR/run_metadata/source_provenance.tsv"

while IFS=$'\t' read -r sample r1 r2; do
  [[ "$sample" == "sample" || -z "${sample:-}" ]] && continue
  bash "$script_dir/targeted/run_star_mapping.sh" "$config_file" "$sample" "$r1" "$r2"
done < "$normalized_manifest"

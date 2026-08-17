#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: bash run_star_mapping.sh <config.sh> <sample> <R1.fq.gz> <R2.fq.gz>" >&2
}

if [[ $# -ne 4 ]]; then
  usage
  exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/../../.." && pwd)"
config_file="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
config_dir="$(dirname "$config_file")"
sample="$2"
r1="$3"
r2="$4"

# shellcheck disable=SC1090
source "$config_file"

: "${OUT_DIR:?OUT_DIR is required}"
: "${STAR:?STAR is required}"
: "${STAR_INDEX:?STAR_INDEX is required}"

THREADS="${THREADS:-8}"
FUSION_TWOPASS_MODE="${FUSION_TWOPASS_MODE:-None}"
CHIM_SEGMENT_MIN="${CHIM_SEGMENT_MIN:-12}"
CHIM_JUNCTION_OVERHANG_MIN="${CHIM_JUNCTION_OVERHANG_MIN:-12}"
CHIM_MAIN_SEGMENT_MULT_NMAX="${CHIM_MAIN_SEGMENT_MULT_NMAX:-1}"
TAG_BAM_WITH_CB_UB="${TAG_BAM_WITH_CB_UB:-1}"
REQUIRE_CB_UB_TAGS="${REQUIRE_CB_UB_TAGS:-1}"
PYTHON="${PYTHON:-python3}"
SAMTOOLS="${SAMTOOLS:-samtools}"
READ_FILES_COMMAND="${READ_FILES_COMMAND:-gzip -cd}"
PREDECOMPRESS_FASTQ="${PREDECOMPRESS_FASTQ:-0}"

resolve_config_path() {
  if [[ "$1" = /* ]]; then
    printf '%s\n' "$1"
  else
    printf '%s\n' "$config_dir/$1"
  fi
}

resolve_executable() {
  local value="$1"
  if [[ "$value" == */* ]]; then
    value="$(resolve_config_path "$value")"
    [[ -x "$value" ]] || { echo "[ERROR] Executable not found: $value" >&2; exit 1; }
    printf '%s\n' "$value"
  else
    command -v "$value" || { echo "[ERROR] Command not found: $value" >&2; exit 1; }
  fi
}

OUT_DIR="$(resolve_config_path "$OUT_DIR")"
STAR_INDEX="$(resolve_config_path "$STAR_INDEX")"
STAR="$(resolve_executable "$STAR")"
PYTHON="$(resolve_executable "$PYTHON")"
SAMTOOLS="$(resolve_executable "$SAMTOOLS")"
# shellcheck disable=SC2206
read_files_command=($READ_FILES_COMMAND)
[[ "${#read_files_command[@]}" -gt 0 ]] || {
  echo "[ERROR] READ_FILES_COMMAND must not be empty" >&2
  exit 1
}
read_files_command[0]="$(resolve_executable "${read_files_command[0]}")"

if [[ -z "$OUT_DIR" || "$OUT_DIR" == "/" ]]; then
  echo "[ERROR] Refusing unsafe OUT_DIR: $OUT_DIR" >&2
  exit 1
fi
if [[ ! "$sample" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
  echo "[ERROR] Unsafe sample ID: $sample" >&2
  exit 1
fi
if [[ ! -s "$r1" || ! -s "$r2" ]]; then
  echo "[ERROR] Missing FASTQ pair for $sample" >&2
  exit 1
fi
if [[ "$TAG_BAM_WITH_CB_UB" != "0" && "$TAG_BAM_WITH_CB_UB" != "1" ]]; then
  echo "[ERROR] TAG_BAM_WITH_CB_UB must be 0 or 1" >&2
  exit 1
fi
if [[ "$PREDECOMPRESS_FASTQ" != "0" && "$PREDECOMPRESS_FASTQ" != "1" ]]; then
  echo "[ERROR] PREDECOMPRESS_FASTQ must be 0 or 1" >&2
  exit 1
fi
if [[ "$REQUIRE_CB_UB_TAGS" != "0" && "$REQUIRE_CB_UB_TAGS" != "1" ]]; then
  echo "[ERROR] REQUIRE_CB_UB_TAGS must be 0 or 1" >&2
  exit 1
fi
if [[ "$FUSION_TWOPASS_MODE" != "Basic" && "$FUSION_TWOPASS_MODE" != "None" ]]; then
  echo "[ERROR] FUSION_TWOPASS_MODE must be Basic or None" >&2
  exit 1
fi
for numeric_setting in THREADS CHIM_SEGMENT_MIN CHIM_JUNCTION_OVERHANG_MIN CHIM_MAIN_SEGMENT_MULT_NMAX; do
  value="${!numeric_setting}"
  if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
    echo "[ERROR] $numeric_setting must be a positive integer" >&2
    exit 1
  fi
done

content_identity() {
  sha256sum "$1" | awk '{print $1}'
}

index_identities=()
for index_file in genomeParameters.txt Genome SA SAindex chrName.txt; do
  index_path="$STAR_INDEX/$index_file"
  [[ -s "$index_path" ]] || {
    echo "[ERROR] STAR index file not found: $index_path" >&2
    exit 1
  }
  if [[ -z "${FUSION_STAR_INDEX_SIGNATURE:-}" ]]; then
    index_identities+=("$index_file=$(content_identity "$index_path")")
  fi
done
if [[ -n "${FUSION_STAR_INDEX_SIGNATURE:-}" ]]; then
  star_index_signature="$FUSION_STAR_INDEX_SIGNATURE"
else
  star_index_signature="$(
    "$PYTHON" -c 'import hashlib,sys; print(hashlib.sha256("\0".join(sys.argv[1:]).encode()).hexdigest())' \
      "${index_identities[@]}"
  )"
fi

tagger="$repo_root/scripts/tag_bam_with_cb_umi.py"
[[ -s "$tagger" ]] || { echo "[ERROR] Shared BAM tagger not found: $tagger" >&2; exit 1; }

fusion_root="$OUT_DIR/fusion_mapped"
sample_dir="$fusion_root/$sample"
prefix="$sample_dir/${sample}_"
bam_out="${prefix}Aligned.sortedByCoord.out.bam"
chim_junction="${prefix}Chimeric.out.junction"
log_final="${prefix}Log.final.out"
log_out="${prefix}Log.out"
log_progress="${prefix}Log.progress.out"
sj_out="${prefix}SJ.out.tab"
tagged_bam="${prefix}fusion_tagged.bam"
done_file="${prefix}fusion_mapping.done"
mkdir -p "$fusion_root" "$OUT_DIR/run_metadata"

signature="$(
  "$PYTHON" -c 'import hashlib,sys; print(hashlib.sha256("\0".join(sys.argv[1:]).encode()).hexdigest())' \
    "R1=$(content_identity "$r1")" \
    "R2=$(content_identity "$r2")" \
    "STAR_INDEX=$star_index_signature" \
    "RUNNER=$(content_identity "$script_dir/run_star_mapping.sh")" \
    "TAGGER=$(content_identity "$tagger")" \
    "CONFIG=$(content_identity "$config_file")" \
    "STAR=$("$STAR" --version 2>&1 | head -n 1)" \
    "SAMTOOLS=$("$SAMTOOLS" --version 2>&1 | head -n 1)" \
    "TWOPASS=$FUSION_TWOPASS_MODE" \
    "CHIM_SEGMENT_MIN=$CHIM_SEGMENT_MIN" \
    "CHIM_OVERHANG_MIN=$CHIM_JUNCTION_OVERHANG_MIN" \
    "CHIM_MULT_NMAX=$CHIM_MAIN_SEGMENT_MULT_NMAX" \
    "READ_FILES_COMMAND=${read_files_command[*]}" \
    "PREDECOMPRESS_FASTQ=$PREDECOMPRESS_FASTQ" \
    "TAG=$TAG_BAM_WITH_CB_UB" \
    "REQUIRE_TAGS=$REQUIRE_CB_UB_TAGS"
)"

mapping_complete() {
  [[ -s "$done_file" ]] || return 1
  [[ "$(<"$done_file")" == "$signature" ]] || return 1
  [[ -s "$bam_out" && -s "${bam_out}.bai" ]] || return 1
  [[ -e "$chim_junction" && -e "$sj_out" ]] || return 1
  [[ -s "$log_final" && -s "$log_out" && -s "$log_progress" ]] || return 1
  "$SAMTOOLS" quickcheck "$bam_out" >/dev/null 2>&1 || return 1
  if [[ "$TAG_BAM_WITH_CB_UB" == "1" ]]; then
    [[ -s "$tagged_bam" && -s "${tagged_bam}.bai" ]] || return 1
    "$SAMTOOLS" quickcheck "$tagged_bam" >/dev/null 2>&1 || return 1
  fi
}

if mapping_complete; then
  echo "[SKIP] Completed fusion mapping exists for $sample"
  exit 0
fi

# Remove only abandoned hidden work directories older than one day. Current
# and previously completed sample directories are left untouched.
find "$fusion_root" -maxdepth 1 -type d -name ".${sample}.work.*" -mtime +1 \
  -exec rm -rf {} +

work_dir="$(mktemp -d "$fusion_root/.${sample}.work.XXXXXX")"
cleanup() {
  if [[ -n "${work_dir:-}" && -d "$work_dir" ]]; then
    rm -rf "$work_dir"
  fi
}
trap cleanup EXIT
work_prefix="$work_dir/${sample}_"

twopass_args=()
if [[ "$FUSION_TWOPASS_MODE" != "None" ]]; then
  twopass_args=(--twopassMode "$FUSION_TWOPASS_MODE")
fi
echo "[STAR-FUSION] $sample"
star_r1="$r1"
star_r2="$r2"
if [[ "$PREDECOMPRESS_FASTQ" == "1" ]]; then
  star_r1="$work_dir/.star_input_R1.fastq"
  star_r2="$work_dir/.star_input_R2.fastq"
  "${read_files_command[@]}" "$r1" > "$star_r1"
  "${read_files_command[@]}" "$r2" > "$star_r2"
fi
star_command=(
  "$STAR"
  --genomeDir "$STAR_INDEX"
  --readFilesIn "$star_r1" "$star_r2"
)
if [[ "$PREDECOMPRESS_FASTQ" == "0" ]]; then
  star_command+=(--readFilesCommand "${read_files_command[@]}")
fi
star_command+=(
  --outFileNamePrefix "$work_prefix"
  --outSAMtype BAM SortedByCoordinate
  --runThreadN "$THREADS"
)
if [[ "$FUSION_TWOPASS_MODE" != "None" ]]; then
  star_command+=("${twopass_args[@]}")
fi
star_command+=(
  --chimSegmentMin "$CHIM_SEGMENT_MIN"
  --chimJunctionOverhangMin "$CHIM_JUNCTION_OVERHANG_MIN"
  --chimOutType Junctions WithinBAM SoftClip
  --chimOutJunctionFormat 1
  --chimMainSegmentMultNmax "$CHIM_MAIN_SEGMENT_MULT_NMAX"
  --outSAMunmapped Within KeepPairs
  --outSAMattributes NH HI AS nM NM MD XS ch
  --outSAMattrRGline "ID:$sample" "SM:$sample"
)
"${star_command[@]}"
rm -f "$work_dir/.star_input_R1.fastq" "$work_dir/.star_input_R2.fastq"

work_bam="${work_prefix}Aligned.sortedByCoord.out.bam"
work_junction="${work_prefix}Chimeric.out.junction"
work_sj="${work_prefix}SJ.out.tab"
work_log_final="${work_prefix}Log.final.out"
work_log_out="${work_prefix}Log.out"
work_log_progress="${work_prefix}Log.progress.out"
test -s "$work_bam"
test -e "$work_junction"
test -e "$work_sj"
test -s "$work_log_final"
test -s "$work_log_out"
test -s "$work_log_progress"
"$SAMTOOLS" quickcheck "$work_bam"
"$SAMTOOLS" index -@ "$THREADS" "$work_bam"

if [[ "$TAG_BAM_WITH_CB_UB" == "1" ]]; then
  echo "[TAG] Adding CB/UB tags without filtering alignments"
  work_tagged="${work_prefix}fusion_tagged.bam"
  tag_args=(--input "$work_bam" --output "$work_tagged")
  if [[ "$REQUIRE_CB_UB_TAGS" == "1" ]]; then
    tag_args+=(--require-all-tags)
  fi
  "$PYTHON" "$tagger" "${tag_args[@]}"
  "$SAMTOOLS" quickcheck "$work_tagged"
  source_records="$("$SAMTOOLS" view -c "$work_bam")"
  tagged_records="$("$SAMTOOLS" view -c "$work_tagged")"
  if [[ "$source_records" != "$tagged_records" ]]; then
    echo "[ERROR] BAM tagging changed alignment count for $sample" >&2
    exit 1
  fi
  "$SAMTOOLS" index -@ "$THREADS" "$work_tagged"
fi

printf '%s\n' "$signature" > "${work_prefix}fusion_mapping.done"

backup_dir="$fusion_root/.${sample}.previous.$$"
[[ ! -e "$backup_dir" ]] || {
  echo "[ERROR] Refusing to replace existing backup: $backup_dir" >&2
  exit 1
}
if [[ -e "$sample_dir" ]]; then
  mv "$sample_dir" "$backup_dir"
fi
if ! mv "$work_dir" "$sample_dir"; then
  if [[ -e "$backup_dir" ]] && ! mv "$backup_dir" "$sample_dir"; then
    echo "[ERROR] Previous output remains at: $backup_dir" >&2
  fi
  echo "[ERROR] Could not promote completed fusion output for $sample" >&2
  exit 1
fi
work_dir=""
if ! mapping_complete; then
  echo "[ERROR] Promoted fusion output failed validation for $sample" >&2
  rm -rf "$sample_dir"
  if [[ -e "$backup_dir" ]]; then
    mv "$backup_dir" "$sample_dir"
  fi
  exit 1
fi
rm -rf "$backup_dir"
trap - EXIT
echo "[DONE] $sample"

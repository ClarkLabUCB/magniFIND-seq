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
for name in INPUT_BAM CALLED_BARCODES GTF_FILE STAR_INDEX OUT_DIR; do
  [[ -n "${!name:-}" ]] || { echo "[ERROR] Missing setting: $name" >&2; exit 2; }
done
resolve_path() { [[ "$1" = /* ]] && printf '%s\n' "$1" || printf '%s/%s\n' "$config_dir" "$1"; }
INPUT_BAM="$(resolve_path "$INPUT_BAM")"
CALLED_BARCODES="$(resolve_path "$CALLED_BARCODES")"
GTF_FILE="$(resolve_path "$GTF_FILE")"
STAR_INDEX="$(resolve_path "$STAR_INDEX")"
OUT_DIR="$(resolve_path "$OUT_DIR")"
STAR="${STAR:-STAR}"
SAMTOOLS="${SAMTOOLS:-samtools}"
PYTHON="${PYTHON:-python3}"
THREADS="${THREADS:-8}"
READ_LENGTH="${READ_LENGTH:-96}"
CHIM_SEGMENT_MIN="${CHIM_SEGMENT_MIN:-12}"
CHIM_JUNCTION_OVERHANG_MIN="${CHIM_JUNCTION_OVERHANG_MIN:-12}"
CHIM_MAIN_SEGMENT_MULT_NMAX="${CHIM_MAIN_SEGMENT_MULT_NMAX:-1}"
for path in "$INPUT_BAM" "$CALLED_BARCODES" "$GTF_FILE"; do
  [[ -s "$path" ]] || { echo "[ERROR] Missing input: $path" >&2; exit 2; }
done
for command_name in "$STAR" "$SAMTOOLS" "$PYTHON" gzip; do
  command -v "$command_name" >/dev/null 2>&1 || { echo "[ERROR] Command not found: $command_name" >&2; exit 2; }
done
for file in genomeParameters.txt Genome SA SAindex chrName.txt; do
  [[ -s "$STAR_INDEX/$file" ]] || { echo "[ERROR] Incomplete STAR index: $STAR_INDEX/$file" >&2; exit 2; }
done
for name in THREADS READ_LENGTH CHIM_SEGMENT_MIN CHIM_JUNCTION_OVERHANG_MIN CHIM_MAIN_SEGMENT_MULT_NMAX; do
  [[ "${!name}" =~ ^[1-9][0-9]*$ ]] || { echo "[ERROR] $name must be positive" >&2; exit 2; }
done
[[ -n "$OUT_DIR" && "$OUT_DIR" != / ]] || { echo "[ERROR] Unsafe OUT_DIR" >&2; exit 2; }

tmp="${OUT_DIR}.tmp.$$"
backup="${OUT_DIR}.previous.$$"
trap 'rm -rf "$tmp"' EXIT
[[ ! -e "$tmp" && ! -e "$backup" ]] || {
  echo "[ERROR] Temporary or backup path already exists" >&2
  exit 1
}
mkdir -p "$tmp/extracted" "$tmp/star" "$tmp/results" "$tmp/run_metadata"
"$PYTHON" "$code_dir/public_k562/extract_primary_transcripts.py" \
  --input "$INPUT_BAM" --samtools "$SAMTOOLS" --read-length "$READ_LENGTH" \
  --fastq "$tmp/extracted/SRR5082088.primary.fastq.gz" \
  --metadata "$tmp/extracted/SRR5082088.read_metadata.tsv"

"$STAR" \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$tmp/extracted/SRR5082088.primary.fastq.gz" \
  --readFilesCommand gzip -cd \
  --outFileNamePrefix "$tmp/star/SRR5082088_" \
  --outSAMtype BAM Unsorted \
  --runThreadN "$THREADS" \
  --chimSegmentMin "$CHIM_SEGMENT_MIN" \
  --chimJunctionOverhangMin "$CHIM_JUNCTION_OVERHANG_MIN" \
  --chimOutType Junctions WithinBAM SoftClip \
  --chimOutJunctionFormat 1 \
  --chimMainSegmentMultNmax "$CHIM_MAIN_SEGMENT_MULT_NMAX" \
  --outSAMunmapped Within \
  --outSAMattributes NH HI AS nM NM MD XS ch
[[ -e "$tmp/star/SRR5082088_Chimeric.out.junction" ]] || {
  echo "[ERROR] STAR did not create Chimeric.out.junction" >&2
  exit 1
}
"$SAMTOOLS" quickcheck "$tmp/star/SRR5082088_Aligned.out.bam"
"$PYTHON" "$code_dir/public_k562/summarize_chimeric_categories.py" \
  --junctions "$tmp/star/SRR5082088_Chimeric.out.junction" \
  --gtf "$GTF_FILE" \
  --read-metadata "$tmp/extracted/SRR5082088.read_metadata.tsv" \
  --called-barcodes "$CALLED_BARCODES" \
  --out-dir "$tmp/results"
cp "$config" "$tmp/run_metadata/config.sh"
{
  "$STAR" --version
  "$SAMTOOLS" --version | head -n 1
  "$PYTHON" --version
} > "$tmp/run_metadata/tool_versions.txt" 2>&1
sha256sum "$INPUT_BAM" "$CALLED_BARCODES" "$GTF_FILE" \
  "$code_dir/public_k562/extract_primary_transcripts.py" \
  "$code_dir/public_k562/summarize_chimeric_categories.py" > \
  "$tmp/run_metadata/input_checksums.sha256"
bash "$repo_root/scripts/write_source_provenance.sh" \
  "$repo_root" "$tmp/run_metadata/source_provenance.tsv"
[[ ! -e "$backup" ]] || { echo "[ERROR] Backup path already exists: $backup" >&2; exit 1; }
if [[ -e "$OUT_DIR" ]]; then mv "$OUT_DIR" "$backup"; fi
if ! mv "$tmp" "$OUT_DIR"; then
  [[ ! -e "$backup" ]] || mv "$backup" "$OUT_DIR"
  exit 1
fi
rm -rf "$backup"
trap - EXIT
echo "[DONE] Public K562 analysis: $OUT_DIR"

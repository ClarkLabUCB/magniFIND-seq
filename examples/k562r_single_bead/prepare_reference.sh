#!/usr/bin/env bash
set -euo pipefail

example_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$example_dir/../.." && pwd)"
reference_root="${REFERENCE_ROOT:-$repo_root/reference/GRCh38_GENCODEv44}"
sources_dir="$reference_root/sources"
star_index="$reference_root/star_index"
threads="${EXAMPLE_THREADS:-8}"

fasta_name="GRCh38.primary_assembly.genome.fa.gz"
gtf_name="gencode.v44.primary_assembly.annotation.gtf.gz"
fasta_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/$fasta_name"
gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/$gtf_name"
fasta_md5="9c3fc2ca260a767530dddb0f26721a6b"
gtf_md5="b182a9f3b134b9cc2da566a2b1692557"

for command_name in curl gzip md5sum STAR; do
  if ! command -v "$command_name" >/dev/null 2>&1; then
    echo "[ERROR] Required command not found: $command_name" >&2
    exit 2
  fi
done
if [[ -z "$reference_root" || "$reference_root" == "/" ]]; then
  echo "[ERROR] Refusing unsafe reference directory: $reference_root" >&2
  exit 1
fi
if [[ ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
  echo "[ERROR] EXAMPLE_THREADS must be a positive integer" >&2
  exit 1
fi

mkdir -p "$sources_dir" "$reference_root"
temporary_paths=()
cleanup() {
  local path
  for path in "${temporary_paths[@]}"; do
    if [[ -e "$path" ]]; then
      rm -rf "$path"
    fi
  done
}
trap cleanup EXIT

download_and_verify() {
  local url="$1"
  local destination="$2"
  local expected_md5="$3"
  local partial="${destination}.part"
  if [[ ! -s "$destination" ]]; then
    echo "[REFERENCE] Download $(basename "$destination")"
    temporary_paths+=("$partial")
    curl --fail --location --retry 5 --continue-at - \
      --output "$partial" "$url"
    mv "$partial" "$destination"
  fi
  local actual_md5
  actual_md5="$(md5sum "$destination" | awk '{print $1}')"
  if [[ "$actual_md5" != "$expected_md5" ]]; then
    echo "[ERROR] MD5 mismatch for $destination" >&2
    echo "[ERROR] Expected $expected_md5; observed $actual_md5" >&2
    exit 1
  fi
  gzip -t "$destination"
}

download_and_verify "$fasta_url" "$sources_dir/$fasta_name" "$fasta_md5"
download_and_verify "$gtf_url" "$sources_dir/$gtf_name" "$gtf_md5"

fasta="$reference_root/GRCh38.primary_assembly.genome.fa"
gtf="$reference_root/gencode.v44.primary_assembly.annotation.gtf"
if [[ ! -s "$fasta" ]]; then
  fasta_tmp="${fasta}.tmp.$$"
  temporary_paths+=("$fasta_tmp")
  echo "[REFERENCE] Decompress genome FASTA"
  gzip -cd "$sources_dir/$fasta_name" > "$fasta_tmp"
  mv "$fasta_tmp" "$fasta"
fi
if [[ ! -s "$gtf" ]]; then
  gtf_tmp="${gtf}.tmp.$$"
  temporary_paths+=("$gtf_tmp")
  echo "[REFERENCE] Decompress GTF"
  gzip -cd "$sources_dir/$gtf_name" > "$gtf_tmp"
  mv "$gtf_tmp" "$gtf"
fi

index_complete=1
for index_file in genomeParameters.txt Genome SA SAindex chrName.txt; do
  [[ -s "$star_index/$index_file" ]] || index_complete=0
done
if [[ "$index_complete" -eq 0 ]]; then
  if [[ -e "$star_index" ]]; then
    echo "[ERROR] STAR index exists but is incomplete: $star_index" >&2
    echo "[ERROR] Move or remove it before rebuilding." >&2
    exit 1
  fi
  index_tmp="$(mktemp -d "$reference_root/.star_index.XXXXXX")"
  temporary_paths+=("$index_tmp")
  echo "[REFERENCE] Build STAR index"
  STAR \
    --runMode genomeGenerate \
    --genomeDir "$index_tmp" \
    --genomeFastaFiles "$fasta" \
    --sjdbGTFfile "$gtf" \
    --sjdbOverhang 99 \
    --runThreadN "$threads"
  mv "$index_tmp" "$star_index"
else
  echo "[SKIP] Complete STAR index already exists"
fi

metadata_tmp="$reference_root/reference_metadata.tsv.tmp.$$"
temporary_paths+=("$metadata_tmp")
{
  printf 'field\tvalue\n'
  printf 'reference_name\tGRCh38_GENCODEv44\n'
  printf 'gencode_release\t44\n'
  printf 'genome_url\t%s\n' "$fasta_url"
  printf 'genome_md5\t%s\n' "$fasta_md5"
  printf 'gtf_url\t%s\n' "$gtf_url"
  printf 'gtf_md5\t%s\n' "$gtf_md5"
  printf 'STAR\t%s\n' "$(STAR --version 2>&1 | head -n 1)"
  printf 'sjdbOverhang\t99\n'
} > "$metadata_tmp"
mv "$metadata_tmp" "$reference_root/reference_metadata.tsv"

echo "[PASS] Reference prepared: $reference_root"

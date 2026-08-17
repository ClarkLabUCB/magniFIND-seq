#!/usr/bin/env bash
set -euo pipefail

example_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$example_dir/../.." && pwd)"
fixture_dir="$repo_root/tests/fastq_integration"
output_dir="$example_dir/output"
reference_dir="$output_dir/reference"
star_index="$reference_dir/star_index"
padded_reference="$reference_dir/mini_genome.padded.fa"

missing=0
for command_name in STAR trim_galore samtools featureCounts umi_tools Rscript python3 gzip sha256sum; do
  if ! command -v "$command_name" >/dev/null 2>&1; then
    echo "[ERROR] Required command not found: $command_name" >&2
    missing=1
  fi
done
if [[ "$missing" -ne 0 ]]; then
  echo "[ERROR] Create and activate the pinned root environment first." >&2
  exit 2
fi
Rscript -e 'suppressPackageStartupMessages({library(Matrix); library(Seurat)})' >/dev/null

if [[ "$(uname -s)" == "Darwin" ]]; then
  echo "[WARN] Linux is the supported platform for the complete STAR example." >&2
  echo "[WARN] Enabling the disk-backed FASTQ decompression fallback on macOS." >&2
  export PREDECOMPRESS_FASTQ="${PREDECOMPRESS_FASTQ:-1}"
fi

mkdir -p "$reference_dir"
if [[ ! -s "$padded_reference" ]]; then
  echo "[EXAMPLE] Prepare the deterministic miniature reference"
  python3 "$fixture_dir/prepare_test_reference.py" \
    --input "$fixture_dir/data/reference/mini_genome.fa" \
    --output "$padded_reference" \
    --target-length 100000
fi

index_complete=1
for index_file in genomeParameters.txt Genome SA SAindex chrName.txt; do
  if [[ ! -s "$star_index/$index_file" ]]; then
    index_complete=0
  fi
done

temporary_index=""
cleanup() {
  if [[ -n "$temporary_index" && -d "$temporary_index" ]]; then
    rm -rf "$temporary_index"
  fi
}
trap cleanup EXIT

if [[ "$index_complete" -eq 0 ]]; then
  if [[ -e "$star_index" ]]; then
    echo "[ERROR] The example STAR index exists but is incomplete: $star_index" >&2
    echo "[ERROR] Move or remove that directory before rerunning the example." >&2
    exit 1
  fi
  temporary_index="$(mktemp -d "$reference_dir/.star_index.XXXXXX")"
  echo "[EXAMPLE] Build the miniature STAR index"
  STAR \
    --runMode genomeGenerate \
    --genomeDir "$temporary_index" \
    --genomeFastaFiles "$padded_reference" \
    --sjdbGTFfile "$fixture_dir/data/reference/mini_genes.gtf" \
    --runThreadN "${EXAMPLE_THREADS:-1}" \
    --genomeSAindexNbases 7 \
    --genomeChrBinNbits 12 \
    --outFileNamePrefix "$reference_dir/index_build_"
  mv "$temporary_index" "$star_index"
  temporary_index=""
else
  echo "[SKIP] Complete miniature STAR index already exists"
fi

echo "[EXAMPLE] Run gene-expression quantification"
bash "$repo_root/01_preprocessing/code/run_gene_expression.sh" \
  "$example_dir/gene_expression.config.sh"

echo "[EXAMPLE] Run fusion mapping and BCR::ABL1 scan"
bash "$repo_root/02_fusion_analysis/code/run_targeted_mapping.sh" \
  "$example_dir/fusion.config.sh"
bash "$repo_root/02_fusion_analysis/code/run_targeted_scan.sh" \
  "$example_dir/fusion.config.sh"

echo "[EXAMPLE] Verify and summarize outputs"
python3 "$example_dir/verify_outputs.py"

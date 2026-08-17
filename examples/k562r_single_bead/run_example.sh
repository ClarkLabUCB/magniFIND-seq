#!/usr/bin/env bash
set -euo pipefail

example_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$example_dir/../.." && pwd)"
reference_root="${REFERENCE_ROOT:-$repo_root/reference/GRCh38_GENCODEv44}"

for index_file in genomeParameters.txt Genome SA SAindex chrName.txt; do
  if [[ ! -s "$reference_root/star_index/$index_file" ]]; then
    echo "[ERROR] Missing STAR reference index: $reference_root/star_index/$index_file" >&2
    echo "[ERROR] Run examples/k562r_single_bead/prepare_reference.sh first." >&2
    exit 2
  fi
done
if [[ ! -s "$reference_root/gencode.v44.primary_assembly.annotation.gtf" ]]; then
  echo "[ERROR] Missing GENCODE v44 GTF under $reference_root" >&2
  exit 2
fi

echo "[EXAMPLE] Run K562-r single-bead gene-expression workflow"
bash "$repo_root/01_preprocessing/code/run_gene_expression.sh" \
  "$example_dir/gene_expression.config.sh"

echo "[EXAMPLE] Run exploratory K562-r fusion mapping and scan"
bash "$repo_root/02_fusion_analysis/code/run_targeted_mapping.sh" \
  "$example_dir/fusion.config.sh"
bash "$repo_root/02_fusion_analysis/code/run_targeted_scan.sh" \
  "$example_dir/fusion.config.sh"

echo "[PASS] K562-r single-bead example completed"
echo "[INFO] Outputs: $example_dir/output"

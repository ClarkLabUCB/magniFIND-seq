#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fixture_dir="$repo_root/tests/fastq_integration"
data_dir="$fixture_dir/data"
expected_dir="$fixture_dir/expected"

require_command() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "[ERROR] Required command not found: $1" >&2
    missing=1
  fi
}

missing=0
for cmd in STAR trim_galore samtools featureCounts umi_tools Rscript python3 gzip diff; do
  require_command "$cmd"
done
if [[ "$missing" -ne 0 ]]; then
  echo "[ERROR] Create and activate the pinned root environment.yml first." >&2
  exit 2
fi
Rscript -e 'suppressPackageStartupMessages({library(Matrix); library(Seurat)})' >/dev/null

tmp_dir="$(mktemp -d)"
if [[ "${KEEP_TEST_OUTPUT:-0}" == "1" ]]; then
  trap 'echo "[INFO] Preserved test output: $tmp_dir"' EXIT
else
  trap 'rm -rf "$tmp_dir"' EXIT
fi

test_reference="$tmp_dir/mini_genome.padded.fa"
star_index="$tmp_dir/star_index"
gene_out="$tmp_dir/gene_expression"
fusion_out="$tmp_dir/fusion_detection"
mkdir -p "$star_index" "$gene_out" "$fusion_out"

predecompress_fastq=0
if [[ "$(uname -s)" == "Darwin" ]]; then
  # The current Bioconda osx-64 STAR build cannot spawn readFilesCommand under
  # Rosetta on Apple Silicon; exercise the documented disk-backed fallback.
  predecompress_fastq=1
fi

echo "[TEST] Prepare deterministic padded reference"
python3 "$fixture_dir/prepare_test_reference.py" \
  --input "$data_dir/reference/mini_genome.fa" \
  --output "$test_reference" \
  --target-length 100000

echo "[TEST] Build mini STAR index"
STAR \
  --runMode genomeGenerate \
  --genomeDir "$star_index" \
  --genomeFastaFiles "$test_reference" \
  --sjdbGTFfile "$data_dir/reference/mini_genes.gtf" \
  --runThreadN 1 \
  --genomeSAindexNbases 7 \
  --genomeChrBinNbits 12 \
  --outFileNamePrefix "$tmp_dir/index_build_"

sample="k562_r_fastq_fixture"
gene_sample_metadata="$tmp_dir/gene_sample_metadata.tsv"
cat > "$gene_sample_metadata" <<EOF
sample	condition
$sample	K562-R
EOF

gene_config="$tmp_dir/gene_expression.config.sh"
cat > "$gene_config" <<EOF
TEST_ROOT="$tmp_dir"
PROJECT_NAME="fastq_integration_gene_expression"
SAMPLES_TSV="$data_dir/samples.tsv"
SAMPLE_METADATA="$gene_sample_metadata"
OUT_DIR="\${TEST_ROOT}/gene_expression"
STAR="$(command -v STAR)"
STAR_INDEX="\${TEST_ROOT}/star_index"
GTF_FILE="$data_dir/reference/mini_genes.gtf"
THREADS=1
TSO_SEQUENCE="TTGCGCAATG"
UMI_LENGTH=8
TRIM_GALORE_OPTS="--paired --quality 1 --length 20"
FEATURE_TYPE="exon"
GENE_ATTRIBUTE="gene_id"
VALIDATE_FASTQ_CONTENTS=1
COUNT_READ_PAIRS=0
UMI_TOOLS_METHOD=""
UMI_TOOLS_PAIRED=0
STAR_OUT_SAM_ATTRIBUTES="NH HI AS nM XS"
PREDECOMPRESS_FASTQ=$predecompress_fastq
EOF

echo "[TEST] Run gene-expression quantification"
bash "$repo_root/01_preprocessing/code/run_gene_expression.sh" "$gene_config"

fusion_samples="$tmp_dir/fusion_samples.tsv"
cat > "$fusion_samples" <<EOF
sample	r1	r2
$sample	$gene_out/umi_trimmed/$sample/${sample}_R1_umi_out_val_1.fq.gz	$gene_out/umi_trimmed/$sample/${sample}_R2_umi_out_val_2.fq.gz
EOF

fusion_config="$tmp_dir/fusion.config.sh"
cat > "$fusion_config" <<EOF
TEST_ROOT="$tmp_dir"
PROJECT_NAME="fastq_integration_fusion"
SAMPLES_TSV="\${TEST_ROOT}/fusion_samples.tsv"
OUT_DIR="\${TEST_ROOT}/fusion_detection"
STAR="$(command -v STAR)"
STAR_INDEX="\${TEST_ROOT}/star_index"
THREADS=1
FUSION_TWOPASS_MODE="None"
CHIM_SEGMENT_MIN=12
CHIM_JUNCTION_OVERHANG_MIN=12
CHIM_MAIN_SEGMENT_MULT_NMAX=1
TAG_BAM_WITH_CB_UB=1
REQUIRE_CB_UB_TAGS=1
BCR_CHROM="chr22"
BCR_START=1
BCR_END=1000
BCR_STRAND="+"
ABL1_CHROM="chr9"
ABL1_START=1
ABL1_END=1000
ABL1_STRAND="+"
MIN_UNIQUE_UMI=1
REQUIRE_MAPPING_DONE=1
FAIL_ON_MISSING_JUNCTIONS=1
FAIL_ON_MALFORMED_RECORDS=1
PREDECOMPRESS_FASTQ=$predecompress_fastq
EOF

echo "[TEST] Run fusion mapping and scan"
bash "$repo_root/02_fusion_analysis/code/run_targeted_mapping.sh" "$fusion_config"
bash "$repo_root/02_fusion_analysis/code/run_targeted_scan.sh" "$fusion_config"

echo "[TEST] Compare exact gene counts, features, barcodes, and fusion summary"
python3 - "$expected_dir" "$gene_out" "$fusion_out" "$sample" <<'PY'
import csv
import gzip
import sys
from pathlib import Path

expected_dir, gene_out, fusion_out = map(Path, sys.argv[1:4])
sample = sys.argv[4]

with (expected_dir / "gene_counts.tsv").open() as handle:
    expected_counts = {
        row["gene"]: int(row["count"])
        for row in csv.DictReader(handle, delimiter="\t")
    }
with gzip.open(gene_out / "gene_by_sample_counts.tsv.gz", "rt") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    actual_counts = {row[reader.fieldnames[0]]: int(row[sample]) for row in reader}
if actual_counts != expected_counts:
    raise SystemExit(f"Gene counts differ: {actual_counts} != {expected_counts}")

expected_features = (expected_dir / "features.tsv").read_text()
actual_features = (gene_out / "features.tsv").read_text()
if actual_features != expected_features:
    raise SystemExit("features.tsv differs from exact expected output")
if (gene_out / "barcodes.tsv").read_text() != f"{sample}\n":
    raise SystemExit("barcodes.tsv differs from exact expected output")

summary = fusion_out / "fusion_scan" / "fastq_integration_fusion_BCR_ABL1_expected_orientation_summary_all_samples.tsv"
if summary.read_text() != (expected_dir / "fusion_summary.tsv").read_text():
    raise SystemExit("Fusion summary differs from exact expected output")

junctions = fusion_out / "fusion_scan" / "fastq_integration_fusion_BCR_ABL1_expected_orientation_junctions.tsv"
with junctions.open() as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
projection = [
    {
        "sample": row["sample"],
        "report_order": row["report_order"],
        "bcr_chrom": row["canonical_bcr_chrom"],
        "bcr_pos": row["canonical_bcr_pos"],
        "bcr_strand": row["canonical_bcr_strand"],
        "abl1_chrom": row["canonical_abl1_chrom"],
        "abl1_pos": row["canonical_abl1_pos"],
        "abl1_strand": row["canonical_abl1_strand"],
        "UB": row["UB"],
    }
    for row in rows
]
with (expected_dir / "fusion_junction_projection.tsv").open() as handle:
    expected_projection = list(csv.DictReader(handle, delimiter="\t"))
if projection != expected_projection:
    raise SystemExit(f"Fusion junction evidence differs: {projection} != {expected_projection}")

versions = fusion_out / "run_metadata" / "fusion_tool_versions.txt"
parameters = fusion_out / "fusion_scan" / "fastq_integration_fusion_fusion_scan_parameters.tsv"
if not versions.is_file() or "chim_segment_min\t12" not in versions.read_text():
    raise SystemExit("Fusion mapping provenance is missing or incomplete")
if not parameters.is_file() or "MIN_UNIQUE_UMI\t1" not in parameters.read_text():
    raise SystemExit("Fusion scan parameter provenance is missing or incomplete")
PY

echo "[TEST] Validate Seurat contents and BAM alignment preservation"
Rscript - "$gene_out/seurat/fastq_integration_gene_expression.seurat.rds" <<'RS'
suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
obj <- readRDS(args[[1]])
stopifnot(nrow(obj) == 2, ncol(obj) == 1)
stopifnot(identical(colnames(obj), "k562_r_fastq_fixture"))
stopifnot(obj$sample_id[[1]] == "k562_r_fastq_fixture")
stopifnot(obj$condition[[1]] == "K562-R")
stopifnot(length(obj@misc$pipeline_provenance) > 0)
stopifnot(length(obj@misc$session_info) > 0)
feature_meta <- obj[["RNA"]][[]]
stopifnot(all(c("gene_id", "gene_name", "gene_type") %in% colnames(feature_meta)))
stopifnot(setequal(feature_meta$gene_name, c("BCR_mock", "ABL1_mock")))
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
stopifnot(counts["BCR-mock", 1] == 2, counts["ABL1-mock", 1] == 1)
RS

fusion_bam="$fusion_out/fusion_mapped/$sample/${sample}_Aligned.sortedByCoord.out.bam"
tagged_bam="$fusion_out/fusion_mapped/$sample/${sample}_fusion_tagged.bam"
[[ "$(samtools view -c "$fusion_bam")" == "$(samtools view -c "$tagged_bam")" ]]
[[ "$(samtools view -c -f 4 "$fusion_bam")" == "$(samtools view -c -f 4 "$tagged_bam")" ]]
[[ "$(samtools view -c -f 2048 "$fusion_bam")" == "$(samtools view -c -f 2048 "$tagged_bam")" ]]
[[ "$(samtools view -c "$fusion_bam")" == "9" ]]
[[ "$(samtools view -c -f 2048 "$fusion_bam")" -gt 0 ]]
samtools view "$tagged_bam" | awk 'BEGIN { ok=0 } /CB:Z:k562_r_fastq_fixture/ && /UB:Z:/ { ok=1 } END { exit !ok }'

echo "[TEST] Verify completed stages are skipped on rerun"
gene_done="$gene_out/run_metadata/seurat_object.done"
fusion_done="$fusion_out/fusion_mapped/$sample/${sample}_fusion_mapping.done"
file_mtime_ns() {
  python3 -c 'import os,sys; print(os.stat(sys.argv[1]).st_mtime_ns)' "$1"
}
gene_stamp="$(file_mtime_ns "$gene_done")"
fusion_stamp="$(file_mtime_ns "$fusion_done")"
bash "$repo_root/01_preprocessing/code/run_gene_expression.sh" "$gene_config"
bash "$repo_root/02_fusion_analysis/code/run_targeted_mapping.sh" "$fusion_config"
[[ "$gene_stamp" == "$(file_mtime_ns "$gene_done")" ]]
[[ "$fusion_stamp" == "$(file_mtime_ns "$fusion_done")" ]]

echo "[PASS] FASTQ-to-Seurat and FASTQ-to-fusion integration outputs matched exactly"

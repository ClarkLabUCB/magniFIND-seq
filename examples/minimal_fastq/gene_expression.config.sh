#!/usr/bin/env bash

# Ready-to-run configuration for the fully synthetic minimal example.
example_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$example_dir/../.." && pwd)"

PROJECT_NAME="minimal_fastq_example"
SAMPLES_TSV="$repo_root/tests/fastq_integration/data/samples.tsv"
SAMPLE_METADATA="$example_dir/sample_metadata.tsv"
OUT_DIR="$example_dir/output/gene_expression"

STAR="STAR"
STAR_INDEX="$example_dir/output/reference/star_index"
GTF_FILE="$repo_root/tests/fastq_integration/data/reference/mini_genes.gtf"
THREADS="${EXAMPLE_THREADS:-1}"
READ_FILES_COMMAND="gzip -cd"
PREDECOMPRESS_FASTQ="${PREDECOMPRESS_FASTQ:-0}"

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

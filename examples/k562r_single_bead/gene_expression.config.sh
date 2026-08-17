#!/usr/bin/env bash

example_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$example_dir/../.." && pwd)"
reference_root="${REFERENCE_ROOT:-$repo_root/reference/GRCh38_GENCODEv44}"

PROJECT_NAME="k562r_single_bead_example"
SAMPLES_TSV="$example_dir/data/samples.tsv"
SAMPLE_METADATA="$example_dir/data/sample_metadata.tsv"
OUT_DIR="$example_dir/output/gene_expression"

STAR="STAR"
STAR_INDEX="$reference_root/star_index"
GTF_FILE="$reference_root/gencode.v44.primary_assembly.annotation.gtf"
THREADS="${EXAMPLE_THREADS:-8}"
READ_FILES_COMMAND="gzip -cd"
PREDECOMPRESS_FASTQ="${PREDECOMPRESS_FASTQ:-0}"

TSO_SEQUENCE="TTGCGCAATG"
UMI_LENGTH=8
TRIM_GALORE_OPTS="--illumina --paired --fastqc"
FEATURE_TYPE="exon"
GENE_ATTRIBUTE="gene_id"
VALIDATE_FASTQ_CONTENTS=1
COUNT_READ_PAIRS=0
UMI_TOOLS_METHOD=""
UMI_TOOLS_PAIRED=0
STAR_OUT_SAM_ATTRIBUTES="NH HI AS nM XS"

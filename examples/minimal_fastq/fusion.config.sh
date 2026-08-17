#!/usr/bin/env bash

# Ready-to-run configuration for the fully synthetic minimal example.
example_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_NAME="minimal_fastq_example"
SAMPLES_TSV="$example_dir/fusion_samples.tsv"
OUT_DIR="$example_dir/output/fusion_detection"

STAR="STAR"
STAR_INDEX="$example_dir/output/reference/star_index"
THREADS="${EXAMPLE_THREADS:-1}"
READ_FILES_COMMAND="gzip -cd"
PREDECOMPRESS_FASTQ="${PREDECOMPRESS_FASTQ:-0}"

FUSION_TWOPASS_MODE="None"
CHIM_SEGMENT_MIN=12
CHIM_JUNCTION_OVERHANG_MIN=12
CHIM_MAIN_SEGMENT_MULT_NMAX=1
TAG_BAM_WITH_CB_UB=1
REQUIRE_CB_UB_TAGS=1

# Coordinates are for the synthetic chr22/chr9 fixture, not GRCh38.
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

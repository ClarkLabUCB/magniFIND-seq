#!/usr/bin/env bash

# Project label used in output file names.
PROJECT_NAME="magnifind_seq"

# Tab-separated sample manifest with columns:
# sample, r1, r2
SAMPLES_TSV="targeted_samples.tsv"

# Output directory for this fusion detection module.
OUT_DIR="../output/fusion_detection"

# STAR executable and genome index.
STAR="STAR"
STAR_INDEX="/path/to/STAR_index"

# Threads used per sample.
THREADS=8

# FASTQ decompression command passed to STAR as separate words.
READ_FILES_COMMAND="gzip -cd"

# Set to 1 only when STAR cannot spawn the decompression command. This uses
# additional temporary disk space for uncompressed FASTQ files.
PREDECOMPRESS_FASTQ=0

# STAR chimeric alignment settings.
FUSION_TWOPASS_MODE="None"
CHIM_SEGMENT_MIN=12
CHIM_JUNCTION_OVERHANG_MIN=12
CHIM_MAIN_SEGMENT_MULT_NMAX=1

# Optional BAM tagging. The input read names are expected to contain
# |CB:<sample>|UB:<UMI>, as produced by the 01_preprocessing workflow.
TAG_BAM_WITH_CB_UB=1
REQUIRE_CB_UB_TAGS=1

# BCR and ABL1 genomic intervals for GRCh38. Adjust these if a different
# genome build, contig naming convention, or custom reference is used.
BCR_CHROM="chr22"
BCR_START=23179704
BCR_END=23318037
BCR_STRAND="+"

ABL1_CHROM="chr9"
ABL1_START=130713946
ABL1_END=130887675
ABL1_STRAND="+"

# Candidate calls are reported as:
# - not_detected: no expected-orientation BCR::ABL1 records
# - candidate_low_umi: >=1 record but unique UB below this cutoff
# - candidate_ge_min_umi: unique UB >= this cutoff
MIN_UNIQUE_UMI=2

# A completed mapping marker is required for production scans. Set this to 0
# only for a documented precomputed Chimeric.out.junction fixture.
REQUIRE_MAPPING_DONE=1

# Missing/incomplete mappings and malformed STAR junction records are written
# as not_evaluated and also make the scanner exit nonzero.
FAIL_ON_MISSING_JUNCTIONS=1
FAIL_ON_MALFORMED_RECORDS=1

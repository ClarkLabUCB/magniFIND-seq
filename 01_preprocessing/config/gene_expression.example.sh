#!/usr/bin/env bash

# Project label used for the Seurat object and output filenames.
PROJECT_NAME="magnifind_seq"

# Tab-separated manifest with columns: sample, r1, r2. Relative paths in this
# config are resolved from the config file directory; FASTQ paths inside the
# manifest are resolved from the manifest directory.
SAMPLES_TSV="samples.tsv"

# Optional tab-separated metadata with a required sample column. Additional
# columns (for example condition or batch) are copied into Seurat metadata.
SAMPLE_METADATA=""

# Output directory.
OUT_DIR="../output/gene_expression"

# Reference files and executables.
STAR="STAR"
STAR_INDEX="/path/to/star/index"
GTF_FILE="/path/to/genes.gtf"

# Runtime settings.
THREADS=8

# FASTQ decompression command passed to STAR as separate words.
READ_FILES_COMMAND="gzip -cd"

# Set to 1 only when STAR cannot spawn the decompression command. This uses
# additional temporary disk space for uncompressed FASTQ files.
PREDECOMPRESS_FASTQ=0

# UMI extraction settings.
TSO_SEQUENCE="TTGCGCAATG"
UMI_LENGTH=8

# Production options: explicit Illumina adapter, paired-end trimming, and FastQC.
TRIM_GALORE_OPTS="--illumina --paired --fastqc"

# Gene assignment settings. The aggregation code currently supports gene_id
# only; using another GTF attribute is rejected to prevent mislabeled matrices.
FEATURE_TYPE="exon"
GENE_ATTRIBUTE="gene_id"

# Validate every paired FASTQ record before processing. Set to 0 only after a
# separate full validation has already been completed.
VALIDATE_FASTQ_CONTENTS=1

# Production featureCounts used -p -B without --countReadPairs.
COUNT_READ_PAIRS=0

# The production umi_tools command relied on its default method and did not use
# --paired. Leave the method empty to avoid adding --method.
UMI_TOOLS_METHOD=""
UMI_TOOLS_PAIRED=0

# STAR attributes used by the production gene-expression mapping command.
STAR_OUT_SAM_ATTRIBUTES="NH HI AS nM XS"

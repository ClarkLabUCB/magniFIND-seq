#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE="${PYTHONDONTWRITEBYTECODE:-1}"

if [[ $# -ne 1 ]]; then
  echo "Usage: bash 01_preprocessing/code/run_gene_expression.sh <config.sh>" >&2
  exit 1
fi

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$CODE_DIR/../.." && pwd)"
CONFIG="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
CONFIG_DIR="$(dirname "$CONFIG")"

# shellcheck disable=SC1090
source "$CONFIG"

: "${PROJECT_NAME:?Missing PROJECT_NAME}"
: "${SAMPLES_TSV:?Missing SAMPLES_TSV}"
: "${OUT_DIR:?Missing OUT_DIR}"
: "${STAR:?Missing STAR}"
: "${STAR_INDEX:?Missing STAR_INDEX}"
: "${GTF_FILE:?Missing GTF_FILE}"

THREADS="${THREADS:-8}"
TSO_SEQUENCE="${TSO_SEQUENCE:-TTGCGCAATG}"
UMI_LENGTH="${UMI_LENGTH:-8}"
TRIM_GALORE_OPTS="${TRIM_GALORE_OPTS:---illumina --paired --fastqc}"
FEATURE_TYPE="${FEATURE_TYPE:-exon}"
GENE_ATTRIBUTE="${GENE_ATTRIBUTE:-gene_id}"
COUNT_READ_PAIRS="${COUNT_READ_PAIRS:-0}"
UMI_TOOLS_METHOD="${UMI_TOOLS_METHOD:-}"
UMI_TOOLS_PAIRED="${UMI_TOOLS_PAIRED:-0}"
STAR_OUT_SAM_ATTRIBUTES="${STAR_OUT_SAM_ATTRIBUTES:-NH HI AS nM XS}"
VALIDATE_FASTQ_CONTENTS="${VALIDATE_FASTQ_CONTENTS:-1}"
SAMPLE_METADATA="${SAMPLE_METADATA:-}"
PYTHON="${PYTHON:-python3}"
TRIM_GALORE="${TRIM_GALORE:-trim_galore}"
SAMTOOLS="${SAMTOOLS:-samtools}"
FEATURECOUNTS="${FEATURECOUNTS:-featureCounts}"
UMI_TOOLS="${UMI_TOOLS:-umi_tools}"
RSCRIPT="${RSCRIPT:-Rscript}"
READ_FILES_COMMAND="${READ_FILES_COMMAND:-gzip -cd}"
PREDECOMPRESS_FASTQ="${PREDECOMPRESS_FASTQ:-0}"

resolve_config_path() {
  local value="$1"
  if [[ "$value" = /* ]]; then
    printf '%s\n' "$value"
  else
    printf '%s\n' "$CONFIG_DIR/$value"
  fi
}

resolve_executable() {
  local value="$1"
  if [[ "$value" == */* ]]; then
    value="$(resolve_config_path "$value")"
    [[ -x "$value" ]] || { echo "[ERROR] Executable not found: $value" >&2; exit 1; }
    printf '%s\n' "$value"
  else
    command -v "$value" || { echo "[ERROR] Command not found: $value" >&2; exit 1; }
  fi
}

SAMPLES_TSV="$(resolve_config_path "$SAMPLES_TSV")"
OUT_DIR="$(resolve_config_path "$OUT_DIR")"
STAR_INDEX="$(resolve_config_path "$STAR_INDEX")"
GTF_FILE="$(resolve_config_path "$GTF_FILE")"
if [[ -n "$SAMPLE_METADATA" ]]; then
  SAMPLE_METADATA="$(resolve_config_path "$SAMPLE_METADATA")"
fi
STAR="$(resolve_executable "$STAR")"
PYTHON="$(resolve_executable "$PYTHON")"
TRIM_GALORE="$(resolve_executable "$TRIM_GALORE")"
SAMTOOLS="$(resolve_executable "$SAMTOOLS")"
FEATURECOUNTS="$(resolve_executable "$FEATURECOUNTS")"
UMI_TOOLS="$(resolve_executable "$UMI_TOOLS")"
RSCRIPT="$(resolve_executable "$RSCRIPT")"

[[ -s "$GTF_FILE" ]] || { echo "[ERROR] GTF not found: $GTF_FILE" >&2; exit 1; }
if [[ -n "$SAMPLE_METADATA" && ! -s "$SAMPLE_METADATA" ]]; then
  echo "[ERROR] Sample metadata not found: $SAMPLE_METADATA" >&2
  exit 1
fi
if [[ -z "$OUT_DIR" || "$OUT_DIR" == "/" ]]; then
  echo "[ERROR] Refusing unsafe OUT_DIR: $OUT_DIR" >&2
  exit 1
fi
if [[ "$GENE_ATTRIBUTE" != "gene_id" ]]; then
  echo "[ERROR] GENE_ATTRIBUTE currently supports gene_id only" >&2
  exit 1
fi
if [[ "$COUNT_READ_PAIRS" != "0" && "$COUNT_READ_PAIRS" != "1" ]]; then
  echo "[ERROR] COUNT_READ_PAIRS must be 0 or 1" >&2
  exit 1
fi
if [[ "$UMI_TOOLS_PAIRED" != "0" && "$UMI_TOOLS_PAIRED" != "1" ]]; then
  echo "[ERROR] UMI_TOOLS_PAIRED must be 0 or 1" >&2
  exit 1
fi
if [[ "$VALIDATE_FASTQ_CONTENTS" != "0" && "$VALIDATE_FASTQ_CONTENTS" != "1" ]]; then
  echo "[ERROR] VALIDATE_FASTQ_CONTENTS must be 0 or 1" >&2
  exit 1
fi
if [[ "$PREDECOMPRESS_FASTQ" != "0" && "$PREDECOMPRESS_FASTQ" != "1" ]]; then
  echo "[ERROR] PREDECOMPRESS_FASTQ must be 0 or 1" >&2
  exit 1
fi
if [[ ! "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
  echo "[ERROR] THREADS must be a positive integer" >&2
  exit 1
fi
if [[ ! "$UMI_LENGTH" =~ ^[1-9][0-9]*$ ]]; then
  echo "[ERROR] UMI_LENGTH must be a positive integer" >&2
  exit 1
fi
if [[ -z "$TSO_SEQUENCE" ]]; then
  echo "[ERROR] TSO_SEQUENCE must not be empty" >&2
  exit 1
fi

content_identity() {
  sha256sum "$1" | awk '{print $1}'
}

star_index_identities=()
for index_file in genomeParameters.txt Genome SA SAindex chrName.txt; do
  index_path="$STAR_INDEX/$index_file"
  [[ -s "$index_path" ]] || { echo "[ERROR] STAR index file not found: $index_path" >&2; exit 1; }
  star_index_identities+=("$index_file=$(content_identity "$index_path")")
done

mkdir -p "$OUT_DIR"/{umi_extracted,umi_trimmed,mapped,seurat,run_metadata}

active_tmp=""
transient_file=""
cleanup() {
  if [[ -n "$active_tmp" && -d "$active_tmp" ]]; then
    rm -rf "$active_tmp"
  fi
  if [[ -n "$transient_file" ]]; then
    rm -f "$transient_file"
  fi
}
trap cleanup EXIT

pipeline_source_sha="$(content_identity "$CODE_DIR/run_gene_expression.sh")"
config_sha="$(content_identity "$CONFIG")"
make_signature() {
  "$PYTHON" -c 'import hashlib,sys; print(hashlib.sha256("\0".join(sys.argv[1:]).encode()).hexdigest())' \
    "PIPELINE_SOURCE=$pipeline_source_sha" "CONFIG=$config_sha" "$@"
}

stage_complete() {
  local done_file="$1"
  local signature="$2"
  shift 2
  [[ -s "$done_file" ]] || return 1
  [[ "$(<"$done_file")" == "$signature" ]] || return 1
  local path
  for path in "$@"; do
    [[ -s "$path" ]] || return 1
  done
}

write_done() {
  local done_file="$1"
  local signature="$2"
  local tmp_done="${done_file}.tmp.$$"
  printf '%s\n' "$signature" > "$tmp_done"
  mv "$tmp_done" "$done_file"
}

normalized_manifest="$OUT_DIR/run_metadata/samples.normalized.tsv"
normalized_tmp="${normalized_manifest}.tmp.$$"
manifest_args=(--manifest "$SAMPLES_TSV" --output "$normalized_tmp")
if [[ "$VALIDATE_FASTQ_CONTENTS" == "1" ]]; then
  manifest_args+=(--check-fastq)
fi
"$PYTHON" "$REPO_ROOT/scripts/validate_sample_manifest.py" "${manifest_args[@]}"
mv "$normalized_tmp" "$normalized_manifest"

cp "$CONFIG" "$OUT_DIR/run_metadata/config.sh"
bash "$REPO_ROOT/scripts/write_source_provenance.sh" \
  "$REPO_ROOT" "$OUT_DIR/run_metadata/source_provenance.tsv"

seurat_version="$("$RSCRIPT" -e 'cat(as.character(packageVersion("Seurat")))')"
matrix_version="$("$RSCRIPT" -e 'cat(as.character(packageVersion("Matrix")))')"

versions_file="$OUT_DIR/run_metadata/tool_versions.txt"
{
  printf 'project\t%s\n' "$PROJECT_NAME"
  printf 'config\t%s\n' "$CONFIG"
  printf 'manifest\t%s\n' "$SAMPLES_TSV"
  printf 'star_index\t%s\n' "$STAR_INDEX"
  printf 'gtf\t%s\n' "$GTF_FILE"
  printf 'STAR\t%s\n' "$("$STAR" --version 2>&1 | head -n 1)"
  printf 'samtools\t%s\n' "$("$SAMTOOLS" --version 2>&1 | head -n 1)"
  printf 'featureCounts\t%s\n' "$("$FEATURECOUNTS" -v 2>&1 | head -n 1)"
  printf 'umi_tools\t%s\n' "$("$UMI_TOOLS" --version 2>&1 | head -n 1)"
  printf 'trim_galore\t%s\n' "$("$TRIM_GALORE" --version 2>&1 | head -n 1)"
  printf 'python\t%s\n' "$("$PYTHON" --version 2>&1)"
  printf 'R\t%s\n' "$("$RSCRIPT" --version 2>&1 | head -n 1)"
  printf 'Seurat\t%s\n' "$seurat_version"
  printf 'Matrix\t%s\n' "$matrix_version"
  printf 'count_read_pairs\t%s\n' "$COUNT_READ_PAIRS"
  printf 'umi_tools_method\t%s\n' "${UMI_TOOLS_METHOD:-default}"
  printf 'umi_tools_paired\t%s\n' "$UMI_TOOLS_PAIRED"
  printf 'star_out_sam_attributes\t%s\n' "$STAR_OUT_SAM_ATTRIBUTES"
} > "${versions_file}.tmp.$$"
mv "${versions_file}.tmp.$$" "$versions_file"

effective_parameters="$OUT_DIR/run_metadata/effective_parameters.tsv"
{
  printf 'parameter\tvalue\n'
  printf 'project_name\t%s\n' "$PROJECT_NAME"
  printf 'threads\t%s\n' "$THREADS"
  printf 'tso_sequence\t%s\n' "$TSO_SEQUENCE"
  printf 'umi_length\t%s\n' "$UMI_LENGTH"
  printf 'trim_galore_opts\t%s\n' "$TRIM_GALORE_OPTS"
  printf 'feature_type\t%s\n' "$FEATURE_TYPE"
  printf 'gene_attribute\t%s\n' "$GENE_ATTRIBUTE"
  printf 'count_read_pairs\t%s\n' "$COUNT_READ_PAIRS"
  printf 'umi_tools_method\t%s\n' "${UMI_TOOLS_METHOD:-default}"
  printf 'umi_tools_paired\t%s\n' "$UMI_TOOLS_PAIRED"
  printf 'star_out_sam_attributes\t%s\n' "$STAR_OUT_SAM_ATTRIBUTES"
  printf 'validate_fastq_contents\t%s\n' "$VALIDATE_FASTQ_CONTENTS"
  printf 'read_files_command\t%s\n' "$READ_FILES_COMMAND"
  printf 'predecompress_fastq\t%s\n' "$PREDECOMPRESS_FASTQ"
  printf 'sample_metadata\t%s\n' "$SAMPLE_METADATA"
  printf 'gtf_sha256\t%s\n' "$(content_identity "$GTF_FILE")"
  printf 'star_index_sha256\t%s\n' "$(make_signature "${star_index_identities[@]}")"
} > "${effective_parameters}.tmp.$$"
mv "${effective_parameters}.tmp.$$" "$effective_parameters"

echo "[INFO] Project: $PROJECT_NAME"
echo "[INFO] Samples: $normalized_manifest"
echo "[INFO] Output:  $OUT_DIR"

# The trusted shell configuration may provide a command plus arguments (for
# example, "gzip -cd"). Preserve these as separate STAR arguments.
# shellcheck disable=SC2206
read_files_command=($READ_FILES_COMMAND)
[[ "${#read_files_command[@]}" -gt 0 ]] || {
  echo "[ERROR] READ_FILES_COMMAND must not be empty" >&2
  exit 1
}
read_files_command[0]="$(resolve_executable "${read_files_command[0]}")"
# shellcheck disable=SC2206
star_out_sam_attributes=($STAR_OUT_SAM_ATTRIBUTES)
[[ "${#star_out_sam_attributes[@]}" -gt 0 ]] || {
  echo "[ERROR] STAR_OUT_SAM_ATTRIBUTES must not be empty" >&2
  exit 1
}

while IFS=$'\t' read -r sample r1 r2; do
  [[ "$sample" == "sample" || -z "${sample:-}" ]] && continue

  echo
  echo "=============================="
  echo "[SAMPLE] $sample"
  echo "=============================="

  extract_dir="$OUT_DIR/umi_extracted/$sample"
  trim_dir="$OUT_DIR/umi_trimmed/$sample"
  map_dir="$OUT_DIR/mapped/$sample"
  mkdir -p "$extract_dir" "$trim_dir" "$map_dir"

  umi_r1="$extract_dir/${sample}_R1_umi_out.fastq.gz"
  umi_r2="$extract_dir/${sample}_R2_umi_out.fastq.gz"
  non_umi_r1="$extract_dir/${sample}_R1_out.fastq.gz"
  non_umi_r2="$extract_dir/${sample}_R2_out.fastq.gz"
  extract_done="$extract_dir/${sample}_umi_extraction.done"
  trimmed_r1="$trim_dir/${sample}_R1_umi_out_val_1.fq.gz"
  trimmed_r2="$trim_dir/${sample}_R2_umi_out_val_2.fq.gz"
  trimmed_non_umi_r1="$trim_dir/${sample}_R1_out_val_1.fq.gz"
  trimmed_non_umi_r2="$trim_dir/${sample}_R2_out_val_2.fq.gz"
  trim_done="$trim_dir/${sample}_trimming.done"
  bam="$map_dir/${sample}_Aligned.sortedByCoord.out.bam"
  map_done="$map_dir/${sample}_gene_mapping.done"
  tagged_bam="$map_dir/${sample}_tagged.bam"
  tag_done="$map_dir/${sample}_bam_tagging.done"
  assigned_bam="$map_dir/${sample}_assigned_sorted.bam"
  counts="$map_dir/${sample}_counts.tsv.gz"
  count_done="$map_dir/${sample}_umi_counting.done"

  extract_signature="$(make_signature \
    "$(content_identity "$r1")" "$(content_identity "$r2")" \
    "$(content_identity "$CODE_DIR/gene_expression/extract_umis.py")" \
    "$(content_identity "$CODE_DIR/gene_expression/umi_methods.py")" \
    "TSO=$TSO_SEQUENCE" "UMI_LENGTH=$UMI_LENGTH" "SAMPLE=$sample")"
  if stage_complete "$extract_done" "$extract_signature" \
      "$umi_r1" "$umi_r2" "$non_umi_r1" "$non_umi_r2"; then
    echo "[SKIP] UMI extraction is complete"
  else
    echo "[UMI extraction] $sample"
    rm -f "$extract_done"
    active_tmp="$(mktemp -d "$extract_dir/.umi_extract.XXXXXX")"
    "$PYTHON" "$CODE_DIR/gene_expression/extract_umis.py" \
      --fastq "$r1" \
      --fastq2 "$r2" \
      --output "$active_tmp" \
      --sample_name "$sample" \
      --tso "$TSO_SEQUENCE" \
      --umi-length "$UMI_LENGTH"
    test -s "$active_tmp/${sample}_R1_umi_out.fastq.gz"
    test -s "$active_tmp/${sample}_R2_umi_out.fastq.gz"
    test -s "$active_tmp/${sample}_R1_out.fastq.gz"
    test -s "$active_tmp/${sample}_R2_out.fastq.gz"
    for name in \
      "${sample}_R1_umi_out.fastq.gz" "${sample}_R2_umi_out.fastq.gz" \
      "${sample}_R1_out.fastq.gz" "${sample}_R2_out.fastq.gz"; do
      mv "$active_tmp/$name" "$extract_dir/$name"
    done
    rm -rf "$active_tmp"
    active_tmp=""
    write_done "$extract_done" "$extract_signature"
  fi

  trim_signature="$(make_signature \
    "$(content_identity "$umi_r1")" "$(content_identity "$umi_r2")" \
    "$(content_identity "$non_umi_r1")" "$(content_identity "$non_umi_r2")" \
    "TRIM_GALORE=$TRIM_GALORE" \
    "TRIM_GALORE_VERSION=$("$TRIM_GALORE" --version 2>&1 | head -n 1)" \
    "OPTS=$TRIM_GALORE_OPTS")"
  if stage_complete "$trim_done" "$trim_signature" \
      "$trimmed_r1" "$trimmed_r2" "$trimmed_non_umi_r1" "$trimmed_non_umi_r2"; then
    echo "[SKIP] Adapter/quality trimming is complete"
  else
    echo "[trim_galore] $sample"
    rm -f "$trim_done"
    active_tmp="$(mktemp -d "$trim_dir/.trim.XXXXXX")"
    # shellcheck disable=SC2206
    trim_options=($TRIM_GALORE_OPTS)
    "$TRIM_GALORE" "${trim_options[@]}" -o "$active_tmp" "$umi_r1" "$umi_r2"
    non_umi_has_records="$("$PYTHON" -c \
      'import gzip,sys; f=gzip.open(sys.argv[1],"rb"); print(int(bool(f.read(1)))); f.close()' \
      "$non_umi_r1")"
    if [[ "$non_umi_has_records" == "1" ]]; then
      "$TRIM_GALORE" "${trim_options[@]}" -o "$active_tmp" "$non_umi_r1" "$non_umi_r2"
    else
      echo "[trim_galore] non-UMI pair is empty; writing empty trimmed pair"
      printf '' | gzip -n -c > "$active_tmp/${sample}_R1_out_val_1.fq.gz"
      printf '' | gzip -n -c > "$active_tmp/${sample}_R2_out_val_2.fq.gz"
    fi
    test -s "$active_tmp/${sample}_R1_umi_out_val_1.fq.gz"
    test -s "$active_tmp/${sample}_R2_umi_out_val_2.fq.gz"
    test -s "$active_tmp/${sample}_R1_out_val_1.fq.gz"
    test -s "$active_tmp/${sample}_R2_out_val_2.fq.gz"
    for generated in "$active_tmp"/*; do
      [[ -e "$generated" ]] || continue
      mv "$generated" "$trim_dir/"
    done
    rm -rf "$active_tmp"
    active_tmp=""
    write_done "$trim_done" "$trim_signature"
  fi

  map_signature="$(make_signature \
    "$(content_identity "$trimmed_r1")" "$(content_identity "$trimmed_r2")" \
    "${star_index_identities[@]}" \
    "STAR=$("$STAR" --version 2>&1 | head -n 1)" "MAPQ=255" \
    "OUTSAMATTRIBUTES=${star_out_sam_attributes[*]}" \
    "READ_FILES_COMMAND=${read_files_command[*]}" \
    "PREDECOMPRESS_FASTQ=$PREDECOMPRESS_FASTQ")"
  if stage_complete "$map_done" "$map_signature" "$bam" "${bam}.bai" \
      "$map_dir/${sample}_Log.final.out"; then
    echo "[SKIP] Gene-expression STAR mapping is complete"
  else
    echo "[STAR] $sample"
    rm -f "$map_done"
    rm -rf "$map_dir/${sample}_STARtmp"
    active_tmp="$(mktemp -d "$map_dir/.star.XXXXXX")"
    tmp_prefix="$active_tmp/${sample}_"
    star_r1="$trimmed_r1"
    star_r2="$trimmed_r2"
    if [[ "$PREDECOMPRESS_FASTQ" == "1" ]]; then
      star_r1="$active_tmp/.star_input_R1.fastq"
      star_r2="$active_tmp/.star_input_R2.fastq"
      "${read_files_command[@]}" "$trimmed_r1" > "$star_r1"
      "${read_files_command[@]}" "$trimmed_r2" > "$star_r2"
    fi
    star_command=(
      "$STAR"
      --genomeDir "$STAR_INDEX"
      --readFilesIn "$star_r1" "$star_r2"
    )
    if [[ "$PREDECOMPRESS_FASTQ" == "0" ]]; then
      star_command+=(--readFilesCommand "${read_files_command[@]}")
    fi
    star_command+=(
      --outFileNamePrefix "$tmp_prefix"
      --outSAMtype BAM SortedByCoordinate
      --runThreadN "$THREADS"
      --outSAMattributes "${star_out_sam_attributes[@]}"
      --outSAMattrRGline "ID:${sample}" "SM:${sample}"
    )
    "${star_command[@]}"
    rm -f "$active_tmp/.star_input_R1.fastq" "$active_tmp/.star_input_R2.fastq"

    tmp_bam="${tmp_prefix}Aligned.sortedByCoord.out.bam"
    filtered_bam="$active_tmp/${sample}_filtered.bam"
    "$SAMTOOLS" view -@ "$THREADS" -q 255 -b "$tmp_bam" > "$filtered_bam"
    "$SAMTOOLS" quickcheck "$filtered_bam"
    mv "$filtered_bam" "$tmp_bam"
    "$SAMTOOLS" index -@ "$THREADS" "$tmp_bam"
    for generated in "$active_tmp/${sample}_"*; do
      [[ -e "$generated" ]] || continue
      mv "$generated" "$map_dir/"
    done
    rm -rf "$active_tmp"
    active_tmp=""
    test -s "$bam"
    test -s "${bam}.bai"
    write_done "$map_done" "$map_signature"
  fi

  tag_signature="$(make_signature \
    "$(content_identity "$bam")" "$(content_identity "$REPO_ROOT/scripts/tag_bam_with_cb_umi.py")")"
  if stage_complete "$tag_done" "$tag_signature" "$tagged_bam"; then
    echo "[SKIP] BAM CB/UB tagging is complete"
  else
    echo "[BAM tag CB/UB] $sample"
    rm -f "$tag_done"
    tagged_tmp="${tagged_bam}.tmp.$$"
    transient_file="$tagged_tmp"
    "$PYTHON" "$REPO_ROOT/scripts/tag_bam_with_cb_umi.py" \
      --input "$bam" \
      --output "$tagged_tmp" \
      --require-all-tags
    "$SAMTOOLS" quickcheck "$tagged_tmp"
    input_records="$("$SAMTOOLS" view -c "$bam")"
    tagged_records="$("$SAMTOOLS" view -c "$tagged_tmp")"
    [[ "$input_records" == "$tagged_records" ]] || {
      echo "[ERROR] BAM tagging changed alignment count for $sample" >&2
      rm -f "$tagged_tmp"
      exit 1
    }
    mv "$tagged_tmp" "$tagged_bam"
    transient_file=""
    write_done "$tag_done" "$tag_signature"
  fi

  count_signature="$(make_signature \
    "$(content_identity "$tagged_bam")" "$(content_identity "$GTF_FILE")" \
    "FEATURE_TYPE=$FEATURE_TYPE" "GENE_ATTRIBUTE=$GENE_ATTRIBUTE" \
    "COUNT_READ_PAIRS=$COUNT_READ_PAIRS" \
    "UMI_METHOD=${UMI_TOOLS_METHOD:-default}" \
    "UMI_PAIRED=$UMI_TOOLS_PAIRED" \
    "FEATURECOUNTS=$("$FEATURECOUNTS" -v 2>&1 | head -n 1)" \
    "UMI_TOOLS=$("$UMI_TOOLS" --version 2>&1 | head -n 1)")"
  if stage_complete "$count_done" "$count_signature" "$assigned_bam" "${assigned_bam}.bai" "$counts"; then
    echo "[SKIP] Gene assignment and UMI counting are complete"
  else
    echo "[featureCounts + umi_tools] $sample"
    rm -f "$count_done"
    active_tmp="$(mktemp -d "$map_dir/.count.XXXXXX")"
    namesort_bam="$active_tmp/${sample}_namesort.bam"
    featurecounts_bam="${namesort_bam}.featureCounts.bam"
    gene_assigned="$active_tmp/${sample}_gene_assigned.txt"
    assigned_tmp="$active_tmp/${sample}_assigned_sorted.bam"
    counts_plain="$active_tmp/${sample}_counts.tsv"
    counts_tmp="$active_tmp/${sample}_counts.tsv.gz"
    count_log_tmp="$active_tmp/${sample}_count.log"

    "$SAMTOOLS" sort -@ "$THREADS" -n -o "$namesort_bam" "$tagged_bam"
    pair_args=(-p -B)
    if [[ "$COUNT_READ_PAIRS" == "1" ]]; then
      pair_args+=(--countReadPairs)
    fi
    "$FEATURECOUNTS" \
      -a "$GTF_FILE" \
      -F GTF \
      -t "$FEATURE_TYPE" \
      -g "$GENE_ATTRIBUTE" \
      -o "$gene_assigned" \
      -R BAM \
      -T "$THREADS" \
      "${pair_args[@]}" \
      "$namesort_bam"

    "$SAMTOOLS" sort -@ "$THREADS" "$featurecounts_bam" -o "$assigned_tmp"
    "$SAMTOOLS" index -@ "$THREADS" "$assigned_tmp"
    umi_tools_args=(
      count
      --extract-umi-method=tag
      --per-cell
      --per-gene
      --gene-tag=XT
      --assigned-status-tag=XS
      --cell-tag=CB
      --umi-tag=UB
    )
    if [[ -n "$UMI_TOOLS_METHOD" ]]; then
      umi_tools_args+=(--method "$UMI_TOOLS_METHOD")
    fi
    if [[ "$UMI_TOOLS_PAIRED" == "1" ]]; then
      umi_tools_args+=(--paired)
    fi
    "$UMI_TOOLS" "${umi_tools_args[@]}" \
      -I "$assigned_tmp" \
      -S "$counts_plain" \
      --log2stderr 2> "$count_log_tmp"

    test -s "$counts_plain"
    gzip -n -c "$counts_plain" > "$counts_tmp"
    gzip -t "$counts_tmp"
    "$SAMTOOLS" quickcheck "$assigned_tmp"
    mv "$assigned_tmp" "$assigned_bam"
    mv "${assigned_tmp}.bai" "${assigned_bam}.bai"
    mv "$counts_tmp" "$counts"
    mv "$gene_assigned" "$map_dir/${sample}_gene_assigned.txt"
    mv "${gene_assigned}.summary" "$map_dir/${sample}_gene_assigned.txt.summary"
    mv "$count_log_tmp" "$map_dir/${sample}_count.log"
    rm -rf "$active_tmp"
    active_tmp=""
    write_done "$count_done" "$count_signature"
  fi
done < "$normalized_manifest"

count_identities=()
while IFS=$'\t' read -r sample _ _; do
  [[ "$sample" == "sample" || -z "${sample:-}" ]] && continue
  count_identities+=("$(content_identity "$OUT_DIR/mapped/$sample/${sample}_counts.tsv.gz")")
done < "$normalized_manifest"

aggregate_done="$OUT_DIR/run_metadata/aggregation.done"
aggregate_signature="$(make_signature \
  "${count_identities[@]}" "$(content_identity "$GTF_FILE")" \
  "$(content_identity "$CODE_DIR/gene_expression/aggregate_umi_counts.py")" \
  "PYTHON=$("$PYTHON" --version 2>&1)")"
aggregate_outputs=(
  "$OUT_DIR/matrix.mtx"
  "$OUT_DIR/features.tsv"
  "$OUT_DIR/barcodes.tsv"
  "$OUT_DIR/gene_by_sample_counts.tsv.gz"
)
if stage_complete "$aggregate_done" "$aggregate_signature" "${aggregate_outputs[@]}"; then
  echo "[SKIP] Count aggregation is complete"
else
  echo "[Aggregate counts]"
  rm -f "$aggregate_done"
  active_tmp="$(mktemp -d "$OUT_DIR/.aggregate.XXXXXX")"
  "$PYTHON" "$CODE_DIR/gene_expression/aggregate_umi_counts.py" \
    --counts-root "$OUT_DIR/mapped" \
    --manifest "$normalized_manifest" \
    --gtf "$GTF_FILE" \
    --out-dir "$active_tmp"
  for name in matrix.mtx features.tsv barcodes.tsv gene_by_sample_counts.tsv.gz; do
    test -s "$active_tmp/$name"
    mv "$active_tmp/$name" "$OUT_DIR/$name"
  done
  rm -rf "$active_tmp"
  active_tmp=""
  write_done "$aggregate_done" "$aggregate_signature"
fi

seurat_out="$OUT_DIR/seurat/${PROJECT_NAME}.seurat.rds"
seurat_done="$OUT_DIR/run_metadata/seurat_object.done"
seurat_signature="$(make_signature \
  "$(content_identity "$OUT_DIR/matrix.mtx")" \
  "$(content_identity "$OUT_DIR/features.tsv")" \
  "$(content_identity "$OUT_DIR/barcodes.tsv")" \
  "$(content_identity "$CODE_DIR/gene_expression/create_seurat_object.R")" \
  "PROJECT=$PROJECT_NAME" "R=$("$RSCRIPT" --version 2>&1 | head -n 1)" \
  "SEURAT=$seurat_version" "MATRIX=$matrix_version" \
  "SAMPLE_METADATA=${SAMPLE_METADATA:+$(content_identity "$SAMPLE_METADATA")}")"
if stage_complete "$seurat_done" "$seurat_signature" "$seurat_out"; then
  echo "[SKIP] Seurat object is complete"
else
  echo "[Create Seurat object]"
  rm -f "$seurat_done"
  seurat_tmp="${seurat_out}.tmp.$$"
  transient_file="$seurat_tmp"
  seurat_args=(
    --matrix "$OUT_DIR/matrix.mtx"
    --features "$OUT_DIR/features.tsv"
    --barcodes "$OUT_DIR/barcodes.tsv"
    --project "$PROJECT_NAME"
    --provenance "$versions_file"
    --out "$seurat_tmp"
  )
  if [[ -n "$SAMPLE_METADATA" ]]; then
    seurat_args+=(--sample-metadata "$SAMPLE_METADATA")
  fi
  "$RSCRIPT" "$CODE_DIR/gene_expression/create_seurat_object.R" "${seurat_args[@]}"
  "$RSCRIPT" -e 'x <- readRDS(commandArgs(TRUE)[1]); stopifnot(nrow(x) > 0, ncol(x) > 0)' "$seurat_tmp"
  mv "$seurat_tmp" "$seurat_out"
  transient_file=""
  write_done "$seurat_done" "$seurat_signature"
fi

echo "[DONE] $seurat_out"

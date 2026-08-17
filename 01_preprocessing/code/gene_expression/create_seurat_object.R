#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

# Seurat 5.0.x Assay5 construction drops dimensions for a one-sample matrix.
# The v3 assay representation is stable for both single- and multi-sample runs.
options(Seurat.object.assay.version = "v3")

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) {
    return(default)
  }
  if (idx == length(args)) {
    stop(paste("Missing value for", flag))
  }
  args[[idx + 1]]
}

matrix_path <- get_arg("--matrix")
features_path <- get_arg("--features")
barcodes_path <- get_arg("--barcodes")
project <- get_arg("--project", "magnifind_seq")
out_path <- get_arg("--out")
provenance_path <- get_arg("--provenance")
sample_metadata_path <- get_arg("--sample-metadata")
min_cells_raw <- get_arg("--min-cells", "0")
min_features_raw <- get_arg("--min-features", "0")
min_cells <- suppressWarnings(as.integer(min_cells_raw))
min_features <- suppressWarnings(as.integer(min_features_raw))

if (is.null(matrix_path) || is.null(features_path) || is.null(barcodes_path) || is.null(out_path)) {
  stop("Usage: create_seurat_object.R --matrix matrix.mtx --features features.tsv --barcodes barcodes.tsv --project NAME --out object.rds")
}
if (!grepl("^[0-9]+$", min_cells_raw) || !grepl("^[0-9]+$", min_features_raw) ||
    is.na(min_cells) || is.na(min_features)) {
  stop("--min-cells and --min-features must be non-negative integers")
}

counts <- as(readMM(matrix_path), "CsparseMatrix")
features <- read.delim(features_path, header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.delim(barcodes_path, header = FALSE, stringsAsFactors = FALSE)

if (nrow(features) != nrow(counts)) {
  stop("features.tsv row count does not match matrix rows")
}
if (nrow(barcodes) != ncol(counts)) {
  stop("barcodes.tsv row count does not match matrix columns")
}
if (anyNA(counts@x) || any(!is.finite(counts@x)) || any(counts@x < 0) ||
    any(abs(counts@x - round(counts@x)) > 1e-8)) {
  stop("Matrix Market input must contain finite, non-negative integer counts")
}

if (ncol(features) < 3) {
  stop("features.tsv must contain gene_id, gene_name, and gene_type columns")
}
feature_table <- features[, 1:3, drop = FALSE]
colnames(feature_table) <- c("gene_id", "gene_name", "gene_type")
if (anyNA(feature_table$gene_id) || any(!nzchar(feature_table$gene_id))) {
  stop("features.tsv contains missing or empty gene IDs")
}
if (anyNA(barcodes[[1]]) || any(!nzchar(barcodes[[1]])) || anyDuplicated(barcodes[[1]])) {
  stop("barcodes.tsv contains missing, empty, or duplicate sample IDs")
}
gene_ids <- make.unique(as.character(feature_table$gene_id))
# Seurat replaces underscores in feature names with hyphens. Apply that
# transformation before object construction so the annotation table uses the
# exact same, collision-safe row names while retaining the original gene_id in
# feature metadata.
seurat_feature_ids <- make.unique(gsub("_", "-", gene_ids, fixed = TRUE))
rownames(counts) <- seurat_feature_ids
colnames(counts) <- barcodes[[1]]
rownames(feature_table) <- seurat_feature_ids

obj <- CreateSeuratObject(
  counts = counts,
  project = project,
  min.cells = min_cells,
  min.features = min_features
)

# CreateSeuratObject can remove features through min.cells and cells through
# min.features. Reindex the original annotation after that filtering instead
# of assigning shorter row names to the unfiltered feature table.
feature_table <- feature_table[rownames(obj), , drop = FALSE]
if (nrow(feature_table) != nrow(obj) || anyNA(feature_table$gene_id)) {
  stop("Feature metadata could not be aligned to the filtered Seurat object")
}
obj[["RNA"]] <- AddMetaData(obj[["RNA"]], metadata = feature_table)
obj@misc$feature_table <- feature_table
obj$sample_id <- colnames(obj)

if (!is.null(sample_metadata_path)) {
  sample_metadata <- read.delim(
    sample_metadata_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!"sample" %in% colnames(sample_metadata)) {
    stop("Sample metadata must contain a 'sample' column")
  }
  if (anyNA(sample_metadata$sample) || any(!nzchar(sample_metadata$sample)) ||
      anyDuplicated(sample_metadata$sample)) {
    stop("Sample metadata contains missing, empty, or duplicate sample IDs")
  }
  missing_samples <- setdiff(colnames(obj), sample_metadata$sample)
  if (length(missing_samples)) {
    stop("Sample metadata is missing samples: ", paste(missing_samples, collapse = ", "))
  }
  metadata_columns <- setdiff(colnames(sample_metadata), "sample")
  conflicting <- intersect(metadata_columns, colnames(obj@meta.data))
  if (length(conflicting)) {
    stop("Sample metadata columns conflict with existing metadata: ",
         paste(conflicting, collapse = ", "))
  }
  sample_metadata <- sample_metadata[match(colnames(obj), sample_metadata$sample), , drop = FALSE]
  for (column in metadata_columns) {
    obj[[column]] <- sample_metadata[[column]]
  }
}
obj@misc$pipeline_provenance <- if (!is.null(provenance_path)) {
  readLines(provenance_path, warn = FALSE)
} else {
  character()
}
obj@misc$session_info <- capture.output(sessionInfo())

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, out_path)

message("[INFO] Seurat object saved to: ", out_path)
message("[INFO] Genes: ", nrow(obj))
message("[INFO] Samples/cells: ", ncol(obj))

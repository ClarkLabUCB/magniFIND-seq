#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: create_fixture_seurat.R <output.rds>")

gene_symbols <- c(
  "GENEA1", "GENEA2", "GENEB1", "GENEB2", "GENESTABLE1",
  "GENESTABLE2", "MT-FIXTURE", "RPL-FIXTURE"
)
gene_ids <- paste0("fixture_gene_", seq_along(gene_symbols))
counts <- matrix(
  c(
    50, 45, 55, 48, 5, 4, 6, 5,
    25, 24, 28, 26, 3, 2, 4, 3,
    4, 5, 3, 4, 50, 48, 55, 52,
    2, 3, 2, 2, 25, 24, 28, 26,
    30, 31, 29, 30, 30, 31, 29, 30,
    10, 11, 9, 10, 10, 11, 9, 10,
    80, 85, 78, 82, 8, 9, 7, 8,
    7, 8, 6, 7, 70, 72, 68, 71
  ),
  nrow = 8,
  byrow = TRUE,
  dimnames = list(
    gene_ids,
    c("R_cell_1", "R_cell_2", "R_cell_3", "R_cell_4",
      "S_cell_1", "S_cell_2", "S_cell_3", "S_cell_4")
  )
)
pbmc_counts <- matrix(
  c(
    8, 7,
    6, 5,
    7, 8,
    5, 6,
    25, 24,
    12, 13,
    4, 5,
    6, 5
  ),
  nrow = 8,
  byrow = TRUE,
  dimnames = list(gene_ids, c("PBMC_cell_1", "PBMC_cell_2"))
)
counts <- cbind(counts, pbmc_counts)
object <- CreateSeuratObject(counts = as(counts, "dgCMatrix"), project = "fixture")
feature_metadata <- data.frame(
  gene_id = gene_ids,
  gene_name = gene_symbols,
  gene_type = "fixture",
  row.names = gene_ids,
  stringsAsFactors = FALSE
)
object[["RNA"]] <- AddMetaData(object[["RNA"]], metadata = feature_metadata)
object$condition <- c(rep("Resistant", 4), rep("Sensitive", 4), rep("PBMC", 2))
object$analysis_keep <- TRUE
saveRDS(object, args[[1]], compress = "xz")

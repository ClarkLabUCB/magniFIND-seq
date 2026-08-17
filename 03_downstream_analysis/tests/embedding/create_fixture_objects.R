#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(Matrix); library(Seurat) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: create_fixture_objects.R <out-dir>")
dir.create(args[1], recursive = TRUE, showWarnings = FALSE)
set.seed(7)
genes <- paste0("G", seq_len(40))
make_object <- function(prefix, shift) {
  counts <- matrix(rpois(40 * 6, lambda = 4), nrow = 40,
                   dimnames = list(paste0("id", seq_len(40)), paste0(prefix, seq_len(6))))
  counts[seq_len(8), ] <- counts[seq_len(8), ] + shift
  object <- CreateSeuratObject(counts = as(counts, "dgCMatrix"), min.cells = 0)
  feature_meta <- data.frame(gene_name = genes, row.names = rownames(object))
  object[["RNA"]] <- AddMetaData(object[["RNA"]], metadata = feature_meta)
  object
}
saveRDS(make_object("P", 0), file.path(args[1], "pbmc.rds"))
saveRDS(make_object("K", 8), file.path(args[1], "k562r.rds"))

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(ggplot2)
})

# Seurat's uwot wrapper derives its worker count from the active future plan.
# A sequential plan makes nearest-neighbor/UMAP execution reproducible across
# hosts instead of relying on machine-specific worker availability.
future::plan(future::sequential)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  index <- match(name, args)
  if (is.na(index)) return(default)
  if (index == length(args)) stop("Missing value for ", name)
  args[[index + 1]]
}
manifest_path <- normalizePath(get_arg("--input-manifest"), mustWork = TRUE)
candidate_path <- normalizePath(get_arg("--fusion-candidates"), mustWork = TRUE)
out_dir <- normalizePath(get_arg("--out-dir"), mustWork = TRUE)
k562_r_group <- get_arg("--k562-r-group")
seed <- as.integer(get_arg("--seed", "20260730"))
pca_components <- as.integer(get_arg("--pca-components", "30"))
joint_dims <- as.integer(get_arg("--joint-dims", "20"))
k562_dims <- as.integer(get_arg("--k562-dims", "15"))
joint_variable_features <- as.integer(get_arg("--joint-variable-features", "3000"))
k562_variable_features <- as.integer(get_arg("--k562-variable-features", "2000"))
numeric_settings <- c(
  seed, pca_components, joint_dims, k562_dims,
  joint_variable_features, k562_variable_features
)
if (!nzchar(k562_r_group) || anyNA(numeric_settings) || pca_components < 2 ||
    joint_dims < 2 || k562_dims < 2 || joint_variable_features < 2 ||
    k562_variable_features < 2 || joint_dims > pca_components ||
    k562_dims > pca_components) {
  stop("Invalid group, seed, dimension, or variable-feature setting")
}

manifest <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("dataset", "path", "format", "group", "assay", "feature_column")
if (!identical(colnames(manifest), required) || nrow(manifest) < 2 ||
    anyNA(manifest[, required[1:5]]) || any(!nzchar(manifest$dataset)) ||
    any(!nzchar(manifest$path)) || any(!nzchar(manifest$group)) ||
    anyDuplicated(manifest$dataset) || any(!manifest$format %in% c("seurat_rds", "10x_dir"))) {
  stop("Input manifest must have unique, valid dataset/path/format/group/assay rows")
}
seurat_rows <- manifest$format == "seurat_rds"
if (any(!nzchar(manifest$assay[seurat_rows]) | manifest$assay[seurat_rows] == ".")) {
  stop("Each seurat_rds row must name an explicit assay")
}
manifest_dir <- dirname(manifest_path)
resolve_input <- function(value) {
  path <- path.expand(value)
  if (!grepl("^/", path)) path <- file.path(manifest_dir, path)
  normalizePath(path, mustWork = TRUE)
}

collapse_features <- function(counts, names) {
  valid <- !is.na(names) & nzchar(names)
  counts <- counts[valid, , drop = FALSE]
  names <- names[valid]
  if (!length(names)) stop("No valid feature names in an input dataset")
  unique_names <- unique(names)
  if (anyDuplicated(names)) {
    collapse <- sparseMatrix(
      i = match(names, unique_names), j = seq_along(names), x = 1,
      dims = c(length(unique_names), length(names))
    )
    counts <- collapse %*% counts
  }
  rownames(counts) <- unique_names
  counts
}

load_dataset <- function(row) {
  path <- resolve_input(row$path)
  if (row$format == "10x_dir") {
    counts <- Read10X(path, gene.column = 2, unique.features = FALSE)
    if (is.list(counts)) {
      if (!"Gene Expression" %in% names(counts)) stop("10x input lacks Gene Expression")
      counts <- counts[["Gene Expression"]]
    }
    feature_names <- rownames(counts)
  } else {
    source <- readRDS(path)
    if (!inherits(source, "Seurat")) stop("RDS is not a Seurat object: ", path)
    if (!row$assay %in% Assays(source)) {
      stop("Assay not found in ", row$dataset, ": ", row$assay)
    }
    counts <- tryCatch(
      LayerData(source, assay = row$assay, layer = "counts"),
      error = function(e) GetAssayData(source, assay = row$assay, slot = "counts")
    )
    feature_names <- rownames(counts)
    if (nzchar(row$feature_column) && row$feature_column != ".") {
      feature_meta <- source[[row$assay]][[]]
      if (!row$feature_column %in% colnames(feature_meta)) {
        stop("Feature column not found in ", row$dataset, ": ", row$feature_column)
      }
      feature_names <- as.character(feature_meta[rownames(counts), row$feature_column])
    }
  }
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  if (anyNA(values) || any(!is.finite(values)) || any(values < 0) ||
      any(abs(values - round(values)) > 1e-8)) {
    stop("Counts must be finite, non-negative integers in dataset: ", row$dataset)
  }
  counts <- collapse_features(counts, feature_names)
  object <- CreateSeuratObject(counts = counts, project = row$dataset, min.cells = 0)
  object$original_sample <- colnames(object)
  object$dataset <- row$dataset
  object$group <- row$group
  RenameCells(object, add.cell.id = row$dataset)
}

objects <- lapply(seq_len(nrow(manifest)), function(i) load_dataset(manifest[i, ]))
names(objects) <- manifest$dataset
joint <- merge(objects[[1]], y = objects[-1], merge.data = FALSE)
joint <- JoinLayers(joint, assay = "RNA")

candidates <- read.delim(candidate_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!"sample" %in% colnames(candidates) || anyDuplicated(candidates$sample)) {
  stop("Fusion candidate table must contain unique sample IDs")
}
if ("expected_BCR_ABL1_candidate" %in% colnames(candidates)) {
  candidate_value <- suppressWarnings(as.integer(candidates$expected_BCR_ABL1_candidate))
} else if ("expected_BCR_ABL1_candidate_call" %in% colnames(candidates)) {
  candidate_value <- as.integer(startsWith(candidates$expected_BCR_ABL1_candidate_call, "candidate_"))
} else {
  stop("Fusion candidate table lacks a candidate column")
}
if (anyNA(candidate_value) || any(!candidate_value %in% 0:1)) {
  stop("Fusion candidate values must be binary")
}
candidate_lookup <- setNames(candidate_value, candidates$sample)
k562_samples <- unique(joint$original_sample[joint$group == k562_r_group])
missing_k562_candidates <- setdiff(k562_samples, candidates$sample)
if (length(missing_k562_candidates)) {
  stop(
    "Fusion candidate table is missing K562-r samples: ",
    paste(missing_k562_candidates, collapse = ", ")
  )
}
joint$fusion_candidate <- unname(candidate_lookup[joint$original_sample])
joint$fusion_candidate[is.na(joint$fusion_candidate)] <- 0L

embed <- function(object, label, variable_features, umap_dims) {
  if (ncol(object) < 4 || nrow(object) < 3) stop(label, " needs at least 4 cells and 3 genes")
  set.seed(seed)
  object <- NormalizeData(object, verbose = FALSE)
  object <- FindVariableFeatures(
    object, selection.method = "vst",
    nfeatures = min(variable_features, nrow(object)), verbose = FALSE
  )
  object <- ScaleData(object, features = VariableFeatures(object), verbose = FALSE)
  npcs <- min(
    pca_components, ncol(object) - 1L, length(VariableFeatures(object)) - 1L
  )
  if (npcs < 2) stop("Insufficient PCA dimensions for ", label)
  dims <- seq_len(min(umap_dims, npcs))
  if (length(dims) < 2) stop("Insufficient UMAP dimensions for ", label)
  object <- RunPCA(object, features = VariableFeatures(object), npcs = npcs, verbose = FALSE)
  object <- RunUMAP(
    object, reduction = "pca", dims = dims,
    n.neighbors = min(30L, ncol(object) - 1L), seed.use = seed,
    umap.method = "uwot", metric = "cosine",
    verbose = FALSE
  )
  object
}

joint <- embed(
  joint, "joint embedding", joint_variable_features, joint_dims
)
k562_cells <- colnames(joint)[joint$group == k562_r_group]
if (length(k562_cells) < 4) stop("Fewer than four K562-r cells were available")
k562_source <- subset(joint, cells = k562_cells)
k562_source <- JoinLayers(k562_source, assay = "RNA")
k562_counts <- LayerData(k562_source, assay = "RNA", layer = "counts")
k562_metadata <- k562_source@meta.data[, c(
  "original_sample", "dataset", "group", "fusion_candidate"
), drop = FALSE]
k562 <- CreateSeuratObject(
  counts = k562_counts, project = "K562r_only", min.cells = 0,
  meta.data = k562_metadata
)
k562 <- embed(
  k562, "K562-r-only embedding", k562_variable_features, k562_dims
)

results <- file.path(out_dir, "results")
figures <- file.path(out_dir, "figures")
dir.create(results, recursive = TRUE, showWarnings = FALSE)
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
saveRDS(joint, file.path(results, "joint_pbmc_k562r.seurat.rds"), compress = "xz")
saveRDS(k562, file.path(results, "k562r_only.seurat.rds"), compress = "xz")

write_coordinates <- function(object, path) {
  coords <- Embeddings(object, "umap")
  output <- data.frame(
    cell = rownames(coords), original_sample = object$original_sample,
    dataset = object$dataset, group = object$group,
    fusion_candidate = object$fusion_candidate,
    UMAP_1 = coords[, 1], UMAP_2 = coords[, 2], stringsAsFactors = FALSE
  )
  write.table(output, path, sep = "\t", quote = FALSE, row.names = FALSE)
  output
}
joint_data <- write_coordinates(joint, file.path(results, "joint_umap_coordinates.tsv"))
k562_data <- write_coordinates(k562, file.path(results, "k562r_only_umap_coordinates.tsv"))

plot_group <- ggplot(joint_data, aes(UMAP_1, UMAP_2, color = group)) +
  geom_point(size = 0.7, alpha = 0.8) + theme_classic() + coord_equal() +
  labs(title = "PBMC and K562-r joint embedding (no batch correction)", color = NULL)
plot_candidate <- function(data, title) {
  data$fusion_candidate <- factor(data$fusion_candidate, levels = c(0, 1))
  ggplot(data, aes(UMAP_1, UMAP_2, color = fusion_candidate)) +
    geom_point(size = 0.8, alpha = 0.85) + coord_equal() + theme_classic() +
    scale_color_manual(values = c("0" = "#D1D5DB", "1" = "#DC2626")) +
    labs(title = title, color = "BCR::ABL1\ncandidate")
}
ggsave(file.path(figures, "joint_umap_by_group.pdf"), plot_group, width = 7, height = 6)
ggsave(file.path(figures, "joint_umap_fusion_candidates.pdf"),
       plot_candidate(joint_data, "Fusion candidates in joint embedding"), width = 7, height = 6)
ggsave(file.path(figures, "k562r_only_umap_fusion_candidates.pdf"),
       plot_candidate(k562_data, "Fusion candidates in K562-r"), width = 7, height = 6)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("[DONE] embeddings written to ", out_dir)

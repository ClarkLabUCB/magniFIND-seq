#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
})

parse_args <- function(x) {
  if (length(x) %% 2 != 0 || any(!startsWith(x[seq(1, length(x), 2)], "--"))) {
    stop("Arguments must be supplied as --name value pairs")
  }
  keys <- sub("^--", "", x[seq(1, length(x), 2)])
  values <- x[seq(2, length(x), 2)]
  if (anyDuplicated(keys)) stop("Duplicate command-line argument")
  setNames(as.list(values), keys)
}

arg <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c(
  "bulk-matrix", "bulk-sample-metadata", "seurat-rds", "out-dir",
  "project-name", "condition-a", "condition-b", "bulk-count-suffix",
  "bulk-gene-id-column", "bulk-gene-name-column",
  "bulk-gene-biotype-column", "sc-assay", "sc-counts-layer",
  "sc-condition-column", "fdr-threshold", "min-abs-log2fc",
  "top-n-per-direction", "signature-require-bulk-deg",
  "require-full-signature", "min-bulk-replicates-per-condition",
  "sc-min-total-counts", "sc-cpm-pseudocount",
  "sc-de-test", "sc-de-logfc-threshold", "sc-de-min-pct",
  "exclude-gene-regex", "module-score-seed", "module-score-ctrl",
  "module-score-nbin"
)
missing <- setdiff(required, names(arg))
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))

number_arg <- function(name, lower = -Inf, strict = FALSE) {
  value <- suppressWarnings(as.numeric(arg[[name]]))
  if (!is.finite(value) || if (strict) value <= lower else value < lower) {
    stop("Invalid numeric value for --", name, ": ", arg[[name]])
  }
  value
}

integer_arg <- function(name, lower = 0L) {
  value <- suppressWarnings(as.integer(arg[[name]]))
  if (is.na(value) || value < lower || as.character(value) != arg[[name]]) {
    stop("Invalid integer value for --", name, ": ", arg[[name]])
  }
  value
}

bulk_path <- normalizePath(arg[["bulk-matrix"]], mustWork = TRUE)
bulk_meta_path <- normalizePath(arg[["bulk-sample-metadata"]], mustWork = TRUE)
seurat_path <- normalizePath(arg[["seurat-rds"]], mustWork = TRUE)
out_dir <- normalizePath(arg[["out-dir"]], mustWork = TRUE)
project_name <- arg[["project-name"]]
condition_a <- arg[["condition-a"]]
condition_b <- arg[["condition-b"]]
if (!nzchar(condition_a) || !nzchar(condition_b) || identical(condition_a, condition_b)) {
  stop("Condition A and condition B must be distinct non-empty labels")
}
if (any(grepl("[\t\r\n]", c(condition_a, condition_b)))) {
  stop("Condition labels must not contain tab or newline characters")
}

fdr_threshold <- number_arg("fdr-threshold", 0, strict = TRUE)
if (fdr_threshold > 1) stop("--fdr-threshold must be at most 1")
min_abs_log2fc <- number_arg("min-abs-log2fc", 0)
top_n <- integer_arg("top-n-per-direction", 1L)
signature_require_bulk_deg <- integer_arg("signature-require-bulk-deg", 0L)
if (!signature_require_bulk_deg %in% 0:1) {
  stop("--signature-require-bulk-deg must be 0 or 1")
}
require_full_signature <- integer_arg("require-full-signature", 0L)
if (!require_full_signature %in% 0:1) {
  stop("--require-full-signature must be 0 or 1")
}
min_bulk_replicates <- integer_arg("min-bulk-replicates-per-condition", 1L)
sc_min_total <- integer_arg("sc-min-total-counts", 1L)
sc_pseudocount <- number_arg("sc-cpm-pseudocount", 0, strict = TRUE)
sc_de_test <- arg[["sc-de-test"]]
if (!nzchar(sc_de_test) || grepl("[[:space:]]", sc_de_test)) {
  stop("--sc-de-test must be a non-empty Seurat test name without whitespace")
}
sc_de_logfc_threshold <- number_arg("sc-de-logfc-threshold", 0)
sc_de_min_pct <- number_arg("sc-de-min-pct", 0)
if (sc_de_min_pct > 1) stop("--sc-de-min-pct must be at most 1")
module_score_seed <- integer_arg("module-score-seed", 0L)
module_score_ctrl <- integer_arg("module-score-ctrl", 1L)
module_score_nbin <- integer_arg("module-score-nbin", 1L)

results_dir <- file.path(out_dir, "results")
figures_dir <- file.path(out_dir, "figures")
metadata_dir <- file.path(out_dir, "run_metadata")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

safe_cor_stats <- function(data, scope) {
  n <- nrow(data)
  pearson <- spearman <- p_value <- slope <- intercept <- NA_real_
  if (n >= 3 && sd(data$bulk_log2FC) > 0 && sd(data$sc_log2FC) > 0) {
    test <- cor.test(data$bulk_log2FC, data$sc_log2FC, method = "pearson")
    pearson <- unname(test$estimate)
    p_value <- test$p.value
    spearman <- cor(data$bulk_log2FC, data$sc_log2FC, method = "spearman")
    fit <- lm(sc_log2FC ~ bulk_log2FC, data = data)
    intercept <- unname(coef(fit)[1])
    slope <- unname(coef(fit)[2])
  }
  concordance <- if (n) {
    mean(sign(data$bulk_log2FC) == sign(data$sc_log2FC))
  } else {
    NA_real_
  }
  data.frame(
    scope = scope,
    n_genes = n,
    pearson_r = pearson,
    pearson_p_value = p_value,
    spearman_rho = spearman,
    direction_concordance = concordance,
    linear_slope = slope,
    linear_intercept = intercept,
    stringsAsFactors = FALSE
  )
}

# Bulk differential expression ------------------------------------------------
bulk <- read.delim(bulk_path, check.names = FALSE, stringsAsFactors = FALSE)
bulk_meta <- read.delim(
  bulk_meta_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!identical(colnames(bulk_meta), c("sample", "condition"))) {
  stop("Bulk sample metadata must have exactly two columns: sample and condition")
}
if (anyNA(bulk_meta) || any(!nzchar(bulk_meta$sample)) || anyDuplicated(bulk_meta$sample)) {
  stop("Bulk sample metadata contains missing, empty, or duplicate sample IDs")
}
if (!all(bulk_meta$condition %in% c(condition_a, condition_b))) {
  stop("Bulk metadata contains a condition other than condition A or B")
}
if (!all(c(condition_a, condition_b) %in% bulk_meta$condition)) {
  stop("Both configured conditions must be represented in bulk metadata")
}
bulk_replicate_counts <- table(factor(
  bulk_meta$condition, levels = c(condition_a, condition_b)
))
if (any(bulk_replicate_counts < min_bulk_replicates)) {
  stop(
    "Each bulk condition must contain at least ", min_bulk_replicates,
    " replicates; observed ", condition_a, "=", bulk_replicate_counts[[condition_a]],
    ", ", condition_b, "=", bulk_replicate_counts[[condition_b]]
  )
}

gene_columns <- unlist(arg[c(
  "bulk-gene-id-column", "bulk-gene-name-column", "bulk-gene-biotype-column"
)])
if (!all(gene_columns %in% colnames(bulk))) {
  stop("Bulk matrix is missing gene columns: ", paste(setdiff(gene_columns, colnames(bulk)), collapse = ", "))
}
count_columns <- paste0(bulk_meta$sample, arg[["bulk-count-suffix"]])
if (!all(count_columns %in% colnames(bulk))) {
  stop("Bulk matrix is missing count columns: ", paste(setdiff(count_columns, colnames(bulk)), collapse = ", "))
}

bulk_counts <- as.matrix(bulk[, count_columns, drop = FALSE])
storage.mode(bulk_counts) <- "double"
if (anyNA(bulk_counts) || any(!is.finite(bulk_counts)) || any(bulk_counts < 0) ||
    any(abs(bulk_counts - round(bulk_counts)) > 1e-8)) {
  stop("Bulk count columns must contain finite, non-negative integer counts")
}
if (any(colSums(bulk_counts) == 0)) stop("A bulk sample has zero total counts")
colnames(bulk_counts) <- bulk_meta$sample

bulk_genes <- data.frame(
  gene_id = bulk[[arg[["bulk-gene-id-column"]]]],
  gene_name = bulk[[arg[["bulk-gene-name-column"]]]],
  gene_biotype = bulk[[arg[["bulk-gene-biotype-column"]]]],
  stringsAsFactors = FALSE
)
if (anyNA(bulk_genes$gene_id) || any(!nzchar(bulk_genes$gene_id))) {
  stop("Bulk gene IDs cannot be missing or empty")
}
rownames(bulk_counts) <- make.unique(as.character(bulk_genes$gene_id))

bulk_group <- factor(bulk_meta$condition, levels = c(condition_b, condition_a))
y <- DGEList(counts = bulk_counts, group = bulk_group, genes = bulk_genes)
keep <- filterByExpr(y, group = bulk_group)
if (sum(keep) < 2) stop("Fewer than two genes passed edgeR filterByExpr")
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
design <- model.matrix(~ bulk_group)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = 2)
bulk_result <- topTags(qlf, n = Inf, sort.by = "PValue")$table

normalized_cpm <- cpm(y, normalized.lib.sizes = TRUE)
bulk_result$condition_a_mean_CPM <- rowMeans(
  normalized_cpm[rownames(bulk_result), bulk_group == condition_a, drop = FALSE]
)
bulk_result$condition_b_mean_CPM <- rowMeans(
  normalized_cpm[rownames(bulk_result), bulk_group == condition_b, drop = FALSE]
)
bulk_result$bulk_comparison <- paste0(condition_a, "_vs_", condition_b)
bulk_result <- bulk_result[, c(
  "gene_id", "gene_name", "gene_biotype", "logFC", "logCPM", "F", "PValue",
  "FDR", "condition_a_mean_CPM", "condition_b_mean_CPM", "bulk_comparison"
)]
write_tsv(bulk_result, file.path(results_dir, "bulk_edger_all_tested_genes.tsv"))

bulk_deg <- subset(bulk_result, FDR < fdr_threshold & abs(logFC) >= min_abs_log2fc)
write_tsv(bulk_deg, file.path(results_dir, "bulk_edger_deg.tsv"))

# The production K562 analysis consumed the archived differential-expression
# table produced by the upstream bulk workflow. When supplied, standardize that
# table and use it for signature ranking while retaining the raw-count edgeR run
# above as an independently reproducible diagnostic.
bulk_signature_result <- bulk_result
bulk_signature_source <- "recomputed_edger"
if ("bulk-de-table" %in% names(arg)) {
  de_path <- normalizePath(arg[["bulk-de-table"]], mustWork = TRUE)
  de_sep <- if (grepl("\\.csv(\\.gz)?$", de_path, ignore.case = TRUE)) "," else "\t"
  production_de <- read.table(
    de_path, header = TRUE, sep = de_sep, check.names = FALSE,
    stringsAsFactors = FALSE, quote = "\"", comment.char = ""
  )
  de_arg_names <- c(
    "bulk-de-gene-id-column", "bulk-de-gene-name-column",
    "bulk-de-gene-biotype-column", "bulk-de-logfc-column",
    "bulk-de-pvalue-column", "bulk-de-fdr-column",
    "bulk-de-logfc-multiplier"
  )
  missing_de_args <- setdiff(de_arg_names, names(arg))
  if (length(missing_de_args)) {
    stop("Missing production bulk-DE arguments: ", paste(missing_de_args, collapse = ", "))
  }
  de_columns <- unlist(arg[de_arg_names[1:6]])
  if (!all(de_columns %in% colnames(production_de))) {
    stop(
      "Production bulk-DE table is missing columns: ",
      paste(setdiff(de_columns, colnames(production_de)), collapse = ", ")
    )
  }
  logfc_multiplier <- suppressWarnings(as.numeric(arg[["bulk-de-logfc-multiplier"]]))
  if (!is.finite(logfc_multiplier) || logfc_multiplier == 0) {
    stop("--bulk-de-logfc-multiplier must be a finite non-zero number")
  }
  bulk_signature_result <- data.frame(
    gene_id = as.character(production_de[[arg[["bulk-de-gene-id-column"]]]]),
    gene_name = as.character(production_de[[arg[["bulk-de-gene-name-column"]]]]),
    gene_biotype = as.character(production_de[[arg[["bulk-de-gene-biotype-column"]]]]),
    logFC = suppressWarnings(as.numeric(production_de[[arg[["bulk-de-logfc-column"]]]])) *
      logfc_multiplier,
    PValue = suppressWarnings(as.numeric(production_de[[arg[["bulk-de-pvalue-column"]]]])),
    FDR = suppressWarnings(as.numeric(production_de[[arg[["bulk-de-fdr-column"]]]])),
    stringsAsFactors = FALSE
  )
  invalid_de <-
    is.na(bulk_signature_result$gene_id) | !nzchar(bulk_signature_result$gene_id) |
    is.na(bulk_signature_result$gene_name) | !nzchar(bulk_signature_result$gene_name) |
    !is.finite(bulk_signature_result$logFC) |
    !is.finite(bulk_signature_result$PValue) |
    !is.finite(bulk_signature_result$FDR) |
    bulk_signature_result$PValue < 0 | bulk_signature_result$PValue > 1 |
    bulk_signature_result$FDR < 0 | bulk_signature_result$FDR > 1
  if (any(invalid_de)) {
    stop("Production bulk-DE table contains invalid identifiers, fold changes, or probabilities")
  }
  bulk_signature_result$bulk_comparison <- paste0(condition_a, "_vs_", condition_b)
  bulk_signature_source <- "production_bulk_de_table"
  write_tsv(
    bulk_signature_result,
    file.path(results_dir, "production_bulk_de_standardized.tsv")
  )
}

write_tsv(
  transform(
    bulk_meta,
    library_size = colSums(bulk_counts),
    effective_library_size = y$samples$lib.size * y$samples$norm.factors
  ),
  file.path(metadata_dir, "bulk_sample_metadata.tsv")
)

pdf(file.path(figures_dir, "bulk_mds.pdf"), width = 7, height = 7)
plotMDS(
  y,
  labels = bulk_meta$sample,
  col = ifelse(bulk_group == condition_a, "#B91C1C", "#1D4ED8"),
  main = paste(condition_a, "vs", condition_b)
)
dev.off()

pdf(file.path(figures_dir, "bulk_bcv.pdf"), width = 7, height = 7)
plotBCV(y)
dev.off()

# Single-cell descriptive pseudobulk ------------------------------------------
seurat_obj <- readRDS(seurat_path)
sc_assay <- arg[["sc-assay"]]
sc_layer <- arg[["sc-counts-layer"]]
sc_condition_column <- arg[["sc-condition-column"]]
if (!sc_assay %in% Assays(seurat_obj)) stop("Seurat assay not found: ", sc_assay)
if (!sc_condition_column %in% colnames(seurat_obj@meta.data)) {
  stop("Seurat metadata column not found: ", sc_condition_column)
}

cell_meta <- seurat_obj@meta.data
# SC_KEEP_COLUMN defines the complete scoring population. Differential
# expression and pseudobulk summaries use only condition A/B within it, while
# additional groups (for example PBMCs) remain eligible for module scoring.
score_keep <- rep(TRUE, nrow(cell_meta))
if ("sc-keep-column" %in% names(arg)) {
  keep_column <- arg[["sc-keep-column"]]
  if (!keep_column %in% colnames(cell_meta)) stop("Seurat keep column not found: ", keep_column)
  keep_value <- cell_meta[[keep_column]]
  if (!is.logical(keep_value)) {
    keep_text <- tolower(trimws(as.character(keep_value)))
    allowed_keep <- c("true", "t", "1", "yes", "false", "f", "0", "no")
    invalid_keep <- !is.na(keep_text) & !keep_text %in% allowed_keep
    if (any(invalid_keep)) {
      stop("Seurat keep column contains values that are not logical: ",
           paste(unique(keep_text[invalid_keep]), collapse = ", "))
    }
    keep_value <- keep_text %in% c("true", "t", "1", "yes")
  }
  score_keep <- !is.na(keep_value) & keep_value
}
score_cells <- rownames(cell_meta)[score_keep]
if (!length(score_cells)) stop("No single cells passed the configured keep filter")
score_conditions <- as.character(cell_meta[score_cells, sc_condition_column])
if (anyNA(score_conditions) || any(!nzchar(score_conditions))) {
  stop("Retained scoring cells must have non-empty condition labels")
}
de_keep <- score_conditions %in% c(condition_a, condition_b)
cells <- score_cells[de_keep]
cell_conditions <- score_conditions[de_keep]
condition_counts <- table(factor(cell_conditions, levels = c(condition_a, condition_b)))
if (any(condition_counts == 0)) stop("Both configured conditions must have retained single cells")

sc_counts <- tryCatch(
  LayerData(seurat_obj, assay = sc_assay, layer = sc_layer),
  error = function(e) GetAssayData(seurat_obj, assay = sc_assay, slot = sc_layer)
)
sc_counts <- sc_counts[, score_cells, drop = FALSE]
count_values <- if (inherits(sc_counts, "sparseMatrix")) sc_counts@x else as.vector(sc_counts)
if (anyNA(count_values) || any(!is.finite(count_values)) || any(count_values < 0) ||
    any(abs(count_values - round(count_values)) > 1e-8)) {
  stop("Single-cell count layer must contain finite, non-negative integer counts")
}
if (any(colSums(sc_counts) == 0)) {
  stop("Each retained single cell must have at least one count")
}

sc_gene <- rownames(sc_counts)
if ("sc-feature-column" %in% names(arg)) {
  feature_column <- arg[["sc-feature-column"]]
  feature_meta <- seurat_obj[[sc_assay]][[]]
  if (!feature_column %in% colnames(feature_meta)) {
    stop("Seurat feature metadata column not found: ", feature_column)
  }
  sc_gene <- as.character(feature_meta[rownames(sc_counts), feature_column])
}
valid_gene <- !is.na(sc_gene) & nzchar(sc_gene)
sc_gene <- sc_gene[valid_gene]
sc_counts <- sc_counts[valid_gene, , drop = FALSE]
if (!length(sc_gene)) stop("No valid single-cell gene identifiers were available")
unique_gene <- unique(sc_gene)
if (anyDuplicated(sc_gene)) {
  collapse <- sparseMatrix(
    i = match(sc_gene, unique_gene),
    j = seq_along(sc_gene),
    x = 1,
    dims = c(length(unique_gene), length(sc_gene))
  )
  sc_counts <- collapse %*% sc_counts
}
rownames(sc_counts) <- unique_gene
sc_counts_de <- sc_counts[, cells, drop = FALSE]

sum_for_condition <- function(label) {
  Matrix::rowSums(sc_counts_de[, cell_conditions == label, drop = FALSE])
}
counts_a <- sum_for_condition(condition_a)
counts_b <- sum_for_condition(condition_b)
library_a <- sum(counts_a)
library_b <- sum(counts_b)
if (library_a == 0 || library_b == 0) stop("A single-cell pseudobulk library is empty")
cpm_a <- counts_a / library_a * 1e6
cpm_b <- counts_b / library_b * 1e6
detection_a <- Matrix::rowMeans(sc_counts_de[, cell_conditions == condition_a, drop = FALSE] > 0)
detection_b <- Matrix::rowMeans(sc_counts_de[, cell_conditions == condition_b, drop = FALSE] > 0)

sc_summary <- data.frame(
  gene = rownames(sc_counts),
  condition_a_counts = as.numeric(counts_a),
  condition_b_counts = as.numeric(counts_b),
  condition_a_CPM = as.numeric(cpm_a),
  condition_b_CPM = as.numeric(cpm_b),
  condition_a_detection_rate = as.numeric(detection_a),
  condition_b_detection_rate = as.numeric(detection_b),
  descriptive_sc_log2FC = log2(as.numeric(cpm_a) + sc_pseudocount) -
    log2(as.numeric(cpm_b) + sc_pseudocount),
  sc_comparison = paste0(condition_a, "_vs_", condition_b),
  stringsAsFactors = FALSE
)
write_tsv(sc_summary, file.path(results_dir, "single_cell_descriptive_pseudobulk.tsv"))
write_tsv(
  data.frame(
    condition = c(condition_a, condition_b),
    retained_cells = as.integer(condition_counts),
    pseudobulk_library_size = c(library_a, library_b)
  ),
  file.path(metadata_dir, "single_cell_group_summary.tsv")
)

# Single-cell differential expression ---------------------------------------
# Concordance in the manuscript is defined using Seurat::FindMarkers effect
# sizes, not the descriptive two-column pseudobulk calculation above. Build a
# dedicated object from the retained raw counts so the input object is never
# mutated and normalization is explicit.
score_obj <- CreateSeuratObject(
  counts = sc_counts,
  project = paste0(project_name, "_single_cell_de_and_score"),
  min.cells = 0,
  min.features = 0
)
score_obj$condition <- score_conditions
score_obj <- NormalizeData(
  score_obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)
Idents(score_obj) <- "condition"
sc_de <- FindMarkers(
  object = score_obj,
  ident.1 = condition_a,
  ident.2 = condition_b,
  assay = "RNA",
  slot = "data",
  test.use = sc_de_test,
  logfc.threshold = sc_de_logfc_threshold,
  min.pct = sc_de_min_pct,
  only.pos = FALSE,
  verbose = FALSE
)
fc_columns <- intersect(c("avg_log2FC", "avg_logFC"), colnames(sc_de))
if (length(fc_columns) != 1L) {
  stop("FindMarkers output did not contain exactly one recognized fold-change column")
}
sc_de$gene <- rownames(sc_de)
sc_de$sc_log2FC <- sc_de[[fc_columns]]
if (anyNA(sc_de$gene) || any(!nzchar(sc_de$gene)) || anyDuplicated(sc_de$gene)) {
  stop("FindMarkers returned invalid or duplicate gene identifiers")
}
if (any(!is.finite(sc_de$sc_log2FC))) {
  stop("FindMarkers returned a non-finite single-cell fold change")
}
sc_de <- sc_de[, c("gene", "sc_log2FC", setdiff(colnames(sc_de), c("gene", "sc_log2FC")))]
write_tsv(sc_de, file.path(results_dir, "single_cell_findmarkers.tsv"))

# Match genes and select the FDR-ranked signatures ----------------------------
bulk_match <- bulk_signature_result
bulk_match$gene <- bulk_match$gene_name
bulk_match$abs_bulk_log2FC <- abs(bulk_match$logFC)
bulk_match$bulk_log2FC <- bulk_match$logFC
bulk_match <- bulk_match[
  !is.na(bulk_match$gene) & nzchar(bulk_match$gene),
  ,
  drop = FALSE
]
bulk_match <- bulk_match[order(bulk_match$FDR, -bulk_match$abs_bulk_log2FC), ]
bulk_match <- bulk_match[!duplicated(bulk_match$gene), ]

matched <- merge(bulk_match, sc_summary, by = "gene", all = FALSE)
matched <- merge(matched, sc_de, by = "gene", all = FALSE)
matched$sc_total_counts <- matched$condition_a_counts + matched$condition_b_counts
matched$bulk_log2FC <- matched$logFC
matched$comparison <- paste0(condition_a, "_vs_", condition_b)
matched$bulk_DEG <- matched$FDR < fdr_threshold &
  abs(matched$bulk_log2FC) >= min_abs_log2fc
write_tsv(matched, file.path(results_dir, "bulk_single_cell_matched_genes.tsv"))

bulk_deg_candidate <- matched[
  matched$bulk_DEG & matched$sc_total_counts >= sc_min_total,
  , drop = FALSE
]
candidate <- matched[matched$sc_total_counts >= sc_min_total, , drop = FALSE]
if (signature_require_bulk_deg == 1L) candidate <- candidate[candidate$bulk_DEG, ]
exclude_regex <- arg[["exclude-gene-regex"]]
if (nzchar(exclude_regex)) candidate <- candidate[!grepl(exclude_regex, candidate$gene), ]

select_direction <- function(data, positive, label) {
  selected <- if (positive) data[data$bulk_log2FC > 0, ] else data[data$bulk_log2FC < 0, ]
  selected <- selected[order(selected$FDR, -abs(selected$bulk_log2FC)), ]
  selected <- head(selected, top_n)
  if (!nrow(selected)) stop("No eligible genes for ", label)
  if (nrow(selected) < top_n) {
    message <- paste0(
      "Only ", nrow(selected), " eligible genes were available for ", label,
      "; requested ", top_n
    )
    if (require_full_signature == 1L) stop(message) else warning(message)
  }
  selected$signature <- label
  selected$high_condition <- if (positive) condition_a else condition_b
  selected$rank_within_direction <- seq_len(nrow(selected))
  selected
}

selected_a <- select_direction(candidate, TRUE, "condition_a_high")
selected_b <- select_direction(candidate, FALSE, "condition_b_high")
selected <- rbind(selected_a, selected_b)
write_tsv(selected, file.path(results_dir, "bulk_top_genes.tsv"))

# Bulk-derived resistance score ----------------------------------------------
score_feature_keys <- setNames(rownames(score_obj), rownames(sc_counts))
score_features_a <- unname(score_feature_keys[selected_a$gene])
score_features_b <- unname(score_feature_keys[selected_b$gene])
if (anyNA(score_features_a) || anyNA(score_features_b)) {
  stop("A selected signature gene could not be represented in the scoring assay")
}
set.seed(module_score_seed)
score_obj <- AddModuleScore(
  object = score_obj,
  features = list(score_features_a, score_features_b),
  assay = "RNA",
  name = "BulkCore",
  ctrl = module_score_ctrl,
  nbin = module_score_nbin,
  seed = module_score_seed
)
score_data <- data.frame(
  sample = rownames(score_obj@meta.data),
  condition = as.character(score_obj$condition),
  condition_a_high_module_score = score_obj$BulkCore1,
  condition_b_high_module_score = score_obj$BulkCore2,
  resistance_score = score_obj$BulkCore1 - score_obj$BulkCore2,
  stringsAsFactors = FALSE
)
write_tsv(score_data, file.path(results_dir, "single_cell_resistance_scores.tsv"))

summarize_scores <- function(label) {
  values <- score_data$resistance_score[score_data$condition == label]
  data.frame(
    condition = label,
    n_samples = length(values),
    mean_resistance_score = mean(values),
    median_resistance_score = median(values),
    sd_resistance_score = if (length(values) > 1) sd(values) else NA_real_,
    q25_resistance_score = unname(quantile(values, 0.25)),
    q75_resistance_score = unname(quantile(values, 0.75)),
    stringsAsFactors = FALSE
  )
}
score_group_levels <- c(
  condition_b,
  condition_a,
  sort(setdiff(unique(score_data$condition), c(condition_a, condition_b)))
)
score_summary <- do.call(rbind, lapply(score_group_levels, summarize_scores))
write_tsv(score_summary, file.path(results_dir, "resistance_score_group_summary.tsv"))

score_data$condition <- factor(
  score_data$condition,
  levels = score_group_levels
)
extra_colors <- if (length(score_group_levels) > 2) {
  grDevices::hcl.colors(length(score_group_levels) - 2, palette = "Dark 3")
} else {
  character()
}
score_colors <- setNames(
  c("#1D4ED8", "#B91C1C", extra_colors),
  score_group_levels
)
p_score <- ggplot(
  score_data,
  aes(x = condition, y = resistance_score, fill = condition, color = condition)
) +
  geom_hline(yintercept = 0, color = "#6B7280", linewidth = 0.4) +
  geom_violin(trim = FALSE, width = 0.85, alpha = 0.3) +
  geom_boxplot(width = 0.18, outlier.shape = NA, color = "#111827", alpha = 0.75) +
  geom_jitter(width = 0.08, size = 1.3, alpha = 0.65) +
  scale_fill_manual(values = score_colors) +
  scale_color_manual(values = score_colors) +
  labs(
    title = "Bulk-derived resistance signature",
    subtitle = paste0(
      condition_a, "-high module score minus ", condition_b,
      "-high module score; ", nrow(selected_a), " genes per direction"
    ),
    x = NULL,
    y = "Resistance score"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", aspect.ratio = 1)
ggsave(
  file.path(figures_dir, "bulk_signature_resistance_score.pdf"),
  p_score,
  width = 6,
  height = 6,
  bg = "white"
)

metric_input <- function(x) {
  data.frame(bulk_log2FC = x$bulk_log2FC, sc_log2FC = x$sc_log2FC)
}
metrics <- rbind(
  safe_cor_stats(metric_input(matched[matched$sc_total_counts >= sc_min_total, ]), "all_matched_expressed_genes"),
  safe_cor_stats(metric_input(bulk_deg_candidate), "all_bulk_deg_candidates"),
  safe_cor_stats(metric_input(selected), "selected_genes_combined"),
  safe_cor_stats(metric_input(selected_a), "selected_condition_a_high"),
  safe_cor_stats(metric_input(selected_b), "selected_condition_b_high")
)
write_tsv(metrics, file.path(results_dir, "concordance_metrics.tsv"))

parameters <- data.frame(
  parameter = c(
    "project_name", "condition_a", "condition_b", "bulk_count_suffix",
    "sc_assay", "sc_counts_layer", "sc_condition_column", "sc_feature_column",
    "sc_keep_column", "fdr_threshold", "min_abs_log2fc",
    "top_n_per_direction", "signature_require_bulk_deg",
    "require_full_signature", "min_bulk_replicates_per_condition",
    "sc_min_total_counts", "sc_cpm_pseudocount",
    "sc_de_test", "sc_de_logfc_threshold", "sc_de_min_pct",
    "module_score_seed", "module_score_ctrl", "module_score_nbin",
    "bulk_signature_source",
    "exclude_gene_regex", "bulk_genes_before_filter", "bulk_genes_after_filter",
    "recomputed_bulk_deg_count", "signature_source_bulk_deg_count",
    "matched_gene_count", "selected_gene_count"
  ),
  value = c(
    project_name, condition_a, condition_b, arg[["bulk-count-suffix"]],
    sc_assay, sc_layer, sc_condition_column,
    ifelse("sc-feature-column" %in% names(arg), arg[["sc-feature-column"]], ""),
    ifelse("sc-keep-column" %in% names(arg), arg[["sc-keep-column"]], ""),
    fdr_threshold, min_abs_log2fc, top_n, signature_require_bulk_deg,
    require_full_signature, min_bulk_replicates,
    sc_min_total, sc_pseudocount,
    sc_de_test, sc_de_logfc_threshold, sc_de_min_pct,
    module_score_seed, module_score_ctrl, module_score_nbin,
    bulk_signature_source,
    exclude_regex, nrow(bulk), sum(keep), nrow(bulk_deg),
    sum(
      bulk_signature_result$FDR < fdr_threshold &
        abs(bulk_signature_result$logFC) >= min_abs_log2fc
    ),
    nrow(matched), nrow(selected)
  ),
  stringsAsFactors = FALSE
)
write_tsv(parameters, file.path(metadata_dir, "analysis_parameters.tsv"))

# Figures ---------------------------------------------------------------------
selected_lookup <- setNames(selected$signature, selected$gene)
plot_bulk <- bulk_match
plot_bulk$selected_signature <- unname(selected_lookup[plot_bulk$gene])
plot_bulk$display_group <- "Other genes"
plot_bulk$display_group[plot_bulk$FDR < fdr_threshold & plot_bulk$bulk_log2FC >= min_abs_log2fc] <-
  paste0("Other ", condition_a, "-high DEG")
plot_bulk$display_group[plot_bulk$FDR < fdr_threshold & plot_bulk$bulk_log2FC <= -min_abs_log2fc] <-
  paste0("Other ", condition_b, "-high DEG")
plot_bulk$display_group[plot_bulk$selected_signature == "condition_a_high"] <-
  paste0("Selected ", condition_a, "-high")
plot_bulk$display_group[plot_bulk$selected_signature == "condition_b_high"] <-
  paste0("Selected ", condition_b, "-high")
plot_bulk$minus_log10_FDR <- -log10(pmax(plot_bulk$FDR, .Machine$double.xmin))

selected_a_label <- paste0("Selected ", condition_a, "-high")
selected_b_label <- paste0("Selected ", condition_b, "-high")
other_a_label <- paste0("Other ", condition_a, "-high DEG")
other_b_label <- paste0("Other ", condition_b, "-high DEG")
plot_levels <- c("Other genes", other_b_label, other_a_label, selected_b_label, selected_a_label)
plot_colors <- setNames(c("#D1D5DB", "#93C5FD", "#FCA5A5", "#1D4ED8", "#B91C1C"), plot_levels)
plot_bulk$display_group <- factor(plot_bulk$display_group, levels = plot_levels)
volcano_labels <- rbind(head(selected_a, 8), head(selected_b, 8))

p_volcano <- ggplot(plot_bulk, aes(bulk_log2FC, minus_log10_FDR, color = display_group)) +
  geom_point(alpha = 0.8, size = 1.2) +
  geom_vline(xintercept = c(-min_abs_log2fc, min_abs_log2fc), linetype = "dashed", color = "#6B7280") +
  geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "#6B7280") +
  geom_text_repel(data = volcano_labels,
                  aes(x = bulk_log2FC, y = -log10(pmax(FDR, .Machine$double.xmin)), label = gene),
                  inherit.aes = FALSE, size = 3.2, max.overlaps = Inf) +
  scale_color_manual(values = plot_colors, drop = FALSE) +
  labs(
    title = paste("Bulk RNA-seq:", condition_a, "vs", condition_b),
    subtitle = paste("Top", top_n, "FDR-ranked genes in each direction are highlighted"),
    x = "Bulk log2 fold change", y = "-log10(FDR)", color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(aspect.ratio = 1, legend.position = "bottom")
ggsave(file.path(figures_dir, "bulk_volcano.pdf"), p_volcano, width = 8, height = 8, bg = "white")

display_a <- head(selected_a, min(15, nrow(selected_a)))
display_b <- head(selected_b, min(15, nrow(selected_b)))
display <- rbind(display_a, display_b)
dot_data <- rbind(
  data.frame(gene = display$gene, signature = display$signature, condition = condition_a,
             detection_rate = display$condition_a_detection_rate,
             log2_CPM_plus1 = log2(display$condition_a_CPM + 1)),
  data.frame(gene = display$gene, signature = display$signature, condition = condition_b,
             detection_rate = display$condition_b_detection_rate,
             log2_CPM_plus1 = log2(display$condition_b_CPM + 1))
)
dot_data$signature <- factor(
  dot_data$signature,
  levels = c("condition_a_high", "condition_b_high"),
  labels = c(paste("Bulk", condition_a, "high"), paste("Bulk", condition_b, "high"))
)
dot_data$condition <- factor(dot_data$condition, levels = c(condition_a, condition_b))
dot_data$gene <- factor(dot_data$gene, levels = rev(c(display_a$gene, display_b$gene)))
p_dot <- ggplot(dot_data, aes(condition, gene)) +
  geom_point(aes(size = detection_rate, color = log2_CPM_plus1)) +
  facet_wrap(~ signature, scales = "free_y", nrow = 1) +
  scale_size_continuous(name = "Detection rate", range = c(1, 9), limits = c(0, 1)) +
  scale_color_gradientn(name = "Pseudobulk log2(CPM + 1)",
                        colors = c("#E5E7EB", "#F4B4AA", "#EF7D73", "#B91C1C")) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.major.y = element_line(color = "#E5E7EB"))
ggsave(file.path(figures_dir, "bulk_signature_single_cell_dotplot.pdf"), p_dot,
       width = 9, height = 8, bg = "white")

plot_cor <- matched[matched$sc_total_counts >= sc_min_total, ]
plot_cor$selected_signature <- unname(selected_lookup[plot_cor$gene])
plot_cor$display_group <- "Other genes"
plot_cor$display_group[plot_cor$bulk_DEG & plot_cor$bulk_log2FC > 0] <- other_a_label
plot_cor$display_group[plot_cor$bulk_DEG & plot_cor$bulk_log2FC < 0] <- other_b_label
plot_cor$display_group[plot_cor$selected_signature == "condition_a_high"] <- selected_a_label
plot_cor$display_group[plot_cor$selected_signature == "condition_b_high"] <- selected_b_label
plot_cor$display_group <- factor(plot_cor$display_group, levels = plot_levels)
plot_cor <- plot_cor[order(as.numeric(plot_cor$display_group)), ]
selected_metric <- metrics[metrics$scope == "selected_genes_combined", ]
selected_labels <- rbind(head(selected_a, 7), head(selected_b, 7))

p_cor <- ggplot(plot_cor, aes(bulk_log2FC, sc_log2FC, color = display_group)) +
  geom_hline(yintercept = 0, color = "#9CA3AF", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "#9CA3AF", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#CBD5E1") +
  geom_point(alpha = 0.75, size = 1.2) +
  geom_abline(slope = selected_metric$linear_slope, intercept = selected_metric$linear_intercept,
              color = "#111827", linewidth = 1) +
  geom_text_repel(data = selected_labels, aes(x = bulk_log2FC, y = sc_log2FC, label = gene),
                  inherit.aes = FALSE, size = 3.2, max.overlaps = Inf) +
  scale_color_manual(values = plot_colors, drop = FALSE) +
  labs(
    title = "Bulk-single-cell differential-expression concordance",
    subtitle = sprintf(
      "Selected genes: Pearson r = %.2f; Spearman rho = %.2f; direction concordance = %.1f%%; n = %d",
      selected_metric$pearson_r, selected_metric$spearman_rho,
      100 * selected_metric$direction_concordance, selected_metric$n_genes
    ),
    caption = "Statistics are descriptive and conditioned on bulk FDR-based gene ranking.",
    x = paste("Bulk RNA-seq log2FC:", condition_a, "vs", condition_b),
    y = paste("Single-cell FindMarkers log2FC:", condition_a, "vs", condition_b),
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(aspect.ratio = 1, legend.position = "bottom")
ggsave(file.path(figures_dir, "bulk_single_cell_concordance.pdf"), p_cor,
       width = 8.5, height = 8.5, bg = "white")
ggsave(file.path(figures_dir, "bulk_single_cell_concordance.png"), p_cor,
       width = 8.5, height = 8.5, dpi = 300, bg = "white")

capture.output(sessionInfo(), file = file.path(metadata_dir, "R_sessionInfo.txt"))
message("[DONE] Selected ", nrow(selected_a), " ", condition_a, "-high and ",
        nrow(selected_b), " ", condition_b, "-high genes")

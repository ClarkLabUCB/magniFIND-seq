#!/usr/bin/env bash
set -euo pipefail

analysis_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

for cmd in Rscript; do
  command -v "$cmd" >/dev/null 2>&1 || {
    echo "[ERROR] Required command not found: $cmd" >&2
    exit 2
  }
done
Rscript -e 'stopifnot(requireNamespace("Matrix", quietly=TRUE), requireNamespace("Seurat", quietly=TRUE))' >/dev/null

cat > "$tmp/matrix.mtx" <<'EOF'
%%MatrixMarket matrix coordinate integer general
%
3 3 6
1 1 1
1 2 1
1 3 1
2 1 2
2 2 2
3 1 3
EOF
cat > "$tmp/features.tsv" <<'EOF'
gene_1	GENE1	test
gene2	GENE2	test
gene3	GENE3	test
EOF
printf 'cell1\ncell2\ncell3\n' > "$tmp/barcodes.tsv"

echo "[TEST] Align feature metadata after CreateSeuratObject filtering"
Rscript "$analysis_dir/code/gene_expression/create_seurat_object.R" \
  --matrix "$tmp/matrix.mtx" \
  --features "$tmp/features.tsv" \
  --barcodes "$tmp/barcodes.tsv" \
  --project fixture \
  --min-cells 2 \
  --out "$tmp/object.rds"

Rscript - "$tmp/object.rds" <<'RS'
suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
object <- readRDS(args[[1]])
stopifnot(identical(rownames(object), c("gene-1", "gene2")))
feature_meta <- object[["RNA"]][[]]
stopifnot(identical(as.character(feature_meta$gene_name), c("GENE1", "GENE2")))
stopifnot(identical(as.character(feature_meta$gene_id), c("gene_1", "gene2")))
stopifnot(identical(rownames(object@misc$feature_table), rownames(object)))
RS

echo "[PASS] Seurat creation fixture passed"

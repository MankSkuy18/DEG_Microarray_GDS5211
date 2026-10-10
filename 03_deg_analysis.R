# 02_preprocessing_normalization.R
# Preprocessing and normalization for GDS5211

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(limma)
  library(ggplot2)
})

# Load raw GDS object
if (!file.exists("results/data/gds5211_raw.rds")) {
  stop("Raw dataset not found. Please run 01_download_and_explore.R first.")
}

gds <- readRDS("results/data/gds5211_raw.rds")

# Extract expression matrix in a flexible way
extract_expr_matrix <- function(obj) {
  if (methods::is(obj, "ExpressionSet")) {
    return(exprs(obj))
  }
  if (methods::is(obj, "GDS")) {
    tab <- Table(obj)
    numeric_cols <- sapply(tab, is.numeric)
    if (sum(numeric_cols) >= 2) {
      expr <- as.matrix(tab[, numeric_cols, drop = FALSE])
      rownames(expr) <- tab[[1]]
      return(expr)
    }
    # fall back: keep all columns except meta columns if label names exist
    idx <- grep("^(ID_REF|ID|GeneSymbol|Symbol|Sample|X)", colnames(tab), invert = TRUE)
    if (length(idx) >= 2) {
      expr <- as.matrix(tab[, idx, drop = FALSE])
      rownames(expr) <- tab[[1]]
      return(expr)
    }
  }
  stop("Unable to extract expression matrix from GEO object.")
}

expr <- extract_expr_matrix(gds)

# Remove non-numeric values if present
expr <- apply(expr, 2, as.numeric)
expr <- as.matrix(expr)
rownames(expr) <- rownames(expr)

# Keep genes with enough variation
keep <- rowSums(!is.na(expr)) > 0
expr <- expr[keep, , drop = FALSE]

# Handle log transformation if needed
if (max(expr, na.rm = TRUE) <= 20) {
  expr <- log2(expr + 1)
  message("Applied log2(x + 1) transformation.")
}

# Remove genes with very low variance
if (nrow(expr) > 1) {
  gene_var <- apply(expr, 1, var, na.rm = TRUE)
  expr <- expr[gene_var > 0.05, , drop = FALSE]
}

# Save processed data
saveRDS(expr, file = "results/data/gds5211_processed.rds")

# Basic QC plot: density plot before and after processing
pdf("results/plots/density_plot_preprocessing.pdf", width = 9, height = 6)
  plotDensities(expr, legend = FALSE)
  title("Expression Density Across Samples")
dev.off()

# Boxplot summary
pdf("results/plots/boxplot_expression_distribution.pdf", width = 9, height = 6)
  boxplot(expr, outline = FALSE, las = 2, main = "Expression Distribution by Sample")
dev.off()

cat("Preprocessing completed.\n")
cat("Processed matrix dimensions:", dim(expr), "\n")
cat("Saved to results/data/gds5211_processed.rds\n")

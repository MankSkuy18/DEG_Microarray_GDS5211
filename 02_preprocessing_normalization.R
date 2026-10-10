# ============================================================================
# 02_PREPROCESSING_AND_NORMALIZATION.R - Process and normalize expression data
# ============================================================================

cat("\n========================================")
cat("\n=== STEP 2: Preprocessing & Normalization ===")
cat("\n========================================\n")

if (!file.exists("results/data/gds5211_raw.rds")) {
  stop("Raw dataset not found. Please run 01_download_and_explore.R first.")
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(limma)
  library(ggplot2)
  library(dplyr)
})

cat("Loading raw GDS object...\n")
gds <- readRDS("results/data/gds5211_raw.rds")
pheno <- pData(gds)

# ============================================================================
# Extract expression matrix
# ============================================================================

cat("\n--- Extracting Expression Matrix ---\n")

extract_expr_matrix <- function(obj) {
  if (methods::is(obj, "ExpressionSet")) {
    return(exprs(obj))
  }
  if (methods::is(obj, "GDS")) {
    tab <- Table(obj)
    rownames(tab) <- tab[[1]]
    # Find numeric columns (skip ID columns)
    numeric_cols <- sapply(tab[, -1], is.numeric)
    if (sum(numeric_cols) >= 2) {
      expr <- as.matrix(tab[, numeric_cols, drop = FALSE])
      return(expr)
    }
    stop("Could not identify numeric expression columns.")
  }
  stop("Unsupported object class.")
}

tryCatch({
  expr <- extract_expr_matrix(gds)
  cat("✓ Expression matrix extracted\n")
  cat("  Dimensions:", nrow(expr), "genes x", ncol(expr), "samples\n")
}, error = function(e) {
  stop("Failed to extract expression matrix:", e$message)
})

# ============================================================================
# Data cleaning and QC
# ============================================================================

cat("\n--- Quality Control & Filtering ---\n")

# Convert to numeric if needed
expr <- apply(expr, 2, as.numeric)
expr <- as.matrix(expr)

# Remove rows with all NAs or zeros
keep <- rowSums(!is.na(expr) & expr != 0) >= ncol(expr) * 0.5
original_genes <- nrow(expr)
expr <- expr[keep, , drop = FALSE]
cat("Genes removed (all NAs/zeros):", original_genes - nrow(expr), "\n")
cat("Genes remaining:", nrow(expr), "\n")

# Check if log transformation is needed
max_val <- max(expr, na.rm = TRUE)
min_val <- min(expr[expr > 0], na.rm = TRUE)
cat("\nExpression value range:", min_val, "to", max_val, "\n")

if (max_val <= 20) {
  cat("Applying log2(x + 1) transformation...\n")
  expr <- log2(expr + 1)
  cat("✓ Log transformation applied\n")
} else {
  cat("Data appears to be already log-transformed or normalized.\n")
}

# Remove low variance genes
cat("\nFiltering low-variance genes...\n")
if (nrow(expr) > 1) {
  gene_var <- apply(expr, 1, var, na.rm = TRUE)
  var_threshold <- quantile(gene_var, 0.25)
  keep_var <- gene_var >= var_threshold
  expr <- expr[keep_var, , drop = FALSE]
  cat("Genes removed (low variance):", sum(!keep_var), "\n")
  cat("Genes remaining:", nrow(expr), "\n")
}

# ============================================================================
# Quality Control Plots - Before Normalization
# ============================================================================

cat("\n--- Generating QC Plots ---\n")

# 1. Boxplot
cat("Creating boxplot...\n")
png(
  "results/plots/quality_control/01_boxplot_raw_expression.png",
  width = CONFIG$width * 120,
  height = CONFIG$height * 120,
  res = CONFIG$dpi
)
boxplot(expr,
        outline = FALSE,
        las = 2,
        main = "Expression Distribution by Sample (After Filtering)",
        ylab = "log2(Expression + 1)",
        xlab = "Sample")
dev.off()
cat("✓ Saved: 01_boxplot_raw_expression.png\n")

# 2. Density plot
cat("Creating density plot...\n")
p_density <- ggplot(data.frame(expr_long = as.vector(expr)),
                    aes(x = expr_long, fill = "Expression")) +
  geom_density(alpha = 0.6, color = "black") +
  labs(title = "Expression Density Distribution",
       x = "log2(Expression + 1)",
       y = "Density") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave("results/plots/quality_control/02_density_raw_expression.png",
       p_density,
       width = CONFIG$width,
       height = CONFIG$height,
       dpi = CONFIG$dpi)
cat("✓ Saved: 02_density_raw_expression.png\n")

# 3. Heatmap of sample correlation
cat("Creating sample correlation heatmap...\n")
corr_matrix <- cor(expr, use = "complete.obs")

png(
  "results/plots/quality_control/03_sample_correlation.png",
  width = CONFIG$width * 120,
  height = CONFIG$height * 120,
  res = CONFIG$dpi
)
heatmap(
  corr_matrix,
  main = "Sample Correlation Matrix",
  col = colorRampPalette(c("blue", "white", "red"))(50),
  symm = TRUE
)
dev.off()
cat("✓ Saved: 03_sample_correlation.png\n")

# ============================================================================
# Save processed data
# ============================================================================

saveRDS(expr, file = "results/data/gds5211_expression_matrix.rds")
saveRDS(pheno, file = "results/data/gds5211_phenotype.rds")
saveRDS(gds, file = "results/data/gds5211_gds_object.rds")

cat("\n✓ Processed data saved:\n")
cat("  - results/data/gds5211_expression_matrix.rds\n")
cat("  - results/data/gds5211_phenotype.rds\n")
cat("  - results/data/gds5211_gds_object.rds\n")

# Save summary statistics
summary_stats <- data.frame(
  Sample = colnames(expr),
  Mean = colMeans(expr, na.rm = TRUE),
  Median = apply(expr, 2, median, na.rm = TRUE),
  SD = apply(expr, 2, sd, na.rm = TRUE),
  Min = apply(expr, 2, min, na.rm = TRUE),
  Max = apply(expr, 2, max, na.rm = TRUE)
)

write.csv(summary_stats,
          file = "results/tables/01_sample_expression_summary.csv",
          row.names = FALSE)

cat("\n--- Expression Summary Statistics ---\n")
print(summary_stats)

cat("\n========================================")
cat("\n=== STEP 2 COMPLETE ===")
cat("\n========================================\n")

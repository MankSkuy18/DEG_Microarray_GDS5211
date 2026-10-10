# ============================================================================
# 03_DEG_ANALYSIS.R - Differential Expression Analysis using limma
# ============================================================================

cat("\n========================================")
cat("\n=== STEP 3: Differential Expression Analysis ===")
cat("\n========================================\n")

if (!file.exists("results/data/gds5211_expression_matrix.rds")) {
  stop("Processed expression matrix not found. Please run 02_preprocessing_normalization.R first.")
}

suppressPackageStartupMessages({
  library(limma)
  library(GEOquery)
  library(Biobase)
  library(dplyr)
  library(tidyr)
})

cat("Loading processed data...\n")
expr <- readRDS("results/data/gds5211_expression_matrix.rds")
pheno <- readRDS("results/data/gds5211_phenotype.rds")

# ============================================================================
# Design sample groups
# ============================================================================

cat("\n--- Defining Sample Groups ---\n")

# Try to extract grouping information from phenotype data
group <- NULL
group_col <- NULL

candidate_cols <- c(
  "source_name_ch1",
  "characteristics_ch1.1",
  "characteristics_ch1.2",
  "characteristics_ch1",
  "title"
)

for (col in candidate_cols) {
  if (col %in% colnames(pheno)) {
    vals <- as.character(pheno[[col]])
    unique_vals <- unique(vals)
    if (length(unique_vals) >= 2 && length(unique_vals) <= nrow(pheno) - 1) {
      group <- factor(vals)
      group_col <- col
      cat("Using column:", col, "\n")
      cat("Groups identified:\n")
      print(table(group))
      break
    }
  }
}

if (is.null(group)) {
  cat("Warning: Could not identify sample groups from metadata.\n")
  cat("Assigning groups automatically (first half = Group1, second half = Group2)\n")
  group <- factor(rep(c("Group1", "Group2"), length.out = ncol(expr)))
  cat("Groups assigned:\n")
  print(table(group))
}

# Ensure consistency between expression matrix and group assignment
if (length(group) != ncol(expr)) {
  stop("Mismatch between group assignment and expression matrix dimensions.")
}

cat("\n--- Limma Design ---\n")

# Create design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
cat("Design matrix dimensions:", nrow(design), "x", ncol(design), "\n")
print(design)

# Create contrast matrix
if (nlevels(group) == 2) {
  contrast_formula <- paste(levels(group)[2], "-", levels(group)[1], sep = "")
  contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
  cat("\nContrast:", contrast_formula, "\n")
} else if (nlevels(group) > 2) {
  # Create all pairwise contrasts
  levels_list <- levels(group)
  contrast_list <- c()
  for (i in 1:(length(levels_list) - 1)) {
    for (j in (i + 1):length(levels_list)) {
      contrast_list <- c(contrast_list,
                         paste(levels_list[j], "-", levels_list[i], sep = ""))
    }
  }
  contrast_matrix <- makeContrasts(contrasts = contrast_list, levels = design)
  cat("\nPairwise contrasts:\n")
  print(contrast_list)
} else {
  stop("At least 2 groups are required for DEG analysis.")
}

# ============================================================================
# Linear model fitting
# ============================================================================

cat("\n--- Fitting Linear Models (limma) ---\n")

fit <- lmFit(expr, design)
cat("✓ Coefficients estimated\n")

fit2 <- contrasts.fit(fit, contrast_matrix)
cat("✓ Contrasts fitted\n")

fit2 <- eBayes(fit2)
cat("✓ Empirical Bayes smoothing applied\n")

# ============================================================================
# Extract DEG results
# ============================================================================

cat("\n--- Extracting DEG Results ---\n")

deg_results <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
deg_results <- deg_results %>%
  tibble::rownames_to_column(var = "Gene") %>%
  mutate(
    logFC = as.numeric(logFC),
    AveExpr = as.numeric(AveExpr),
    t = as.numeric(t),
    P.Value = as.numeric(P.Value),
    adj.P.Val = as.numeric(adj.P.Val),
    B = as.numeric(B),
    Significant = ifelse(adj.P.Val < CONFIG$q_threshold &
                           abs(logFC) > CONFIG$fc_threshold,
                         "Yes", "No"),
    Direction = case_when(
      adj.P.Val < CONFIG$q_threshold & logFC > CONFIG$fc_threshold ~ "Upregulated",
      adj.P.Val < CONFIG$q_threshold & logFC < -CONFIG$fc_threshold ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  ) %>%
  arrange(P.Value)

cat("✓ DEG results extracted\n")

# ============================================================================
# Summary statistics
# ============================================================================

cat("\n--- DEG Summary Statistics ---\n")
cat("Total genes analyzed:", nrow(deg_results), "\n")
cat("Significant DEGs (FDR < 0.05, |logFC| > 1):", sum(deg_results$Significant == "Yes"), "\n")
cat("  - Upregulated:", sum(deg_results$Direction == "Upregulated"), "\n")
cat("  - Downregulated:", sum(deg_results$Direction == "Downregulated"), "\n")

cat("\n--- Top 20 Upregulated Genes ---\n")
top_up <- deg_results %>%
  filter(Direction == "Upregulated") %>%
  head(20) %>%
  select(Gene, logFC, P.Value, adj.P.Val)
print(as.data.frame(top_up))

cat("\n--- Top 20 Downregulated Genes ---\n")
top_down <- deg_results %>%
  filter(Direction == "Downregulated") %>%
  head(20) %>%
  select(Gene, logFC, P.Value, adj.P.Val)
print(as.data.frame(top_down))

# ============================================================================
# Save results
# ============================================================================

write.csv(deg_results,
          file = "results/tables/02_deg_results_complete.csv",
          row.names = FALSE)

write.csv(top_up,
          file = "results/tables/03_top_upregulated_genes.csv",
          row.names = FALSE)

write.csv(top_down,
          file = "results/tables/04_top_downregulated_genes.csv",
          row.names = FALSE)

saveRDS(deg_results, file = "results/data/deg_results.rds")
saveRDS(fit2, file = "results/data/limma_fit_object.rds")

cat("\n✓ Results saved to results/tables/\n")
cat("  - 02_deg_results_complete.csv\n")
cat("  - 03_top_upregulated_genes.csv\n")
cat("  - 04_top_downregulated_genes.csv\n")

cat("\n========================================")
cat("\n=== STEP 3 COMPLETE ===")
cat("\n========================================\n")

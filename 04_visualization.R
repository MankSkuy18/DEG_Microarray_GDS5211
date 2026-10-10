# 03_deg_analysis.R
# Differential expression analysis using limma

suppressPackageStartupMessages({
  library(limma)
  library(GEOquery)
  library(Biobase)
  library(dplyr)
  library(stringr)
})

if (!file.exists("results/data/gds5211_processed.rds")) {
  stop("Processed dataset not found. Run 02_preprocessing_normalization.R first.")
}

expr <- readRDS("results/data/gds5211_processed.rds")
gds <- readRDS("results/data/gds5211_raw.rds")

# Analyze sample metadata from GEO object
pheno <- pData(gds)

# Infer groups from metadata
candidate_cols <- c("source_name_ch1", "characteristics_ch1", "title", "description")
group <- NULL

for (col in candidate_cols) {
  if (col %in% colnames(pheno)) {
    vals <- pheno[[col]]
    if (length(unique(vals)) >= 2) {
      group <- vals
      break
    }
  }
}

if (is.null(group)) {
  stop("Unable to infer treatment/control groups from the GEO metadata. Please edit this script to define the sample groups manually.")
}

# Convert to a binary comparison when possible
# If the metadata contains keywords like control and treatment, this is inferred automatically.
if (length(unique(group)) == 2) {
  group <- factor(group)
  levels(group) <- c("Control", "Treatment")
} else {
  # In case more than two groups are present, use the first two unique groups only.
  group <- factor(group)
  group <- droplevels(group)
  if (length(levels(group)) < 2) {
    stop("At least two groups are required for differential expression analysis.")
  }
  if (length(levels(group)) > 2) {
    warning("More than two groups detected. The first two groups will be used for the contrast.")
    group <- factor(as.character(group), levels = levels(group)[1:2])
  }
}

# Ensure sample ordering matches expression matrix columns
if (length(group) != ncol(expr)) {
  message("Sample metadata length does not match expression matrix. Attempting to align using sample names.")
  sample_names <- colnames(expr)
  pheno_names <- rownames(pheno)
  if (length(intersect(sample_names, pheno_names)) > 0) {
    pheno_matched <- pheno[match(sample_names, pheno_names), , drop = FALSE]
    group <- pheno_matched[[candidate_cols[1]]]
  }
}

if (length(group) != ncol(expr)) {
  stop("Final alignment failed. Manual editing of group assignment may be required.")
}

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Build contrast matrix for comparison between first and second group
contrast <- paste(levels(group)[2], "-", levels(group)[1], sep = "")
contrast_matrix <- makeContrasts(contrasts = contrast, levels = design)

fit <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

deg_results <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")

deg_results <- deg_results %>%
  tibble::rownames_to_column(var = "Gene") %>%
  mutate(
    logFC = as.numeric(logFC),
    AveExpr = as.numeric(AveExpr),
    P.Value = as.numeric(P.Value),
    adj.P.Val = as.numeric(adj.P.Val),
    significance = case_when(
      adj.P.Val < 0.05 & logFC > 1 ~ "Upregulated",
      adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

write.csv(deg_results, file = "results/tables/deg_results_gds5211.csv", row.names = FALSE)

cat("DEG analysis complete.\n")
cat("Total significant DEGs:", sum(deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 1), "\n")
print(head(deg_results, 10))

saveRDS(deg_results, file = "results/data/deg_results_gds5211.rds")

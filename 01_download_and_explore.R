# ============================================================================
# 01_DOWNLOAD_AND_EXPLORE.R - Download GDS5211 and perform quality checks
# ============================================================================

cat("\n========================================")
cat("\n=== STEP 1: Download and Explore Data ===")
cat("\n========================================\n")

if (!file.exists("00_setup.R")) {
  source("00_setup.R")
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(dplyr)
})

cat("Attempting to download GDS5211 from NCBI GEO...\n")
cat("(This may take a few minutes on first download)\n\n")

tryCatch({
  # Download GEO dataset
  gds <- getGEO(CONFIG$gds_accession, 
                 destdir = "data/raw", 
                 GSElimits = NULL, 
                 getGPL = FALSE)
  
  if (length(gds) == 0) {
    stop("Dataset returned empty. Check GEO access.")
  }
  
  gds <- gds[[1]]
  cat("✓ Dataset downloaded successfully!\n\n")
  
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  stop("Failed to download dataset. Please check internet connection or GEO status.")
})

# ============================================================================
# Extract and display dataset information
# ============================================================================

cat("\n--- Dataset Information ---\n")
cat("Accession:", CONFIG$gds_accession, "\n")
cat("Class:", class(gds), "\n")

# Get phenotype/sample information
pheno <- pData(gds)
cat("Number of samples:", nrow(pheno), "\n")
cat("Sample metadata columns:", ncol(pheno), "\n\n")

cat("--- Sample Information ---\n")
print(head(pheno, 5))

# Extract expression data
if (methods::is(gds, "GDS")) {
  feature_table <- Table(gds)
  cat("\n--- Feature/Gene Information ---\n")
  cat("Number of genes/features:", nrow(feature_table), "\n")
  cat("Data columns:", ncol(feature_table), "\n")
  cat("Column names:\n")
  print(colnames(feature_table))
  
  cat("\nFirst 10 rows of data:\n")
  print(head(feature_table, 10))
}

# Identify sample groups from metadata
cat("\n--- Identifying Sample Groups ---\n")

group_col <- NULL
candidate_cols <- c(
  "source_name_ch1",
  "characteristics_ch1.1", 
  "characteristics_ch1.2",
  "characteristics_ch1",
  "title",
  "organism_ch1"
)

for (col in candidate_cols) {
  if (col %in% colnames(pheno)) {
    vals <- as.character(pheno[[col]])
    unique_vals <- length(unique(vals))
    if (unique_vals >= 2 && unique_vals <= nrow(pheno) - 1) {
      group_col <- col
      cat("Found grouping variable in column:", col, "\n")
      cat("Unique values:", unique_vals, "\n")
      cat("Groups:", paste(unique(vals), collapse = ", "), "\n")
      break
    }
  }
}

if (is.null(group_col)) {
  cat("Warning: Could not automatically identify sample groups.\n")
  cat("Manual assignment will be needed in subsequent steps.\n")
}

# Save raw GDS object
saveRDS(gds, file = "results/data/gds5211_raw.rds")
cat("\n✓ Raw GDS object saved to: results/data/gds5211_raw.rds\n")

cat("\n========================================")
cat("\n=== STEP 1 COMPLETE ===")
cat("\n========================================\n")

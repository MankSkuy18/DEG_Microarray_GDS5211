# 01_download_and_explore.R
# Download and explore the GDS5211 dataset from NCBI GEO

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

# Create folders if needed
if (!dir.exists("data")) dir.create("data", recursive = TRUE)
if (!dir.exists("results")) dir.create("results", recursive = TRUE)
if (!dir.exists("results/plots")) dir.create("results/plots", recursive = TRUE)
if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
if (!dir.exists("results/data")) dir.create("results/data", recursive = TRUE)

cat("Downloading GDS5211 from NCBI GEO...\n")

# Download dataset from GEO
# The dataset is a GEO GDS record, not necessarily a GSE matrix.
# We can fetch it by accession directly.
gds <- getGEO("GDS5211", destdir = "data", GSElimits = NULL, getGPL = FALSE)

if (length(gds) == 0) {
  stop("Dataset GDS5211 could not be downloaded. Please check internet access or GEO access.")
}

gds <- gds[[1]]

cat("Data downloaded successfully.\n")
cat("Dataset class:", class(gds), "\n")

# Extract phenotype metadata
pheno <- pData(gds)
cat("Number of samples:", nrow(pheno), "\n")
print(head(pheno[, c("title", colnames(pheno)[grep("description|characteristics|source_name|sample", tolower(colnames(pheno)))][1:min(5, length(grep("description|characteristics|source_name|sample", tolower(colnames(pheno)))))] )], 5))

# Sample names
cat("Sample names:\n")
print(sampleNames(gds))

# Save object for later steps
saveRDS(gds, file = "results/data/gds5211_raw.rds")

# Show basic feature dimension
if (methods::is(gds, "GDS")) {
  feature_table <- Table(gds)
  cat("Rows (features):", nrow(feature_table), "\n")
  cat("Columns (sample-related fields):", ncol(feature_table), "\n")
  print(head(feature_table))
} else {
  cat("Expression matrix dimensions:", dim(exprs(gds)), "\n")
}

cat("\nSTEP 1 completed. Data saved to results/data/gds5211_raw.rds\n")

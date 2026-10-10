# ============================================================================
# REQUIREMENTS.R - Install all necessary packages for DEG analysis
# ============================================================================

cat("\n=== Installing Dependencies ===")
cat("\nThis script will install required packages for GDS5211 analysis.\n")

# Install BiocManager if not available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager")
}

library(BiocManager)

# Bioconductor packages
bioc_packages <- c(
  "GEOquery",      # For downloading GEO datasets
  "limma",         # For differential expression analysis
  "affy",          # For Affymetrix microarray data
  "ComplexHeatmap" # For advanced heatmap visualization
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "from Bioconductor...\n")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    cat(pkg, "is already installed.\n")
  }
}

# CRAN packages
cran_packages <- c(
  "ggplot2",      # For plotting
  "pheatmap",     # For heatmaps
  "reshape2",     # For data reshaping
  "dplyr",        # For data manipulation
  "tidyr",        # For tidy data
  "stringr",      # For string operations
  "plotly",       # For interactive plots
  "cowplot",      # For plot arrangement
  "RColorBrewer"  # For color palettes
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "from CRAN...\n")
    install.packages(pkg, repos = "https://cran.rstudio.com", dependencies = TRUE)
  } else {
    cat(pkg, "is already installed.\n")
  }
}

cat("\n=== All dependencies installed successfully! ===")
cat("\nYou can now run the DEG analysis pipeline.\n")

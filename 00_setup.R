# ============================================================================
# 00_SETUP.R - Setup directories and initial configuration
# ============================================================================

cat("\n=== Setting up directories ===")

# Create necessary directories
dirs_to_create <- c(
  "data",
  "data/raw",
  "data/processed",
  "results",
  "results/plots",
  "results/plots/quality_control",
  "results/plots/deg_analysis",
  "results/tables",
  "results/data",
  "logs"
)

for (dir in dirs_to_create) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cat("Created directory:", dir, "\n")
  }
}

# Set seed for reproducibility
set.seed(12345)

# Configuration parameters
CONFIG <- list(
  gds_accession = "GDS5211",
  output_format = "png",
  dpi = 300,
  width = 10,
  height = 8,
  fc_threshold = 1,
  p_threshold = 0.05,
  q_threshold = 0.05,
  top_n_genes = 50,
  top_n_heatmap = 30,
  random_seed = 12345
)

cat("\n=== Configuration ===")
cat("GDS Accession:", CONFIG$gds_accession, "\n")
cat("FC Threshold:", CONFIG$fc_threshold, "\n")
cat("P-value Threshold:", CONFIG$p_threshold, "\n")
cat("Q-value Threshold:", CONFIG$q_threshold, "\n")

cat("\nSetup completed.\n")

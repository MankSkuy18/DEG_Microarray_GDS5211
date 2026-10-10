#!/usr/bin/env Rscript
# ============================================================================
# 05_FULL_PIPELINE.R - Complete DEG Analysis Workflow for GDS5211
# ============================================================================
# This script executes the entire analysis pipeline in sequence
# ============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║          GDS5211 Microarray DEG Analysis Pipeline             ║\n")
cat("║                    Complete Workflow                          ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

start_time <- Sys.time()

tryCatch({
  # Load configuration and setup
  cat("Loading configuration...\n")
  source("00_setup.R")
  
  # Check for required packages
  cat("\nChecking dependencies...\n")
  required_packages <- c("GEOquery", "limma", "ggplot2", "pheatmap", "dplyr")
  missing_packages <- required_packages[!sapply(required_packages, require, character.only = TRUE)]
  
  if (length(missing_packages) > 0) {
    cat("\nWARNING: Missing packages:", paste(missing_packages, collapse = ", "), "\n")
    cat("Please run requirements.R to install dependencies.\n")
    stop("Dependencies not installed.")
  }
  
  cat("✓ All dependencies available\n")
  
  # Execute analysis steps
  cat("\n--- Executing Analysis Pipeline ---\n\n")
  
  # Step 1: Download and Explore
  source("01_download_and_explore.R")
  
  # Step 2: Preprocessing and Normalization
  source("02_preprocessing_normalization.R")
  
  # Step 3: Differential Expression Analysis
  source("03_deg_analysis.R")
  
  # Step 4: Visualization
  source("04_visualization.R")
  
  # ========================================================================
  # Generate Summary Report
  # ========================================================================
  
  cat("\n\n========================================")
  cat("\n=== GENERATING FINAL REPORT ===")
  cat("\n========================================\n")
  
  end_time <- Sys.time()
  elapsed_time <- difftime(end_time, start_time, units = "mins")
  
  # Load final results
  deg_results <- readRDS("results/data/deg_results.rds")
  expr <- readRDS("results/data/gds5211_expression_matrix.rds")
  
  # Create summary report
  report <- paste0(
    "\n╔════════════════════════════════════════════════════════════════╗\n",
    "║                    ANALYSIS SUMMARY REPORT                    ║\n",
    "║                  GDS5211 Microarray DEG Study                 ║\n",
    "╚════════════════════════════════════════════════════════════════╝\n",
    "\n",
    "DATASET INFORMATION:\n",
    "  Accession ID: ", CONFIG$gds_accession, "\n",
    "  Source: NCBI Gene Expression Omnibus (GEO)\n",
    "  Platform: Microarray\n",
    "\n",
    "DATA DIMENSIONS:\n",
    "  Samples: ", ncol(expr), "\n",
    "  Genes analyzed: ", nrow(deg_results), "\n",
    "\n",
    "ANALYSIS PARAMETERS:\n",
    "  Log2 Fold Change threshold: > ", CONFIG$fc_threshold, "\n",
    "  P-value threshold: < ", CONFIG$p_threshold, "\n",
    "  FDR (q-value) threshold: < ", CONFIG$q_threshold, "\n",
    "\n",
    "DEG STATISTICS:\n",
    "  Total significant DEGs: ", sum(deg_results$Significant == "Yes"), "\n",
    "  Upregulated genes: ", sum(deg_results$Direction == "Upregulated"), "\n",
    "  Downregulated genes: ", sum(deg_results$Direction == "Downregulated"), "\n",
    "  Not significant: ", sum(deg_results$Direction == "Not Significant"), "\n",
    "\n",
    "OUTPUT FILES GENERATED:\n",
    "\n",
    "TABLES:\n",
    "  • results/tables/01_sample_expression_summary.csv\n",
    "  • results/tables/02_deg_results_complete.csv\n",
    "  • results/tables/03_top_upregulated_genes.csv\n",
    "  • results/tables/04_top_downregulated_genes.csv\n",
    "\n",
    "QUALITY CONTROL PLOTS:\n",
    "  • results/plots/quality_control/01_boxplot_raw_expression.png\n",
    "  • results/plots/quality_control/02_density_raw_expression.png\n",
    "  • results/plots/quality_control/03_sample_correlation.png\n",
    "\n",
    "DEG ANALYSIS PLOTS:\n",
    "  • results/plots/deg_analysis/01_volcano_plot.png\n",
    "  • results/plots/deg_analysis/02_heatmap_top_degs.png\n",
    "  • results/plots/deg_analysis/03_ma_plot.png\n",
    "  • results/plots/deg_analysis/04_boxplot_gene_*.png\n",
    "  • results/plots/deg_analysis/05_pca_plot.png\n",
    "  • results/plots/deg_analysis/06_logfc_distribution.png\n",
    "  • results/plots/deg_analysis/07_pvalue_distribution.png\n",
    "  • results/plots/deg_analysis/08_deg_summary_barplot.png\n",
    "\n",
    "DATA OBJECTS:\n",
    "  • results/data/gds5211_expression_matrix.rds\n",
    "  • results/data/gds5211_phenotype.rds\n",
    "  • results/data/deg_results.rds\n",
    "  • results/data/limma_fit_object.rds\n",
    "\n",
    "EXECUTION TIME:\n",
    "  Total time: ", round(as.numeric(elapsed_time), 2), " minutes\n",
    "\n",
    "╔════════════════════════════════════════════════════════════════╗\n",
    "║                   ANALYSIS COMPLETED                          ║\n",
    "║              All plots and tables have been saved             ║\n",
    "║                  Check results/ directory                     ║\n",
    "╚════════════════════════════════════════════════════════════════╝\n"
  )
  
  cat(report)
  
  # Save report to file
  sink("results/ANALYSIS_REPORT.txt")
  cat(report)
  sink()
  
  cat("\n✓ Report saved to: results/ANALYSIS_REPORT.txt\n")
  
}, error = function(e) {
  cat("\n\nERROR ENCOUNTERED:\n")
  cat("Message:", e$message, "\n")
  cat("\nPlease check the error message above and ensure all dependencies are installed.\n")
  stop(e)
})

cat("\nAnalysis pipeline completed.\n")

# ============================================================================
# 04_VISUALIZATION.R - Create publication-ready plots for DEG results
# ============================================================================

cat("\n========================================")
cat("\n=== STEP 4: Visualization ===")
cat("\n========================================\n")

if (!file.exists("results/data/deg_results.rds")) {
  stop("DEG results not found. Please run 03_deg_analysis.R first.")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(limma)
  library(dplyr)
  library(RColorBrewer)
})

cat("Loading analysis data...\n")
expr <- readRDS("results/data/gds5211_expression_matrix.rds")
deg_results <- readRDS("results/data/deg_results.rds")
fit2 <- readRDS("results/data/limma_fit_object.rds")

# ============================================================================
# 1. VOLCANO PLOT
# ============================================================================

cat("Creating Volcano Plot...\n")

volcano_df <- deg_results %>%
  mutate(
    color = case_when(
      Direction == "Upregulated" ~ "Upregulated",
      Direction == "Downregulated" ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

p_volcano <- ggplot(volcano_df, aes(x = logFC, y = -log10(P.Value), color = color)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_vline(xintercept = c(-CONFIG$fc_threshold, CONFIG$fc_threshold),
             linetype = "dashed",
             color = "gray50",
             linewidth = 0.8) +
  geom_hline(yintercept = -log10(CONFIG$p_threshold),
             linetype = "dashed",
             color = "gray50",
             linewidth = 0.8) +
  scale_color_manual(
    values = c(
      "Upregulated" = "#e41a1c",
      "Downregulated" = "#377eb8",
      "Not Significant" = "#cccccc"
    )
  ) +
  labs(
    title = "Volcano Plot: Differential Expression Analysis",
    subtitle = paste0("GDS5211 - FDR < ", CONFIG$q_threshold, ", |logFC| > ", CONFIG$fc_threshold),
    x = "log2 Fold Change",
    y = "-log10(P-value)",
    color = "Gene Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "right"
  )

ggsave(
  "results/plots/deg_analysis/01_volcano_plot.png",
  p_volcano,
  width = CONFIG$width,
  height = CONFIG$height,
  dpi = CONFIG$dpi
)
cat("✓ Saved: 01_volcano_plot.png\n")

# ============================================================================
# 2. HEATMAP OF TOP DEGs
# ============================================================================

cat("Creating Heatmap of Top DEGs...\n")

sig_genes <- deg_results %>%
  filter(Significant == "Yes") %>%
  arrange(P.Value) %>%
  slice(1:min(CONFIG$top_n_heatmap, n()))

if (nrow(sig_genes) > 0) {
  top_genes <- sig_genes$Gene
  heatmat <- expr[rownames(expr) %in% top_genes, , drop = FALSE]
  
  if (nrow(heatmat) > 0) {
    # Z-score normalization
    heatmat_scaled <- t(scale(t(heatmat)))
    
    # Create annotation
    anno_col <- data.frame(
      Direction = sig_genes$Direction[match(rownames(heatmat_scaled), sig_genes$Gene)]
    )
    rownames(anno_col) <- rownames(heatmat_scaled)
    
    anno_colors <- list(
      Direction = c(
        "Upregulated" = "#e41a1c",
        "Downregulated" = "#377eb8",
        "Not Significant" = "#cccccc"
      )
    )
    
    png(
      "results/plots/deg_analysis/02_heatmap_top_degs.png",
      width = CONFIG$width * 120,
      height = CONFIG$height * 120,
      res = CONFIG$dpi
    )
    pheatmap::pheatmap(
      heatmat_scaled,
      scale = "none",
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "ward.D2",
      main = paste("Heatmap of Top", nrow(heatmat_scaled), "DEGs (Z-score normalized)"),
      color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
      border_color = NA,
      fontsize = 10,
      breaks = seq(-3, 3, length.out = 101)
    )
    dev.off()
    cat("✓ Saved: 02_heatmap_top_degs.png\n")
  }
}

# ============================================================================
# 3. MA PLOT
# ============================================================================

cat("Creating MA Plot...\n")

ma_df <- deg_results %>%
  rename(M = logFC, A = AveExpr) %>%
  mutate(color = Direction)

p_ma <- ggplot(ma_df, aes(x = A, y = M, color = color)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(-CONFIG$fc_threshold, CONFIG$fc_threshold),
             linetype = "dashed",
             color = "gray50",
             linewidth = 0.8) +
  scale_color_manual(
    values = c(
      "Upregulated" = "#e41a1c",
      "Downregulated" = "#377eb8",
      "Not Significant" = "#cccccc"
    )
  ) +
  labs(
    title = "MA Plot: Mean vs Fold Change",
    x = "Average log2(Expression)",
    y = "log2(Fold Change)",
    color = "Gene Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

ggsave(
  "results/plots/deg_analysis/03_ma_plot.png",
  p_ma,
  width = CONFIG$width,
  height = CONFIG$height,
  dpi = CONFIG$dpi
)
cat("✓ Saved: 03_ma_plot.png\n")

# ============================================================================
# 4. BOXPLOTS FOR TOP GENES
# ============================================================================

cat("Creating Boxplots for Top Genes...\n")

top_sig_genes <- deg_results %>%
  filter(Significant == "Yes") %>%
  arrange(P.Value) %>%
  slice(1:min(6, n()))

if (nrow(top_sig_genes) > 0) {
  for (i in 1:nrow(top_sig_genes)) {
    gene <- top_sig_genes$Gene[i]
    gene_expr <- expr[rownames(expr) == gene, ]
    
    if (length(gene_expr) > 0) {
      # Create groups for boxplot
      group <- factor(rep(c("Group1", "Group2"), length.out = length(gene_expr)))
      box_df <- data.frame(
        Gene = gene,
        Expression = as.numeric(gene_expr),
        Group = group
      )
      
      p_box <- ggplot(box_df, aes(x = Group, y = Expression, fill = Group)) +
        geom_boxplot(width = 0.5, alpha = 0.7, outlier.size = 2) +
        geom_jitter(width = 0.15, size = 2.5, alpha = 0.6) +
        labs(
          title = paste("Gene:", gene),
          subtitle = paste0("logFC = ", round(top_sig_genes$logFC[i], 2),
                           ", adj.P.Val = ", format(top_sig_genes$adj.P.Val[i], scientific = TRUE)),
          x = "Group",
          y = "log2(Expression + 1)"
        ) +
        scale_fill_manual(values = c("Group1" = "#377eb8", "Group2" = "#e41a1c")) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = "none"
        )
      
      filename <- paste0(
        "results/plots/deg_analysis/04_boxplot_gene_",
        i,
        "_",
        gsub("\\W+", "_", gene),
        ".png"
      )
      ggsave(filename, p_box, width = 7, height = 6, dpi = CONFIG$dpi)
      cat("✓ Saved:", basename(filename), "\n")
    }
  }
}

# ============================================================================
# 5. PCA PLOT
# ============================================================================

cat("Creating PCA Plot...\n")

pca <- prcomp(t(expr), center = TRUE, scale. = TRUE)
var_explained <- summary(pca)$importance[2, 1:2] * 100

pca_df <- as.data.frame(pca$x[, 1:2])
pca_df$Group <- factor(rep(c("Group1", "Group2"), length.out = nrow(pca_df)))

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_ellipse(
    aes(fill = Group),
    alpha = 0.1,
    linewidth = 1
  ) +
  scale_color_manual(values = c("Group1" = "#377eb8", "Group2" = "#e41a1c")) +
  scale_fill_manual(values = c("Group1" = "#377eb8", "Group2" = "#e41a1c")) +
  scale_shape_manual(values = c(16, 17)) +
  labs(
    title = "Principal Component Analysis",
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance)"),
    color = "Group",
    shape = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

ggsave(
  "results/plots/deg_analysis/05_pca_plot.png",
  p_pca,
  width = CONFIG$width,
  height = CONFIG$height,
  dpi = CONFIG$dpi
)
cat("✓ Saved: 05_pca_plot.png\n")

# ============================================================================
# 6. DISTRIBUTION PLOTS
# ============================================================================

cat("Creating distribution plots...\n")

# logFC distribution
p_fc_dist <- ggplot(deg_results, aes(x = logFC, fill = Direction)) +
  geom_histogram(bins = 50, alpha = 0.7, color = "black") +
  scale_fill_manual(
    values = c(
      "Upregulated" = "#e41a1c",
      "Downregulated" = "#377eb8",
      "Not Significant" = "#cccccc"
    )
  ) +
  labs(
    title = "Distribution of log2 Fold Changes",
    x = "log2 Fold Change",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(
  "results/plots/deg_analysis/06_logfc_distribution.png",
  p_fc_dist,
  width = CONFIG$width,
  height = CONFIG$height,
  dpi = CONFIG$dpi
)
cat("✓ Saved: 06_logfc_distribution.png\n")

# -log10(P-value) distribution
p_pval_dist <- ggplot(deg_results, aes(x = -log10(P.Value), fill = Direction)) +
  geom_histogram(bins = 50, alpha = 0.7, color = "black") +
  scale_fill_manual(
    values = c(
      "Upregulated" = "#e41a1c",
      "Downregulated" = "#377eb8",
      "Not Significant" = "#cccccc"
    )
  ) +
  labs(
    title = "Distribution of -log10(P-values)",
    x = "-log10(P-value)",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(
  "results/plots/deg_analysis/07_pvalue_distribution.png",
  p_pval_dist,
  width = CONFIG$width,
  height = CONFIG$height,
  dpi = CONFIG$dpi
)
cat("✓ Saved: 07_pvalue_distribution.png\n")

# ============================================================================
# 7. SUMMARY STATISTICS PLOT
# ============================================================================

cat("Creating summary statistics plot...\n")

summary_stats <- deg_results %>%
  summarize(
    Total_Genes = n(),
    Upregulated = sum(Direction == "Upregulated"),
    Downregulated = sum(Direction == "Downregulated"),
    Not_Significant = sum(Direction == "Not Significant")
  ) %>%
  pivot_longer(everything(), names_to = "Category", values_to = "Count")

p_summary <- ggplot(summary_stats, aes(x = reorder(Category, -Count), y = Count, fill = Category)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 1) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4) +
  scale_fill_manual(
    values = c(
      "Total_Genes" = "#999999",
      "Upregulated" = "#e41a1c",
      "Downregulated" = "#377eb8",
      "Not_Significant" = "#cccccc"
    )
  ) +
  labs(
    title = "DEG Summary Statistics",
    x = "",
    y = "Number of Genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(
  "results/plots/deg_analysis/08_deg_summary_barplot.png",
  p_summary,
  width = CONFIG$width,
  height = CONFIG$height,
  dpi = CONFIG$dpi
)
cat("✓ Saved: 08_deg_summary_barplot.png\n")

cat("\n✓ All visualization plots saved to results/plots/\n")

cat("\n========================================")
cat("\n=== STEP 4 COMPLETE ===")
cat("\n========================================\n")

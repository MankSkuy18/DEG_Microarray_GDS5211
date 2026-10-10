# 04_visualization.R
# Visualization of DEGs for GDS5211

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(limma)
  library(dplyr)
})

# Load processed expression and DEG results
if (!file.exists("results/data/gds5211_processed.rds")) {
  stop("Processed expression data not found. Run 02_preprocessing_normalization.R first.")
}
if (!file.exists("results/data/deg_results_gds5211.rds")) {
  stop("DEG results not found. Run 03_deg_analysis.R first.")
}

expr <- readRDS("results/data/gds5211_processed.rds")
res <- readRDS("results/data/deg_results_gds5211.rds")

# 1. Volcano Plot
volcano_df <- res %>%
  mutate(sig = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Significant", "Not Significant"))

p_volcano <- ggplot(volcano_df, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  labs(
    title = "Volcano Plot - GDS5211",
    x = "log2 Fold Change",
    y = "-log10(P-value)",
    color = "Status"
  ) +
  theme_minimal(base_size = 12)

ggsave("results/plots/volcano_plot.png", p_volcano, width = 8, height = 6, dpi = 300)

# 2. Heatmap of top DEGs
sig_genes <- res %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  top_n(50, wt = -adj.P.Val)

if (nrow(sig_genes) > 0) {
  top_genes <- sig_genes$Gene
  heatmat <- expr[top_genes, , drop = FALSE]
  if (nrow(heatmat) > 1 && ncol(heatmat) > 1) {
    heatmat <- t(scale(t(heatmat)))
    png("results/plots/heatmap_top_genes.png", width = 1200, height = 1000, res = 150)
      pheatmap::pheatmap(
        heatmat,
        scale = "none",
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        main = "Top DEGs Heatmap",
        color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
      )
    dev.off()
  }
}

# 3. Boxplot for most significant gene
if (nrow(res) > 0) {
  top_gene <- res %>% arrange(P.Value) %>% slice(1) %>% pull(Gene)
  gene_expr <- as.numeric(expr[top_gene, ])
  group <- factor(rep(c("Group1", "Group2"), length.out = length(gene_expr)))
  box_df <- data.frame(Gene = top_gene, Expression = gene_expr, Group = group)
  p_box <- ggplot(box_df, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(width = 0.5, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1.5) +
    labs(title = paste("Boxplot for", top_gene), x = "Group", y = "Expression") +
    theme_minimal(base_size = 12)
  ggsave("results/plots/boxplot_top_gene.png", p_box, width = 7, height = 6, dpi = 300)
}

# 4. PCA plot
pca <- prcomp(t(expr), center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2])
pca_df$Group <- factor(rep(c("Group1", "Group2"), length.out = nrow(pca_df)))

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA of GDS5211 Samples", x = "PC1", y = "PC2") +
  theme_minimal(base_size = 12)
ggsave("results/plots/pca_plot.png", p_pca, width = 7, height = 6, dpi = 300)

# 5. MA plot
png("results/plots/ma_plot.png", width = 1000, height = 700, res = 150)
  limma::plotMA(
    as.matrix(expr),
    main = "MA Plot",
    cex.main = 1.2
  )
dev.off()

# 6. Density plot
p_density <- ggplot(data.frame(Value = as.vector(expr)), aes(x = Value)) +
  geom_density(alpha = 0.5, fill = "steelblue") +
  labs(title = "Expression Density", x = "Expression Value", y = "Density") +
  theme_minimal(base_size = 12)
ggsave("results/plots/density_plot.png", p_density, width = 7, height = 6, dpi = 300)

cat("Visualization plots saved in results/plots\n")

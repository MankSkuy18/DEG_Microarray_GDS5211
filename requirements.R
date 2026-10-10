# Install required packages

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages <- c(
  "GEOquery",
  "limma",
  "ggplot2",
  "pheatmap",
  "reshape2",
  "dplyr",
  "tidyr",
  "stringr",
  "stats",
  "utils"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("GEOquery", "limma", "pheatmap")) {
      BiocManager::install(pkg, ask = FALSE)
    } else {
      install.packages(pkg, repos = "https://cran.rstudio.com")
    }
  }
}

message("Dependencies installed successfully.")

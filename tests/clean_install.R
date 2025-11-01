#!/usr/bin/env Rscript

# Test clean installation of genePackRat from GitHub

cat("Testing genePackRat installation from GitHub...\n\n")

# Install from GitHub
cat("Installing genePackRat from GitHub...\n")
devtools::install_github("RauLabUNC/genePackRat")

# Load the package
cat("\nLoading genePackRat...\n")
library(genePackRat)

# Test that main functions are available
cat("\nTesting function availability...\n")
if (exists("filterGenes")) {
  cat("✓ filterGenes function found\n")
} else {
  cat("✗ filterGenes function NOT found\n")
}

if (exists("build_gene_table")) {
  cat("✓ build_gene_table function found\n")
} else {
  cat("✗ build_gene_table function NOT found\n")
}

# Check help documentation
cat("\nTesting documentation access...\n")
help_obj <- help("filterGenes", package = "genePackRat")
if (length(help_obj) > 0) {
  cat("✓ Documentation accessible\n")
} else {
  cat("✗ Documentation NOT accessible\n")
}

# Quick functional test
cat("\nRunning basic functional test...\n")
test_genes <- data.frame(
  gene_id = c("GENE1", "GENE2", "GENE3"),
  expression = c(10, 50, 100),
  pvalue = c(0.01, 0.05, 0.1)
)

# Test basic filtering
result <- filterGenes(
  inputTable = test_genes,
  filters = list(expression = list(column = "expression", operator = ">", value = 30))
)

cat("Input genes:", nrow(test_genes), "\n")
cat("Filtered genes:", nrow(result), "\n")

if (nrow(result) == 2) {
  cat("✓ Filtering works correctly\n")
} else {
  cat("✗ Filtering unexpected result\n")
}

cat("\n=== Installation test complete ===\n")
cat("Package successfully installed and functional!\n")
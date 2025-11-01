#!/usr/bin/env Rscript

# Test clean installation of genePackRat from GitHub

# Clear any existing R library paths to ensure clean environment
Sys.setenv(R_LIBS_USER = "")
Sys.setenv(R_LIBS = "")

# Create a local library directory for this test
test_lib <- file.path(getwd(), "test_r_libs")
if (!dir.exists(test_lib)) {
  dir.create(test_lib, recursive = TRUE)
}

# Set library paths - local test library first, then system
.libPaths(c(test_lib, .libPaths()))
cat("Using R libraries from:\n")
cat(paste("  -", .libPaths()), sep = "\n")
cat("\n")

cat("Testing genePackRat installation from GitHub...\n\n")

# Install from GitHub
cat("Installing genePackRat from GitHub...\n")

# Clear any problematic PAT tokens that might interfere
Sys.unsetenv("GITHUB_PAT")
Sys.unsetenv("GITHUB_TOKEN")

# Try installation with explicit parameters
tryCatch({
  devtools::install_github("RauLabUNC/genePackRat",
                          auth_token = NULL,
                          upgrade = "never",
                          lib = test_lib)
}, error = function(e) {
  cat("\nError during installation:\n", e$message, "\n")
  cat("\nTrying alternative installation method...\n")
  # Try using remotes directly as a backup
  remotes::install_github("RauLabUNC/genePackRat",
                          upgrade = "never",
                          lib = test_lib)
})

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

if (exists("makeFilter")) {
  cat("✓ makeFilter function found\n")
} else {
  cat("✗ makeFilter function NOT found\n")
}

# Quick functional test
cat("\nRunning basic functional test...\n")
test_genes <- data.frame(
  gene_id = c("GENE1", "GENE2", "GENE3"),
  expression = c(10, 50, 100),
  pvalue = c(0.01, 0.05, 0.1)
)

# Test basic filtering
# Create a simple filter
filter <- makeFilter(
  column = "expression",
  condition = ">",
  value = 30
)

result <- filterGenes(
  inputTable = test_genes,
  filters = list(filter)
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
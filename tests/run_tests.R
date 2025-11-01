#!/usr/bin/env Rscript

# Run genePackRat functionality tests
# This script ensures proper library paths and runs comprehensive tests

# Clear and set library path
Sys.setenv(R_LIBS_USER = "")
Sys.setenv(R_LIBS = "")

# Use test_r_libs if it exists (where we installed the package)
if (dir.exists("test_r_libs")) {
  .libPaths(c("test_r_libs", .libPaths()))
}

cat("========================================\n")
cat("genePackRat Functionality Tests\n")
cat("========================================\n\n")

cat("Library paths:\n")
cat(paste("  -", .libPaths()), sep = "\n")
cat("\n")

# Load the package
library(genePackRat)
cat("✓ Package loaded successfully\n")
cat("Package version:", as.character(packageVersion("genePackRat")), "\n\n")

# Load test data
genes <- read.csv("../inst/testdata/mini_genes.csv", stringsAsFactors = FALSE)
expression <- read.csv("../inst/testdata/mini_expression.csv", stringsAsFactors = FALSE)
gwas <- read.csv("../inst/testdata/mini_gwas.csv", stringsAsFactors = FALSE)

cat("✓ Test data loaded\n")
cat("  - Genes:", nrow(genes), "rows\n")
cat("  - Expression:", nrow(expression), "rows\n")
cat("  - GWAS:", nrow(gwas), "rows\n\n")

# Run tests
test_results <- list()

# Test 1: Basic filtering
cat("Test 1: Basic Filtering\n")
cat("-----------------------\n")
filter1 <- makeFilter(column = "biotype", condition = "==", value = "protein_coding")
result1 <- filterGenes(inputTable = genes, filters = list(filter1))
test_results$basic <- nrow(result1) > 0 && all(result1$biotype == "protein_coding")
cat("  Protein coding filter:", ifelse(test_results$basic, "✓ PASSED", "✗ FAILED"), "\n")
cat("  Found", nrow(result1), "protein coding genes\n\n")

# Test 2: Join with expression data
cat("Test 2: Data Integration\n")
cat("------------------------\n")
gene_table <- build_gene_table(genes_df = genes, expression_df = expression)
test_results$integration <- ncol(gene_table) > ncol(genes)
cat("  Build gene table:", ifelse(test_results$integration, "✓ PASSED", "✗ FAILED"), "\n")
cat("  Columns added:", ncol(gene_table) - ncol(genes), "\n\n")

# Test 3: Numeric filtering on integrated data
cat("Test 3: Expression Filtering\n")
cat("----------------------------\n")
filter_expr <- makeFilter(column = "TPM_brain", condition = ">", value = 100)
high_expr <- filterGenes(inputTable = gene_table, filters = list(filter_expr))
test_results$expression <- nrow(high_expr) > 0 && all(high_expr$TPM_brain > 100)
cat("  High brain expression:", ifelse(test_results$expression, "✓ PASSED", "✗ FAILED"), "\n")
cat("  Found", nrow(high_expr), "genes with TPM_brain > 100\n\n")

# Test 4: Relational filtering
cat("Test 4: Relational Filtering\n")
cat("----------------------------\n")

# Semi join - genes with GWAS hits
gwas_genes <- filterGenes(
  inputTable = genes,
  referenceTable = gwas,
  by = "gene_id",
  joinType = "semi"
)
test_results$semi_join <- nrow(gwas_genes) > 0 && ncol(gwas_genes) == ncol(genes)
cat("  Semi join:", ifelse(test_results$semi_join, "✓ PASSED", "✗ FAILED"), "\n")
cat("  Found", nrow(gwas_genes), "genes with GWAS associations\n")

# Anti join - genes without GWAS hits
no_gwas <- filterGenes(
  inputTable = genes,
  referenceTable = gwas,
  by = "gene_id",
  joinType = "anti"
)
test_results$anti_join <- nrow(no_gwas) + nrow(gwas_genes) == nrow(genes)
cat("  Anti join:", ifelse(test_results$anti_join, "✓ PASSED", "✗ FAILED"), "\n")
cat("  Found", nrow(no_gwas), "genes without GWAS associations\n\n")

# Test 5: Multiple filters
cat("Test 5: Multiple Filters\n")
cat("------------------------\n")
filter_bio <- makeFilter(column = "biotype", condition = "==", value = "protein_coding")
filter_sig <- makeFilter(column = "padj", condition = "<", value = 0.05)
multi_result <- filterGenes(
  inputTable = gene_table,
  filters = list(filter_bio, filter_sig)
)
test_results$multiple <- all(multi_result$biotype == "protein_coding") &&
                         all(multi_result$padj < 0.05)
cat("  Multiple filters:", ifelse(test_results$multiple, "✓ PASSED", "✗ FAILED"), "\n")
cat("  Genes passing both filters:", nrow(multi_result), "\n\n")

# Test 6: Excel export
cat("Test 6: Excel Export\n")
cat("--------------------\n")
output_file <- "test_output.xlsx"
tryCatch({
  wb <- create_gene_workbook(
    gene_data = gene_table,
    output_file = output_file
  )
  test_results$excel <- file.exists(output_file)
  if (test_results$excel) {
    file_size <- file.info(output_file)$size
    cat("  Excel creation:", "✓ PASSED\n")
    cat("  File size:", file_size, "bytes\n")
    unlink(output_file)
  } else {
    cat("  Excel creation:", "✗ FAILED\n")
  }
}, error = function(e) {
  cat("  Excel creation: ✗ FAILED\n")
  cat("  Error:", e$message, "\n")
  test_results$excel <- FALSE
})

cat("\n")
cat("========================================\n")
cat("Test Summary\n")
cat("========================================\n")

# Summary
total <- length(test_results)
passed <- sum(unlist(test_results))
failed <- total - passed

for (name in names(test_results)) {
  status <- ifelse(test_results[[name]], "✓", "✗")
  cat(sprintf("  %-15s %s\n", paste0(name, ":"), status))
}

cat("\n")
cat("Results: ", passed, "/", total, " tests passed\n", sep = "")

if (passed == total) {
  cat("\n✓ All tests passed successfully!\n")
  quit(status = 0)
} else {
  cat("\n✗", failed, "test(s) failed\n")
  quit(status = 1)
}
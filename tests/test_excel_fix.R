#!/usr/bin/env Rscript

# Test Excel fix by loading from source

# Load required libraries
library(data.table)
library(dplyr)
library(openxlsx)

# Source the updated Excel file
source("../R/excel.R")
source("../R/filterGenes.R")
source("../R/join_tables.R")

# Load test data
genes <- read.csv("../inst/testdata/mini_genes.csv", stringsAsFactors = FALSE)
expression <- read.csv("../inst/testdata/mini_expression.csv", stringsAsFactors = FALSE)

cat("Testing Excel export with fixed function...\n\n")

# Build gene table
gene_table <- build_gene_table(
  genes_df = genes,
  expression_df = expression
)

cat("Gene table created with columns:\n")
cat(paste("  -", names(gene_table)), sep = "\n")
cat("\n")

# Test Excel export
output_file <- "test_excel_output.xlsx"

tryCatch({
  wb <- create_gene_workbook(
    gene_data = gene_table,
    output_file = output_file
  )

  if (file.exists(output_file)) {
    file_info <- file.info(output_file)
    cat("\n✓ Excel file created successfully!\n")
    cat("  File:", output_file, "\n")
    cat("  Size:", file_info$size, "bytes\n")

    # Check the workbook structure
    wb_check <- openxlsx::loadWorkbook(output_file)
    sheets <- names(wb_check)
    cat("\n  Sheets created:\n")
    cat(paste("    -", sheets), sep = "\n")

    # Clean up
    unlink(output_file)
    cat("\n✓ Test passed! Excel export is working.\n")
  } else {
    cat("\n✗ Excel file was not created\n")
  }
}, error = function(e) {
  cat("\n✗ Excel creation failed:\n")
  cat("  Error:", e$message, "\n")
  cat("\n  Traceback:\n")
  traceback()
})
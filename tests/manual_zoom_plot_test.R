#!/usr/bin/env Rscript
# Manual test script for generateLocusZoomPlot with your working setup

cat("=== Testing generateLocusZoomPlot ===\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(plotgardener)
  library(RColorBrewer)
})

# Source the function
source("/proj/raulab/users/brian/packrat/R/packet_core.R")

# Load test data
cat("Loading data...\n")
# Use subset test data stored in repository
scan_data <- readRDS("/proj/raulab/users/brian/packrat/tests/testthat/fixtures/test_scan_data_subset.rds")
threshold_data <- readRDS("/proj/raulab/users/brian/packrat/tests/testthat/fixtures/test_threshold_data_subset.rds")
genes <- fread("/proj/raulab/users/brian/packrat/tests/testthat/fixtures/test_genes_mouse.csv")

# Define locus (using your working parameters)
locus_info <- data.frame(
  chr = 1,
  peak_pos = 101.0,
  start_pos = 90,  # Wider region as in your example
  end_pos = 102,
  trait = 'BW.day.0',
  drug = 'Ctrl',
  max_lod = 7.0
)

cat("Locus: chr", locus_info$chr, " from", locus_info$start_pos, "to", locus_info$end_pos, "Mb\n")
cat("Trait:", locus_info$trait, "Drug:", locus_info$drug, "\n")

# Get genes in locus
genes_in_locus <- genes[chr == 1 &
                        start_bp < locus_info$end_pos * 1e6 &
                        end_bp > locus_info$start_pos * 1e6]
cat("Found", nrow(genes_in_locus), "genes in region\n")

# Select top genes (first 10 as in your example)
n_top_genes <- min(10, nrow(genes_in_locus))
if (n_top_genes > 0) {
  top_genes <- data.frame(
    gene = genes_in_locus$mouse_gene_symbol[1:n_top_genes],
    color = '#e34a33',
    stringsAsFactors = FALSE
  )
  cat("Highlighting", n_top_genes, "genes\n")
} else {
  top_genes <- data.frame(gene = character(0), color = character(0))
  cat("No genes to highlight\n")
}

# Single overlap (as you noted this part needs work)
overlapping_loci <- data.frame(
  chrom = 'chr1',
  start = locus_info$start_pos * 1e6,
  end = locus_info$end_pos * 1e6,
  strand = '-',
  trait = 'BW.day.0',
  drug = 'Ctrl',
  traitXdrug = 'BW.day.0: Ctrl',
  stringsAsFactors = FALSE
)

# Test 1: Basic plot generation
cat("\n--- Test 1: Basic plot generation ---\n")
output_file <- 'tests/test_plot_basic.pdf'

tryCatch({
  result <- generateLocusZoomPlot(
    locus_info = locus_info,
    scan_data = scan_data,
    threshold_data = threshold_data,
    genes_in_locus = genes_in_locus,
    top_genes_in_locus = top_genes,
    overlapping_loci = overlapping_loci,
    output_file = output_file,
    assembly = NULL,  # Use default
    plot_params = NULL  # Use defaults
  )

  if (file.exists(output_file)) {
    file_size_kb <- round(file.size(output_file) / 1024, 1)
    cat("✓ SUCCESS: Plot generated -", output_file, "(", file_size_kb, "KB)\n")
  } else {
    cat("✗ FAIL: Output file not created\n")
  }
}, error = function(e) {
  cat("✗ ERROR:", conditionMessage(e), "\n")
})

# Test 2: Multiple overlapping loci
cat("\n--- Test 2: Multiple overlapping loci ---\n")
overlapping_loci_multi <- data.frame(
  chrom = c('chr1', 'chr1', 'chr1'),
  start = c(90e6, 95e6, 98e6),
  end = c(95e6, 100e6, 102e6),
  strand = rep('-', 3),
  trait = c('BW.day.0', 'HR.0', 'EF.0'),
  drug = c('Ctrl', 'Ctrl', 'Iso'),
  traitXdrug = c('BW.day.0: Ctrl', 'HR.0: Ctrl', 'EF.0: Iso'),
  stringsAsFactors = FALSE
)

output_file2 <- 'tests/test_plot_multi_overlap.pdf'

tryCatch({
  result <- generateLocusZoomPlot(
    locus_info = locus_info,
    scan_data = scan_data,
    threshold_data = threshold_data,
    genes_in_locus = genes_in_locus,
    top_genes_in_locus = top_genes,
    overlapping_loci = overlapping_loci_multi,
    output_file = output_file2,
    assembly = NULL,
    plot_params = NULL
  )

  if (file.exists(output_file2)) {
    file_size_kb <- round(file.size(output_file2) / 1024, 1)
    cat("✓ SUCCESS: Multi-overlap plot -", output_file2, "(", file_size_kb, "KB)\n")
  } else {
    cat("✗ FAIL: Output file not created\n")
  }
}, error = function(e) {
  cat("✗ ERROR:", conditionMessage(e), "\n")
})

# Test 3: Different chromosome and region
cat("\n--- Test 3: Different region (chr 2) ---\n")
locus_info2 <- data.frame(
  chr = 2,
  peak_pos = 50.0,
  start_pos = 48,
  end_pos = 52,
  trait = 'HR.0',
  drug = 'Ctrl',
  max_lod = 6.5
)

# Use empty genes for chr 2 test
genes_in_locus2 <- data.table()
top_genes2 <- data.frame(gene = character(0), color = character(0))
overlapping_loci2 <- data.frame(
  chrom = 'chr2',
  start = 48e6,
  end = 52e6,
  strand = '-',
  trait = 'HR.0',
  drug = 'Ctrl',
  traitXdrug = 'HR.0: Ctrl',
  stringsAsFactors = FALSE
)

output_file3 <- 'tests/test_plot_chr2.pdf'

tryCatch({
  # Check if this trait exists in the data
  if (!"HR.0_Ctrl" %in% names(scan_data)) {
    cat("✗ SKIP: HR.0_Ctrl not in scan_data\n")
  } else {
    result <- generateLocusZoomPlot(
      locus_info = locus_info2,
      scan_data = scan_data,
      threshold_data = threshold_data,
      genes_in_locus = genes_in_locus2,
      top_genes_in_locus = top_genes2,
      overlapping_loci = overlapping_loci2,
      output_file = output_file3,
      assembly = NULL,
      plot_params = NULL
    )

    if (file.exists(output_file3)) {
      file_size_kb <- round(file.size(output_file3) / 1024, 1)
      cat("✓ SUCCESS: Chr2 plot -", output_file3, "(", file_size_kb, "KB)\n")
    } else {
      cat("✗ FAIL: Output file not created\n")
    }
  }
}, error = function(e) {
  cat("✗ ERROR:", conditionMessage(e), "\n")
})

cat("\n=== Testing complete ===\n")
cat("Check the PDF files in tests/ directory\n")
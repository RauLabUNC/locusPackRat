test_that("generateLocusZoomPlot produces a PDF with real data", {
  # Load required libraries
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(plotgardener)
    library(RColorBrewer)
  })

  # Source the function directly since package isn't installed yet
  source("/proj/raulab/users/brian/packrat/R/packet_core.R")

  # Load test data from fixtures
  scan_data <- readRDS(test_path("fixtures/test_scan_data_subset.rds"))
  threshold_data <- readRDS(test_path("fixtures/test_threshold_data_subset.rds"))
  genes <- fread(test_path("fixtures/test_genes_mouse.csv"))

  # Define test locus
  locus_info <- data.frame(
    chr = 1,
    peak_pos = 101.0,
    start_pos = 90,
    end_pos = 102,
    trait = 'BW.day.0',
    drug = 'Ctrl',
    max_lod = 7.0
  )

  # Get genes in locus region
  genes_in_locus <- genes[chr == 1 &
                          start_bp < locus_info$end_pos * 1e6 &
                          end_bp > locus_info$start_pos * 1e6]

  # Select top genes (first 10)
  top_genes <- data.frame(
    gene = genes_in_locus$mouse_gene_symbol[1:min(10, nrow(genes_in_locus))],
    color = '#e34a33',
    stringsAsFactors = FALSE
  )

  # Create overlapping loci data
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

  # Define output file in temp directory
  output_file <- tempfile("test_locus_zoom_", fileext = ".pdf")

  # Generate the plot
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

  # Check that output was created and is a reasonable size
  expect_true(file.exists(output_file))
  expect_gt(file.size(output_file), 5000)  # At least 5KB
  expect_equal(result, output_file)

  # Clean up
  unlink(output_file)
})

test_that("generateLocusZoomPlot handles multiple overlapping loci", {
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(plotgardener)
    library(RColorBrewer)
  })

  # Source the function directly since package isn't installed yet
  source("/proj/raulab/users/brian/packrat/R/packet_core.R")

  # Load test data from fixtures
  scan_data <- readRDS(test_path("fixtures/test_scan_data_subset.rds"))
  threshold_data <- readRDS(test_path("fixtures/test_threshold_data_subset.rds"))
  genes <- fread(test_path("fixtures/test_genes_mouse.csv"))

  locus_info <- data.frame(
    chr = 1, peak_pos = 101.0, start_pos = 90, end_pos = 102,
    trait = 'BW.day.0', drug = 'Ctrl', max_lod = 7.0
  )

  genes_in_locus <- genes[chr == 1 &
                          start_bp < locus_info$end_pos * 1e6 &
                          end_bp > locus_info$start_pos * 1e6]

  top_genes <- data.frame(
    gene = genes_in_locus$mouse_gene_symbol[1:min(10, nrow(genes_in_locus))],
    color = '#e34a33',
    stringsAsFactors = FALSE
  )

  # Test with multiple overlapping loci
  overlapping_loci <- data.frame(
    chrom = c('chr1', 'chr1', 'chr1'),
    start = c(90e6, 95e6, 98e6),
    end = c(95e6, 100e6, 102e6),
    strand = rep('-', 3),
    trait = c('BW.day.0', 'HR.0', 'EF.0'),
    drug = c('Ctrl', 'Ctrl', 'Iso'),
    traitXdrug = c('BW.day.0: Ctrl', 'HR.0: Ctrl', 'EF.0: Iso'),
    stringsAsFactors = FALSE
  )

  output_file <- tempfile("test_multi_overlap_", fileext = ".pdf")

  result <- generateLocusZoomPlot(
    locus_info = locus_info,
    scan_data = scan_data,
    threshold_data = threshold_data,
    genes_in_locus = genes_in_locus,
    top_genes_in_locus = top_genes,
    overlapping_loci = overlapping_loci,
    output_file = output_file,
    assembly = NULL,
    plot_params = NULL
  )

  expect_true(file.exists(output_file))
  expect_gt(file.size(output_file), 5000)

  unlink(output_file)
})

test_that("generateLocusZoomPlot works with no genes in region", {
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(plotgardener)
    library(RColorBrewer)
  })

  # Source the function directly since package isn't installed yet
  source("/proj/raulab/users/brian/packrat/R/packet_core.R")

  # Load test data from fixtures
  scan_data <- readRDS(test_path("fixtures/test_scan_data_subset.rds"))
  threshold_data <- readRDS(test_path("fixtures/test_threshold_data_subset.rds"))

  locus_info <- data.frame(
    chr = 1, peak_pos = 101.0, start_pos = 90, end_pos = 102,
    trait = 'BW.day.0', drug = 'Ctrl', max_lod = 7.0
  )

  # Empty genes
  genes_in_locus <- data.table()
  top_genes <- data.frame(gene = character(0), color = character(0))

  overlapping_loci <- data.frame(
    chrom = 'chr1', start = 90e6, end = 102e6, strand = '-',
    trait = 'BW.day.0', drug = 'Ctrl', traitXdrug = 'BW.day.0: Ctrl',
    stringsAsFactors = FALSE
  )

  output_file <- tempfile("test_no_genes_", fileext = ".pdf")

  # Should still generate a plot without genes
  result <- generateLocusZoomPlot(
    locus_info = locus_info,
    scan_data = scan_data,
    threshold_data = threshold_data,
    genes_in_locus = genes_in_locus,
    top_genes_in_locus = top_genes,
    overlapping_loci = overlapping_loci,
    output_file = output_file,
    assembly = NULL,
    plot_params = NULL
  )

  expect_true(file.exists(output_file))
  expect_gt(file.size(output_file), 5000)

  unlink(output_file)
})
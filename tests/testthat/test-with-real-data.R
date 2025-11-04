test_that("real scan data from paper has correct structure", {
  skip_if_not(
    file.exists("/proj/raulab/users/brian/cc_gwas/data/processed/trait_qtl/all_scans.rds"),
    "Real data not available"
  )

  scan_data <- readRDS("/proj/raulab/users/brian/cc_gwas/data/processed/trait_qtl/all_scans.rds")
  threshold_data <- readRDS("/proj/raulab/users/brian/cc_gwas/data/processed/trait_qtl/all_thresholds.rds")

  # Check first scan
  first_scan_name <- names(scan_data)[1]
  cat("Testing with real scan:", first_scan_name, "\n")

  current_scan <- scan_data[[first_scan_name]]

  # Verify all required fields exist
  expect_true("LOD" %in% names(current_scan))
  expect_true("chr" %in% names(current_scan))
  expect_true("pos" %in% names(current_scan))
  expect_true("loci" %in% names(current_scan))
  expect_true("allele.effects" %in% names(current_scan))

  # Check pos conversion works
  pos_integer <- as.integer(current_scan$pos$Mb * 1e6)
  expect_type(pos_integer, "integer")
  expect_false(any(is.na(pos_integer)))

  # Verify we can construct the plotting data frame
  miqtl_df_for_plot <- data.frame(
    marker = names(current_scan$LOD),
    chr    = as.character(current_scan$chr),
    pos    = as.integer(current_scan$pos$Mb * 1e6),
    lod    = current_scan$LOD,
    stringsAsFactors = FALSE
  )

  expect_s3_class(miqtl_df_for_plot, "data.frame")
  expect_gt(nrow(miqtl_df_for_plot), 0)
  expect_type(miqtl_df_for_plot$pos, "integer")

  cat("✓ Real data structure is compatible\n")
})

test_that("can prepare data for plotManhattan with real data", {
  skip_if_not(
    file.exists("/proj/raulab/users/brian/cc_gwas/data/processed/trait_qtl/all_scans.rds"),
    "Real data not available"
  )

  scan_data <- readRDS("/proj/raulab/users/brian/cc_gwas/data/processed/trait_qtl/all_scans.rds")

  # Pick a scan with chromosome 1 data
  first_scan <- scan_data[[1]]

  # Get chromosome 1 data
  current_chr <- "1"
  miqtl_df_for_plot <- data.frame(
    marker = names(first_scan$LOD),
    chr    = as.character(first_scan$chr),
    pos    = as.integer(first_scan$pos$Mb * 1e6),
    lod    = first_scan$LOD,
    stringsAsFactors = FALSE
  )

  # Filter as in packet_core.R line 202-206
  miqtl_df_for_plot <- miqtl_df_for_plot[
    miqtl_df_for_plot$chr == current_chr &
    !is.na(miqtl_df_for_plot$pos) &
    !is.na(miqtl_df_for_plot$lod),
  ]

  # Add required columns for plotManhattan
  miqtl_df_for_plot$chrom <- paste0("chr", miqtl_df_for_plot$chr)
  miqtl_df_for_plot$p <- 10^(-miqtl_df_for_plot$lod)

  # Verify structure is ready for plotManhattan
  expect_true("chrom" %in% names(miqtl_df_for_plot))
  expect_true("pos" %in% names(miqtl_df_for_plot))
  expect_true("p" %in% names(miqtl_df_for_plot))

  expect_type(miqtl_df_for_plot$pos, "integer")
  expect_type(miqtl_df_for_plot$p, "double")

  cat("✓ Data ready for plotManhattan\n")
})

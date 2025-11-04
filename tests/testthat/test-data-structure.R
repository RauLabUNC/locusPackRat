test_that("scan_data has correct structure for plotting", {
  scan_data <- readRDS(test_path("fixtures/test_scan_data.rds"))

  # Check that HR_Ctrl scan exists
  expect_true("HR_Ctrl" %in% names(scan_data))

  current_scan <- scan_data$HR_Ctrl

  # Check all required fields are present
  expect_true("LOD" %in% names(current_scan))
  expect_true("chr" %in% names(current_scan))
  expect_true("pos" %in% names(current_scan))
  expect_true("loci" %in% names(current_scan))  # This was missing!
  expect_true("allele.effects" %in% names(current_scan))

  # Check data types
  expect_type(current_scan$LOD, "double")
  expect_true(is.numeric(current_scan$chr))
  expect_type(current_scan$pos, "list")
  expect_true("Mb" %in% names(current_scan$pos))
  expect_type(current_scan$loci, "character")
  expect_true(is.matrix(current_scan$allele.effects))

  # Check pos conversion to integer works (line 196 of packet_core.R)
  pos_integer <- as.integer(current_scan$pos$Mb * 1e6)
  expect_type(pos_integer, "integer")
  expect_false(any(is.na(pos_integer)))

  # Check loci matches LOD names
  expect_equal(current_scan$loci, names(current_scan$LOD))
})

test_that("threshold_data has correct structure", {
  threshold_data <- readRDS(test_path("fixtures/test_threshold_data.rds"))

  # Check that HR_Ctrl_threshold exists
  expect_true("HR_Ctrl_threshold" %in% names(threshold_data))

  # Check it's a single numeric value
  expect_type(threshold_data$HR_Ctrl_threshold, "double")
  expect_length(threshold_data$HR_Ctrl_threshold, 1)
})

test_that("miqtl_df_for_plot data frame can be constructed", {
  scan_data <- readRDS(test_path("fixtures/test_scan_data.rds"))
  current_scan_data <- scan_data$HR_Ctrl

  # This is the data frame construction from line 193-199 of packet_core.R
  miqtl_df_for_plot <- data.frame(
    marker = names(current_scan_data$LOD),
    chr    = as.character(current_scan_data$chr),
    pos    = as.integer(current_scan_data$pos$Mb * 1e6),
    lod    = current_scan_data$LOD,
    stringsAsFactors = FALSE
  )

  # Verify structure
  expect_s3_class(miqtl_df_for_plot, "data.frame")
  expect_equal(ncol(miqtl_df_for_plot), 4)
  expect_type(miqtl_df_for_plot$pos, "integer")
  expect_false(any(is.na(miqtl_df_for_plot$pos)))

  # Filter for current chromosome (line 202-206)
  current_chr <- "1"
  miqtl_df_filtered <- miqtl_df_for_plot[
    miqtl_df_for_plot$chr == current_chr &
    !is.na(miqtl_df_for_plot$pos) &
    !is.na(miqtl_df_for_plot$lod),
  ]

  expect_gt(nrow(miqtl_df_filtered), 0)
})

#!/usr/bin/env Rscript
# Script to create subset test data from original RDS files
# This creates smaller test data files that can be included in the repository

cat("Creating subset test data from original RDS files...\n")

# Load original data
cat("Loading original scan data...\n")
scan_data <- readRDS("/proj/raulab/users/brian/cc_gwas/data/processed/trait_qtl/all_scans.rds")

cat("Loading original threshold data...\n")
threshold_data <- readRDS("/proj/raulab/users/brian/cc_gwas/data/processed/trait_qtl/all_thresholds.rds")

# Check structure
cat("\nOriginal data structure:\n")
cat("  Scan data: ", length(scan_data), " scans\n", sep="")
cat("  Threshold data: ", length(threshold_data), " thresholds\n", sep="")
cat("  First 5 scan names: ", paste(head(names(scan_data), 5), collapse=", "), "\n", sep="")
cat("  First 5 threshold names: ", paste(head(names(threshold_data), 5), collapse=", "), "\n", sep="")

# Subset to first 3 items
cat("\nSubsetting to first 3 items...\n")
scan_data_subset <- scan_data[1:3]
threshold_data_subset <- threshold_data[1:3]

cat("Subset data structure:\n")
cat("  Scan data: ", length(scan_data_subset), " scans\n", sep="")
cat("  Threshold data: ", length(threshold_data_subset), " thresholds\n", sep="")
cat("  Scan names: ", paste(names(scan_data_subset), collapse=", "), "\n", sep="")
cat("  Threshold names: ", paste(names(threshold_data_subset), collapse=", "), "\n", sep="")

# Save subset data to fixtures directory
output_dir <- "/proj/raulab/users/brian/packrat/tests/testthat/fixtures"

scan_file <- file.path(output_dir, "test_scan_data_subset.rds")
threshold_file <- file.path(output_dir, "test_threshold_data_subset.rds")

cat("\nSaving subset data to:\n")
cat("  Scans: ", scan_file, "\n", sep="")
cat("  Thresholds: ", threshold_file, "\n", sep="")

saveRDS(scan_data_subset, scan_file)
saveRDS(threshold_data_subset, threshold_file)

# Verify files were created and show sizes
if (file.exists(scan_file)) {
  size_mb <- round(file.size(scan_file) / (1024^2), 2)
  cat("✓ Scan data saved: ", size_mb, " MB\n", sep="")
} else {
  cat("✗ Failed to save scan data\n")
}

if (file.exists(threshold_file)) {
  size_mb <- round(file.size(threshold_file) / (1024^2), 2)
  cat("✓ Threshold data saved: ", size_mb, " MB\n", sep="")
} else {
  cat("✗ Failed to save threshold data\n")
}

cat("\nDone! Test data files created successfully.\n")
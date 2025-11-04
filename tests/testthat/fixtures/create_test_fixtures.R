#!/usr/bin/env Rscript
# Script to create test fixtures for locus zoom plot tests

# Set seed for reproducibility
set.seed(42)

# Create scan_data with proper structure
# Must have: LOD (named vector), chr (vector), pos (list with Mb element), allele.effects (matrix)

# Create marker names
marker_names <- paste0("marker_", 1:101)

# Create LOD scores (named vector)
lod_scores <- c(
  # Increase from 95-101 Mb (peak at 101.2)
  seq(3, 7.8, length.out = 62),
  # Decrease from 101-105 Mb
  seq(7.8, 3, length.out = 39)
) + rnorm(101, mean = 0, sd = 0.3)  # Add small noise
names(lod_scores) <- marker_names

# Create positions in Mb
positions_mb <- seq(95, 105, length.out = 101)

# Create allele effects matrix (101 positions x 8 founders)
allele_effects <- matrix(
  rnorm(101 * 8, mean = 0, sd = 0.5),
  ncol = 8
)
colnames(allele_effects) <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")

# Assemble scan_data list
scan_data <- list(
  HR_Ctrl = list(
    LOD = lod_scores,                      # Named vector
    chr = rep(1, 101),                     # Vector
    pos = list(Mb = positions_mb),         # List with Mb element
    loci = marker_names,                   # Character vector of marker names
    allele.effects = allele_effects        # Matrix
  )
)

# Save scan_data
saveRDS(scan_data, "test_scan_data.rds")
cat("Created test_scan_data.rds\n")

# Create threshold_data
# List with significance threshold for HR_Ctrl (with _threshold suffix)
threshold_data <- list(
  HR_Ctrl_threshold = 4.5
)
saveRDS(threshold_data, "test_threshold_data.rds")
cat("Created test_threshold_data.rds\n")

cat("\nFixtures created successfully!\n")

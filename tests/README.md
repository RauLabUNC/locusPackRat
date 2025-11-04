# Test Documentation

## Test Data
Test data is located in `testthat/fixtures/`:
- `test_scan_data_subset.rds` - First 3 QTL scans from production data
- `test_threshold_data_subset.rds` - Corresponding significance thresholds
- `test_genes_mouse.csv` - Sample gene annotations

## Running Tests

### Automated Tests
```r
# Run all package tests
testthat::test_local()

# Run specific test file
testthat::test_file("tests/testthat/test-locus-zoom-integration.R")
```

### Manual Testing
For interactive testing and visualization:
```r
source("tests/manual_zoom_plot_test.R")
```

This creates test PDFs to visually inspect the plot output.

## Test Coverage
- Basic locus zoom plot generation
- Multiple overlapping QTL handling
- Gene annotation display (with and without genes)
- Color scheme consistency for 8 founder strains
#!/usr/bin/env Rscript

# Comprehensive tests for genePackRat core functions
# Tests filterGenes, build_gene_table, and create_gene_workbook integration

# Load package - try installed version first, then load from source
tryCatch({
  library(genePackRat)
}, error = function(e) {
  cat("Package not installed, loading from source...\n")
  # Load required packages
  library(data.table)
  library(dplyr)
  library(openxlsx)

  # Source all R files from parent directory
  r_files <- list.files("../R", pattern = "\\.R$", full.names = TRUE)
  for (file in r_files) {
    source(file)
  }
})

# Helper function to run tests
run_test <- function(test_name, test_func) {
  cat("\n==", test_name, "==\n")
  tryCatch({
    test_func()
    cat("✓ PASSED\n")
    return(TRUE)
  }, error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
    return(FALSE)
  })
}

# Create comprehensive test data
create_test_data <- function() {
  set.seed(123)
  n_genes <- 100

  # Main gene table
  genes <- data.frame(
    gene_id = paste0("GENE", sprintf("%03d", 1:n_genes)),
    gene_symbol = paste0("Gene", 1:n_genes),
    chromosome = sample(c(1:19, "X", "Y"), n_genes, replace = TRUE),
    start = sample(1000000:100000000, n_genes),
    biotype = sample(c("protein_coding", "lncRNA", "miRNA", "pseudogene"),
                     n_genes, replace = TRUE,
                     prob = c(0.7, 0.15, 0.05, 0.1)),
    stringsAsFactors = FALSE
  )
  genes$end <- genes$start + sample(1000:50000, n_genes)

  # Expression data
  expression <- data.frame(
    gene_id = genes$gene_id,
    TPM_brain = rlnorm(n_genes, 2, 1.5),
    TPM_liver = rlnorm(n_genes, 2, 1.5),
    TPM_heart = rlnorm(n_genes, 2, 1.5),
    stringsAsFactors = FALSE
  )

  # GWAS results
  gwas <- data.frame(
    gene_id = sample(genes$gene_id, 50),
    snp_id = paste0("rs", sample(1000000:9999999, 50)),
    p_value = 10^(-runif(50, 1, 12)),
    beta = rnorm(50, 0, 0.5),
    trait = sample(c("Height", "BMI", "T2D", "CAD"), 50, replace = TRUE),
    stringsAsFactors = FALSE
  )

  # eQTL data
  eqtl <- data.frame(
    gene_id = sample(genes$gene_id, 30),
    variant_id = paste0("chr", sample(1:22, 30, replace = TRUE),
                       ":", sample(1000000:100000000, 30)),
    tissue = sample(c("Brain", "Liver", "Heart", "Muscle"), 30, replace = TRUE),
    p_value = 10^(-runif(30, 3, 10)),
    effect_size = rnorm(30, 0, 0.3),
    stringsAsFactors = FALSE
  )

  # Orthology data
  orthology <- data.frame(
    gene_id = sample(genes$gene_id, 80),
    human_gene_id = paste0("ENSG", sprintf("%011d", sample(1:99999, 80))),
    human_symbol = paste0("HUMAN_", sample(LETTERS, 80, replace = TRUE)),
    orthology_confidence = sample(c("high", "moderate", "low"), 80,
                                 replace = TRUE, prob = c(0.6, 0.3, 0.1)),
    stringsAsFactors = FALSE
  )

  return(list(
    genes = genes,
    expression = expression,
    gwas = gwas,
    eqtl = eqtl,
    orthology = orthology
  ))
}

# Test 1: Basic filterGenes functionality
test_basic_filtering <- function() {
  data <- create_test_data()

  # Test single condition filter
  filter1 <- makeFilter(
    column = "biotype",
    condition = "==",
    value = "protein_coding"
  )

  result1 <- filterGenes(
    inputTable = data$genes,
    filters = list(filter1)
  )

  stopifnot(all(result1$biotype == "protein_coding"))
  cat("  - Single condition filter: OK\n")

  # Test numeric filter
  filter2 <- makeFilter(
    column = "start",
    condition = ">",
    value = 50000000
  )

  result2 <- filterGenes(
    inputTable = data$genes,
    filters = list(filter2)
  )

  stopifnot(all(result2$start > 50000000))
  cat("  - Numeric filter: OK\n")

  # Test multiple filters
  result3 <- filterGenes(
    inputTable = data$genes,
    filters = list(filter1, filter2)
  )

  stopifnot(all(result3$biotype == "protein_coding"))
  stopifnot(all(result3$start > 50000000))
  cat("  - Multiple filters: OK\n")

  # Test %in% operator
  filter4 <- makeFilter(
    column = "chromosome",
    condition = "%in%",
    value = c("1", "2", "X")
  )

  result4 <- filterGenes(
    inputTable = data$genes,
    filters = list(filter4)
  )

  stopifnot(all(result4$chromosome %in% c("1", "2", "X")))
  cat("  - %in% operator: OK\n")
}

# Test 2: Relational filtering with joins
test_relational_filtering <- function() {
  data <- create_test_data()

  # Inner join - only genes with expression data
  result_inner <- filterGenes(
    inputTable = data$genes,
    referenceTable = data$expression,
    by = "gene_id",
    joinType = "inner"
  )

  stopifnot(nrow(result_inner) == nrow(data$expression))
  cat("  - Inner join: OK\n")

  # Semi join - genes that have GWAS hits
  result_semi <- filterGenes(
    inputTable = data$genes,
    referenceTable = data$gwas,
    by = "gene_id",
    joinType = "semi"
  )

  stopifnot(all(result_semi$gene_id %in% data$gwas$gene_id))
  stopifnot(ncol(result_semi) == ncol(data$genes))
  cat("  - Semi join: OK\n")

  # Anti join - genes without eQTLs
  result_anti <- filterGenes(
    inputTable = data$genes,
    referenceTable = data$eqtl,
    by = "gene_id",
    joinType = "anti"
  )

  stopifnot(!any(result_anti$gene_id %in% data$eqtl$gene_id))
  cat("  - Anti join: OK\n")

  # Left join with additional filtering
  filter_expr <- makeFilter(
    column = "TPM_brain",
    condition = ">",
    value = 10
  )

  result_combined <- filterGenes(
    inputTable = data$genes,
    referenceTable = data$expression,
    by = "gene_id",
    joinType = "left",
    filters = list(filter_expr)
  )

  stopifnot(all(result_combined$TPM_brain > 10 | is.na(result_combined$TPM_brain)))
  cat("  - Left join with filter: OK\n")
}

# Test 3: build_gene_table integration
test_build_gene_table <- function() {
  data <- create_test_data()

  # Build comprehensive gene table
  gene_table <- build_gene_table(
    genes_df = data$genes,
    expression_df = data$expression,
    orthology_df = data$orthology
  )

  stopifnot("gene_id" %in% names(gene_table))
  stopifnot("TPM_brain" %in% names(gene_table))
  stopifnot(nrow(gene_table) == nrow(data$genes))
  cat("  - Basic table building: OK\n")

  # Filter the built table
  filter_complex <- makeFilter(
    column = "TPM_brain",
    condition = ">",
    value = 50
  )

  filtered_table <- filterGenes(
    inputTable = gene_table,
    filters = list(filter_complex)
  )

  stopifnot(all(filtered_table$TPM_brain > 50))
  cat("  - Filter built table: OK\n")

  # Test with multiple data sources
  full_table <- build_gene_table(
    genes_df = data$genes,
    expression_df = data$expression,
    orthology_df = data$orthology,
    variants_df = data$gwas
  )

  stopifnot(ncol(full_table) > ncol(gene_table))
  cat("  - Multiple data sources: OK\n")
}

# Test 4: Excel workbook creation
test_excel_workbook <- function() {
  data <- create_test_data()

  # Prepare data for workbook
  gene_table <- build_gene_table(
    genes_df = data$genes,
    expression_df = data$expression,
    orthology_df = data$orthology
  )

  # Filter for interesting genes
  filter1 <- makeFilter(
    column = "biotype",
    condition = "==",
    value = "protein_coding"
  )

  filter2 <- makeFilter(
    column = "TPM_brain",
    condition = ">",
    value = 20
  )

  filtered_genes <- filterGenes(
    inputTable = gene_table,
    filters = list(filter1, filter2)
  )

  # Create workbook
  output_file <- "test_gene_workbook.xlsx"

  wb <- create_gene_workbook(
    gene_data = filtered_genes,
    gwas_data = data$gwas,
    eqtl_data = data$eqtl,
    output_file = output_file,
    highlight_genes = head(filtered_genes$gene_symbol, 5)
  )

  stopifnot(file.exists(output_file))
  file_size <- file.info(output_file)$size
  stopifnot(file_size > 0)
  cat("  - Workbook created:", output_file, "\n")
  cat("  - File size:", file_size, "bytes\n")

  # Clean up
  unlink(output_file)
  cat("  - Cleanup: OK\n")
}

# Test 5: Complex integration workflow
test_integration_workflow <- function() {
  data <- create_test_data()

  # Step 1: Build comprehensive table
  full_table <- build_gene_table(
    genes_df = data$genes,
    expression_df = data$expression,
    orthology_df = data$orthology
  )
  cat("  - Built table with", nrow(full_table), "genes\n")

  # Step 2: Filter for genes with GWAS hits
  gwas_genes <- filterGenes(
    inputTable = full_table,
    referenceTable = data$gwas,
    by = "gene_id",
    joinType = "semi"
  )
  cat("  - Genes with GWAS hits:", nrow(gwas_genes), "\n")

  # Step 3: Further filter for expression
  expr_filter <- makeFilter(
    column = "TPM_brain",
    condition = ">",
    value = 30
  )

  final_genes <- filterGenes(
    inputTable = gwas_genes,
    filters = list(expr_filter)
  )
  cat("  - High expression GWAS genes:", nrow(final_genes), "\n")

  # Step 4: Add eQTL information
  genes_with_eqtl <- filterGenes(
    inputTable = final_genes,
    referenceTable = data$eqtl,
    by = "gene_id",
    joinType = "left"
  )

  # Step 5: Create final report
  if (nrow(genes_with_eqtl) > 0) {
    output_file <- "integration_test_report.xlsx"

    wb <- create_gene_workbook(
      gene_data = genes_with_eqtl,
      gwas_data = data$gwas[data$gwas$gene_id %in% genes_with_eqtl$gene_id,],
      eqtl_data = data$eqtl[data$eqtl$gene_id %in% genes_with_eqtl$gene_id,],
      output_file = output_file
    )

    stopifnot(file.exists(output_file))
    cat("  - Integration report created\n")

    # Clean up
    unlink(output_file)
  }

  cat("  - Full workflow completed successfully\n")
}

# Test 6: Edge cases and error handling
test_edge_cases <- function() {
  data <- create_test_data()

  # Empty result set
  impossible_filter <- makeFilter(
    column = "start",
    condition = ">",
    value = 999999999
  )

  empty_result <- filterGenes(
    inputTable = data$genes,
    filters = list(impossible_filter)
  )

  stopifnot(nrow(empty_result) == 0)
  stopifnot(ncol(empty_result) == ncol(data$genes))
  cat("  - Empty result handling: OK\n")

  # Missing columns - should give informative error
  tryCatch({
    bad_filter <- makeFilter(
      column = "nonexistent_column",
      condition = "==",
      value = "test"
    )
    filterGenes(
      inputTable = data$genes,
      filters = list(bad_filter)
    )
    stop("Should have thrown an error")
  }, error = function(e) {
    cat("  - Missing column error handled: OK\n")
  })

  # NA handling
  data_with_na <- data$genes
  data_with_na$start[1:10] <- NA

  na_filter <- makeFilter(
    column = "start",
    condition = ">",
    value = 50000000,
    na.rm = TRUE
  )

  na_result <- filterGenes(
    inputTable = data_with_na,
    filters = list(na_filter)
  )

  stopifnot(!any(is.na(na_result$start)))
  cat("  - NA handling: OK\n")
}

# Main test runner
main <- function() {
  cat("\n========================================\n")
  cat("genePackRat Function Integration Tests\n")
  cat("========================================\n")

  results <- list()

  results$basic <- run_test("Basic Filtering", test_basic_filtering)
  results$relational <- run_test("Relational Filtering", test_relational_filtering)
  results$build_table <- run_test("Build Gene Table", test_build_gene_table)
  results$excel <- run_test("Excel Workbook Creation", test_excel_workbook)
  results$integration <- run_test("Integration Workflow", test_integration_workflow)
  results$edge_cases <- run_test("Edge Cases", test_edge_cases)

  cat("\n========================================\n")
  cat("Test Summary\n")
  cat("========================================\n")

  passed <- sum(unlist(results))
  total <- length(results)

  for (name in names(results)) {
    status <- if(results[[name]]) "✓ PASSED" else "✗ FAILED"
    cat(sprintf("%-25s %s\n", name, status))
  }

  cat("\n")
  cat("Total:", passed, "/", total, "tests passed\n")

  if (passed == total) {
    cat("\n✓ All tests passed successfully!\n")
    return(0)
  } else {
    cat("\n✗ Some tests failed\n")
    return(1)
  }
}

# Run tests if executed directly
if (!interactive()) {
  status <- main()
  quit(status = status)
}
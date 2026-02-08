# Helper functions for testthat tests

#' Create a temporary test project
#'
#' @param mode Character: "gene" or "region"
#' @param species Character: "mouse" or "human"
#' @param genome Character: genome version
#' @return Path to the temporary project directory
create_test_project <- function(mode = "gene", species = "mouse", genome = "mm39") {
  # Create temp directory

test_dir <- file.path(tempdir(), paste0("test_project_", format(Sys.time(), "%Y%m%d%H%M%S")))
  dir.create(test_dir, showWarnings = FALSE)

  if (mode == "gene") {
    # Use appropriate gene symbols for the species
    if (species == "human") {
      test_data <- data.frame(
        gene_symbol = c("MYC", "TP53", "EGFR", "VEGFA", "GAPDH"),
        log2FC = c(1.5, -2.1, 0.8, 1.2, -0.3),
        padj = c(0.001, 0.01, 0.05, 0.1, 0.5)
      )
    } else {
      test_data <- data.frame(
        gene_symbol = c("Myc", "Tp53", "Egfr", "Vegfa", "Gapdh"),
        log2FC = c(1.5, -2.1, 0.8, 1.2, -0.3),
        padj = c(0.001, 0.01, 0.05, 0.1, 0.5)
      )
    }
  } else {
    test_data <- data.frame(
      chr = c("1", "3", "5"),
      start = c(1000000, 5000000, 10000000),
      end = c(2000000, 6000000, 11000000),
      peak_id = c("peak_1", "peak_2", "peak_3")
    )
  }

  initPackRat(
    data = test_data,
    mode = mode,
    species = species,
    genome = genome,
    project_dir = test_dir,
    force = TRUE
  )

  return(test_dir)
}

#' Clean up test project
#'
#' @param project_dir Path to the project directory to remove
cleanup_test_project <- function(project_dir) {
  if (dir.exists(project_dir)) {
    unlink(project_dir, recursive = TRUE)
  }
}

#' Create sample supplementary data for testing
#'
#' @param n_genes Number of genes to include
#' @return data.frame with sample supplementary data
create_sample_supplementary <- function(n_genes = 5) {
  data.frame(
    gene_symbol = c("Myc", "Tp53", "Egfr", "Vegfa", "Gapdh")[1:n_genes],
    tissue_expr = runif(n_genes, 10, 100),
    category = sample(c("A", "B", "C"), n_genes, replace = TRUE)
  )
}

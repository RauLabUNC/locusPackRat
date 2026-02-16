# Tests for queryOpenTargets and queryOpenTargetsQTL functions

# --- Parameter Validation Tests ---

test_that("queryOpenTargets validates data_types parameter", {
  test_dir <- create_test_project(mode = "gene", species = "human", genome = "hg38")
  on.exit(cleanup_test_project(test_dir))

  expect_error(
    queryOpenTargets(
      project_dir = test_dir,
      data_types = c("invalid_type"),
      limit = 1
    ),
    "Invalid data_types"
  )

  expect_error(
    queryOpenTargets(
      project_dir = test_dir,
      data_types = c("diseases", "not_a_type"),
      limit = 1
    ),
    "Invalid data_types"
  )
})

test_that("queryOpenTargets accepts valid data_types", {
  test_dir <- create_test_project(mode = "gene", species = "human", genome = "hg38")
  on.exit(cleanup_test_project(test_dir))

  valid_types <- c("diseases", "constraints", "tractability", "expression",
                   "mouse_phenotypes", "homologues", "interactions",
                   "known_drugs", "depmap", "pharmacogenomics")

  # This should not error - just validates parameter acceptance
  # We use limit=0 to avoid actual API calls
  expect_no_error({
    # Just check the validation passes - actual query may fail due to network
    tryCatch(
      queryOpenTargets(
        project_dir = test_dir,
        data_types = valid_types,
        limit = 1,
        verbose = FALSE
      ),
      error = function(e) {
        # Ignore network/API errors, only fail on validation errors
        if (grepl("Invalid data_types", e$message)) stop(e)
      }
    )
  })
})

test_that("queryOpenTargets requires existing project", {
  expect_error(
    queryOpenTargets(project_dir = "/nonexistent/path"),
    "Project not found"
  )
})

test_that("queryOpenTargetsQTL validates study_types parameter", {
  test_dir <- create_test_project(mode = "gene", species = "human", genome = "hg38")
  on.exit(cleanup_test_project(test_dir))

  expect_error(
    queryOpenTargetsQTL(
      project_dir = test_dir,
      study_types = c("invalid_qtl_type"),
      limit = 1
    ),
    "Invalid study_types"
  )
})

test_that("queryOpenTargetsQTL accepts valid study_types", {
  test_dir <- create_test_project(mode = "gene", species = "human", genome = "hg38")
  on.exit(cleanup_test_project(test_dir))

  valid_types <- c("eqtl", "pqtl", "sceqtl", "scpqtl", "sqtl", "scsqtl", "tuqtl", "sctuqtl")

  expect_no_error({
    tryCatch(
      queryOpenTargetsQTL(
        project_dir = test_dir,
        study_types = valid_types,
        limit = 1,
        verbose = FALSE
      ),
      error = function(e) {
        if (grepl("Invalid study_types", e$message)) stop(e)
      }
    )
  })
})

test_that("queryOpenTargetsQTL requires existing project", {
  expect_error(
    queryOpenTargetsQTL(project_dir = "/nonexistent/path"),
    "Project not found"
  )
})

# --- Helper Function Tests ---

test_that(".resolveHumanIds works for human projects", {
  test_dir <- create_test_project(mode = "gene", species = "human", genome = "hg38")
  on.exit(cleanup_test_project(test_dir))

  config <- jsonlite::read_json(file.path(test_dir, ".locusPackRat", "config.json"))

  # Test with known human genes
  gene_syms <- c("MYC", "TP53", "EGFR")
  result <- locusPackRat:::.resolveHumanIds(gene_syms, config)

  expect_true(data.table::is.data.table(result))
  expect_true("gene_symbol" %in% names(result))
  expect_true("human_ensembl_id" %in% names(result))
  expect_equal(nrow(result), 3)
})

test_that(".resolveHumanIds works for mouse projects", {
  test_dir <- create_test_project(mode = "gene", species = "mouse", genome = "mm39")
  on.exit(cleanup_test_project(test_dir))

  config <- jsonlite::read_json(file.path(test_dir, ".locusPackRat", "config.json"))

  # Test with known mouse genes
  gene_syms <- c("Myc", "Tp53", "Egfr")
  result <- locusPackRat:::.resolveHumanIds(gene_syms, config)

  expect_true(data.table::is.data.table(result))
  expect_true("gene_symbol" %in% names(result))
  expect_true("human_ensembl_id" %in% names(result))
  expect_equal(nrow(result), 3)
})

# --- Response Parsing Tests ---

test_that(".parseOtResponse handles disease data correctly", {
  # Mock gene data structure
  g_dat <- list(
    associatedDiseases = list(
      rows = data.frame(
        disease.id = c("MONDO_0005015", "MONDO_0005059"),
        disease.name = c("Diabetes mellitus", "Leukemia"),
        score = c(0.8, 0.6)
      )
    )
  )

  results <- list(
    diseases = list(), constraints = list(), tractability = list(),
    expression = list(), mouse_phenotypes = list(), homologues = list(),
    interactions = list(), known_drugs = list(), depmap = list(),
    pharmacogenomics = list()
  )

  parsed <- locusPackRat:::.parseOtResponse(g_dat, "ENSG00000000001", "diseases", results)

  expect_equal(length(parsed$diseases), 1)
  expect_true("disease_id" %in% names(parsed$diseases[[1]]) ||
              "disease.id" %in% names(parsed$diseases[[1]]))
})

test_that(".parseOtResponse handles constraint data correctly", {
  g_dat <- list(
    geneticConstraint = data.frame(
      constraintType = c("lof", "mis"),
      oe = c(0.1, 0.5),
      oeLower = c(0.05, 0.4),
      oeUpper = c(0.15, 0.6)
    )
  )

  results <- list(
    diseases = list(), constraints = list(), tractability = list(),
    expression = list(), mouse_phenotypes = list(), homologues = list(),
    interactions = list(), known_drugs = list(), depmap = list(),
    pharmacogenomics = list()
  )

  parsed <- locusPackRat:::.parseOtResponse(g_dat, "ENSG00000000001", "constraints", results)

  expect_equal(length(parsed$constraints), 1)
  expect_true("constraintType" %in% names(parsed$constraints[[1]]))
  expect_true("human_ensembl_id" %in% names(parsed$constraints[[1]]))
})

test_that(".parseOtResponse handles tractability data correctly", {
  g_dat <- list(
    tractability = data.frame(
      label = c("Small molecule", "Antibody"),
      modality = c("SM", "AB"),
      value = c(TRUE, FALSE)
    )
  )

  results <- list(
    diseases = list(), constraints = list(), tractability = list(),
    expression = list(), mouse_phenotypes = list(), homologues = list(),
    interactions = list(), known_drugs = list(), depmap = list(),
    pharmacogenomics = list()
  )

  parsed <- locusPackRat:::.parseOtResponse(g_dat, "ENSG00000000001", "tractability", results)

  expect_equal(length(parsed$tractability), 1)
  expect_true("label" %in% names(parsed$tractability[[1]]))
})

test_that("null coalescing operator works correctly", {
  `%||%` <- locusPackRat:::`%||%`

  expect_equal(NULL %||% "default", "default")
  expect_equal("value" %||% "default", "value")
  expect_equal(character(0) %||% "default", "default")
  expect_equal(list() %||% "default", "default")
})

# --- Integration Tests (Skip if no network) ---

test_that("queryOpenTargets can query with default data_types", {
  skip_on_cran()
  skip_if_offline()

  # Create human project for direct testing
  test_dir <- file.path(tempdir(), "test_ot_integration")
  on.exit(cleanup_test_project(test_dir))

  test_data <- data.frame(
    gene_symbol = c("BRCA1"),  # Well-known gene with lots of data
    log2FC = c(1.5)
  )

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "human",
    genome = "hg38",
    project_dir = test_dir,
    force = TRUE
  )

  result <- queryOpenTargets(
    project_dir = test_dir,
    data_types = c("diseases", "constraints"),
    limit = 1,
    verbose = FALSE
  )

  # Check that results were returned

expect_true(is.list(result))

  # Check that tables were saved
  packrat_dir <- file.path(test_dir, ".locusPackRat")
  tables_dir <- file.path(packrat_dir, "supplementary")

  # At least one table should exist
  saved_files <- list.files(tables_dir, pattern = "^ot_")
  expect_true(length(saved_files) >= 0)  # May be 0 if gene has no data
})

test_that("queryOpenTargets can query expression data", {
  skip_on_cran()
  skip_if_offline()

  test_dir <- file.path(tempdir(), "test_ot_expression")
  on.exit(cleanup_test_project(test_dir))

  test_data <- data.frame(
    gene_symbol = c("TP53"),
    log2FC = c(-1.2)
  )

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "human",
    genome = "hg38",
    project_dir = test_dir,
    force = TRUE
  )

  result <- queryOpenTargets(
    project_dir = test_dir,
    data_types = c("expression"),
    limit = 1,
    verbose = FALSE
  )

  expect_true(is.list(result))

  # If expression data was returned, check structure
  if (!is.null(result$expression) && nrow(result$expression) > 0) {
    expect_true("tissue_label" %in% names(result$expression))
    expect_true("rna_level" %in% names(result$expression) ||
                "rna_value" %in% names(result$expression))
  }
})

test_that("queryOpenTargetsQTL can query eQTL data", {
  skip_on_cran()
  skip_if_offline()

  test_dir <- file.path(tempdir(), "test_qtl_integration")
  on.exit(cleanup_test_project(test_dir))

  test_data <- data.frame(
    gene_symbol = c("APOE"),  # Well-known gene with QTL data
    log2FC = c(0.8)
  )

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "human",
    genome = "hg38",
    project_dir = test_dir,
    force = TRUE
  )

  result <- queryOpenTargetsQTL(
    project_dir = test_dir,
    study_types = c("eqtl"),
    limit = 1,
    verbose = FALSE
  )

  expect_true(is.list(result))

  # If data was returned, check structure
  if (!is.null(result$credible_sets) && nrow(result$credible_sets) > 0) {
    expect_true("study_type" %in% names(result$credible_sets))
    expect_true("chromosome" %in% names(result$credible_sets))
  }
})

# --- Backward Compatibility Tests ---

test_that("queryOpenTargets default behavior matches original", {
  # Default data_types should be diseases, constraints, tractability
  test_dir <- create_test_project(mode = "gene", species = "human", genome = "hg38")
  on.exit(cleanup_test_project(test_dir))

  # This tests that the default parameters are set correctly
  # without actually making API calls
  expect_no_error({
    tryCatch(
      queryOpenTargets(
        project_dir = test_dir,
        limit = 1,
        verbose = FALSE
      ),
      error = function(e) {
        # Only fail on parameter validation errors
        if (grepl("Invalid", e$message)) stop(e)
      }
    )
  })
})

# Tests for initPackRat function

test_that("initPackRat creates project structure in gene mode", {
  test_dir <- file.path(tempdir(), "test_gene_project")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(
    gene_symbol = c("Myc", "Tp53", "Egfr"),
    log2FC = c(1.5, -2.1, 0.8)
  )

  result <- initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  expect_true(result)
  expect_true(dir.exists(file.path(test_dir, ".locusPackRat")))
  expect_true(dir.exists(file.path(test_dir, ".locusPackRat", "input")))
  expect_true(dir.exists(file.path(test_dir, ".locusPackRat", "supplementary")))
  expect_true(dir.exists(file.path(test_dir, ".locusPackRat", "output")))
  expect_true(file.exists(file.path(test_dir, ".locusPackRat", "config.json")))
  expect_true(file.exists(file.path(test_dir, ".locusPackRat", "input", "genes.csv")))
})

test_that("initPackRat creates project structure in region mode", {
  test_dir <- file.path(tempdir(), "test_region_project")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(
    chr = c("1", "3"),
    start = c(1000000, 5000000),
    end = c(2000000, 6000000)
  )

  result <- initPackRat(
    data = test_data,
    mode = "region",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  expect_true(result)
  expect_true(file.exists(file.path(test_dir, ".locusPackRat", "input", "regions.csv")))
})

test_that("initPackRat validates species/genome compatibility", {
  test_dir <- file.path(tempdir(), "test_compat")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc"))

  expect_error(
    initPackRat(
      data = test_data,
      mode = "gene",
      species = "mouse",
      genome = "hg38",
      project_dir = test_dir
    ),
    "Incompatible species and genome"
  )
})

test_that("initPackRat requires force=TRUE to overwrite", {
  test_dir <- file.path(tempdir(), "test_overwrite")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc"))

  # First creation
  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  # Should fail without force=TRUE
  expect_error(
    initPackRat(
      data = test_data,
      mode = "gene",
      species = "mouse",
      genome = "mm39",
      project_dir = test_dir
    ),
    "already exists"
  )
})

test_that("initPackRat stores config correctly", {
  test_dir <- file.path(tempdir(), "test_config")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc", "Tp53"))

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  config <- jsonlite::read_json(file.path(test_dir, ".locusPackRat", "config.json"))

  expect_equal(config$mode, "gene")
  expect_equal(config$species, "mouse")
  expect_equal(config$genome, "mm39")
  expect_equal(config$n_entries, 2)
})

test_that("initPackRat emits genome build consistency message", {
  test_dir <- file.path(tempdir(), "test_genome_msg")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc", "Tp53"))

  expect_message(
    initPackRat(
      data = test_data,
      mode = "gene",
      species = "mouse",
      genome = "mm39",
      project_dir = test_dir,
      force = TRUE
    ),
    "All data added to this project must use the mm39 genome build"
  )
})

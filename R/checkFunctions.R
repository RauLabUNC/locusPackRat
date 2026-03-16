#' Check locusPackRat Project Status
#'
#' Diagnoses a locusPackRat project: validates the project directory, prints
#' project metadata, checks file integrity, and reports optional dependency
#' availability. Useful for troubleshooting or confirming a project is ready
#' for analysis.
#'
#' @param project_dir Character: Path to project directory containing .locusPackRat
#' @param verbose Logical: If TRUE (default), check optional dependencies and
#'   report output staleness
#'
#' @return Invisibly returns a list with:
#'   \item{valid}{Logical: TRUE if no critical issues found}
#'   \item{project_dir}{Normalized project directory path}
#'   \item{config}{Parsed config.json contents}
#'   \item{tables}{data.table from listPackRatTables()}
#'   \item{issues}{Character vector of warnings/problems found}
#'   \item{suggestions}{Character vector of suggested next steps}
#'
#' @examples
#' \dontrun{
#' # Check current directory
#' lpr_check()
#'
#' # Check a specific project
#' lpr_check("my_analysis")
#'
#' # Quick check without dependency info
#' lpr_check("my_analysis", verbose = FALSE)
#' }
#'
#' @importFrom jsonlite read_json
#' @importFrom data.table fread
#'
#' @export
lpr_check <- function(project_dir = ".", verbose = TRUE) {

  issues <- character(0)
  suggestions <- character(0)
  config <- list()
  tables <- data.table::data.table()

  # --- Project existence ---
  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) {
    message("-- locusPackRat Project Check --")
    message(sprintf("  Directory: %s", normalizePath(project_dir, mustWork = FALSE)))
    message("  [FAIL] No .locusPackRat directory found")
    message("  Run initPackRat() to create a project here.")
    return(invisible(list(
      valid = FALSE,
      project_dir = normalizePath(project_dir, mustWork = FALSE),
      config = config,
      tables = tables,
      issues = "No .locusPackRat directory found",
      suggestions = "Run initPackRat() to create a project"
    )))
  }

  # --- Config integrity ---
  config_file <- file.path(packrat_dir, "config.json")
  if (!file.exists(config_file)) {
    issues <- c(issues, "config.json not found - project may be corrupted")
  } else {
    config <- tryCatch(
      jsonlite::read_json(config_file),
      error = function(e) {
        issues <<- c(issues, sprintf("config.json is unreadable: %s", e$message))
        list()
      }
    )
  }

  # --- Header ---
  message("-- locusPackRat Project Check --")
  message(sprintf("  Directory: %s", normalizePath(project_dir, mustWork = FALSE)))

  if (length(config) > 0) {
    message(sprintf("  Mode: %s | Species: %s | Genome: %s",
                    config$mode %||% "?",
                    config$species %||% "?",
                    config$genome %||% "?"))
    message(sprintf("  Initialized: %s | Package: v%s",
                    config$initialization_date %||% "?",
                    config$package_version %||% "?"))

    # Version check
    current_ver <- tryCatch(
      as.character(utils::packageVersion("locusPackRat")),
      error = function(e) NULL
    )
    if (!is.null(current_ver) && !is.null(config$package_version) &&
        config$package_version != current_ver) {
      issues <- c(issues, sprintf(
        "Project created with v%s, current package is v%s",
        config$package_version, current_ver
      ))
    }
  }

  # --- Input files ---
  input_dir <- file.path(packrat_dir, "input")
  if (dir.exists(input_dir)) {
    genes_file <- file.path(input_dir, "genes.csv")
    regions_file <- file.path(input_dir, "regions.csv")

    if (file.exists(genes_file)) {
      n_genes <- nrow(data.table::fread(genes_file, select = 1L))
      message(sprintf("  Entries: %d genes", n_genes))
    } else if (file.exists(regions_file)) {
      regions_dt <- data.table::fread(regions_file)
      n_regions <- nrow(regions_dt)
      # Count genes within regions if available
      if ("genes" %in% names(regions_dt)) {
        all_genes <- unlist(strsplit(regions_dt$genes, ", "))
        message(sprintf("  Entries: %d regions, %d genes",
                        n_regions, length(unique(all_genes))))
      } else {
        message(sprintf("  Entries: %d regions", n_regions))
      }
    } else {
      issues <- c(issues, "No genes.csv or regions.csv found in input/")
    }
  } else {
    issues <- c(issues, "input/ directory missing")
  }

  # --- Supplementary tables ---
  supp_dir <- file.path(packrat_dir, "supplementary")
  supp_config_names <- if (!is.null(config$supplementary_tables)) {
    names(config$supplementary_tables)
  } else {
    character(0)
  }

  supp_files_on_disk <- character(0)
  if (dir.exists(supp_dir)) {
    supp_files_on_disk <- sub("\\.csv$", "",
                              list.files(supp_dir, pattern = "\\.csv$"))
  }

  if (length(supp_files_on_disk) > 0) {
    # Use listPackRatTables for the table summary (suppress its messages)
    tables <- tryCatch(
      suppressMessages(listPackRatTables(project_dir = project_dir)),
      error = function(e) data.table::data.table()
    )

    message(sprintf("\n  Supplementary tables (%d):", length(supp_files_on_disk)))
    for (sf in supp_files_on_disk) {
      row <- if (nrow(tables) > 0 && sf %in% tables$table_name) {
        tables[tables$table_name == sf, ]
      } else {
        NULL
      }
      if (!is.null(row) && nrow(row) == 1) {
        message(sprintf("    %-25s | %s rows | linked by %s",
                        sf,
                        format(row$n_rows, big.mark = ","),
                        row$link_by %||% "?"))
      } else {
        message(sprintf("    %-25s | (not in config)", sf))
      }
    }

    # Check consistency: config vs disk
    orphans <- setdiff(supp_files_on_disk, supp_config_names)
    missing <- setdiff(supp_config_names, supp_files_on_disk)
    if (length(orphans) > 0) {
      issues <- c(issues, sprintf(
        "Orphan files in supplementary/ (not in config): %s",
        paste(orphans, collapse = ", ")
      ))
    }
    if (length(missing) > 0) {
      issues <- c(issues, sprintf(
        "Config references missing files: %s",
        paste(missing, collapse = ", ")
      ))
    }
  } else {
    message("\n  Supplementary tables: none")
    suggestions <- c(suggestions, "Add data with addRatTable() or lpr_add_table()")
  }

  # --- Output files ---
  output_dir <- file.path(packrat_dir, "output")
  if (dir.exists(output_dir)) {
    output_files <- list.files(output_dir, full.names = TRUE)
    if (length(output_files) > 0) {
      message(sprintf("\n  Output files (%d):", length(output_files)))
      for (of in output_files) {
        finfo <- file.info(of)
        message(sprintf("    %s (%s, %s)",
                        basename(of),
                        .formatFileSize(finfo$size),
                        format(finfo$mtime, "%Y-%m-%d")))
      }

      # Check staleness
      if (verbose && length(supp_files_on_disk) > 0) {
        supp_paths <- file.path(supp_dir,
                                paste0(supp_files_on_disk, ".csv"))
        newest_supp <- max(file.info(supp_paths)$mtime, na.rm = TRUE)
        oldest_output <- min(file.info(output_files)$mtime, na.rm = TRUE)
        if (oldest_output < newest_supp) {
          issues <- c(issues,
                      "Output is older than supplementary data - re-run makeGeneSheet() to update")
        }
      }
    } else {
      message("\n  Output files: none")
      suggestions <- c(suggestions,
                       "Run makeGeneSheet() or lpr_export() to generate output")
    }
  }

  # --- Optional dependency checks ---
  if (verbose) {
    message("\n  Dependency checks:")
    .check_dep <- function(pkg, purpose) {
      if (requireNamespace(pkg, quietly = TRUE)) {
        message(sprintf("    [ok]   %s", pkg))
      } else {
        message(sprintf("    [skip] %s (needed for %s)", pkg, purpose))
      }
    }
    .check_dep("plotgardener", "locus zoom plots")
    .check_dep("httr", "API queries")

    # Check genome-specific TxDb
    if (!is.null(config$genome) && !is.null(config$species)) {
      txdb_map <- list(
        mm39 = "TxDb.Mmusculus.UCSC.mm39.knownGene",
        mm10 = "TxDb.Mmusculus.UCSC.mm10.knownGene",
        hg38 = "TxDb.Hsapiens.UCSC.hg38.knownGene",
        hg19 = "TxDb.Hsapiens.UCSC.hg19.knownGene"
      )
      txdb <- txdb_map[[config$genome]]
      if (!is.null(txdb)) {
        .check_dep(txdb, "locus zoom gene track")
      }
    }

    # Check zip availability
    has_zip_pkg <- requireNamespace("zip", quietly = TRUE)
    has_sys_zip <- nzchar(Sys.which("zip")) || nzchar(Sys.getenv("R_ZIPCMD"))
    if (has_zip_pkg || has_sys_zip) {
      message("    [ok]   zip (for buildPacket)")
    } else {
      message("    [skip] zip (needed for buildPacket)")
    }
  }

  # --- Summary ---
  if (length(issues) > 0) {
    message("\n  Issues:")
    for (iss in issues) {
      message(sprintf("    [!] %s", iss))
    }
  }
  if (length(suggestions) > 0) {
    message("\n  Suggestions:")
    for (sug in suggestions) {
      message(sprintf("    -> %s", sug))
    }
  }
  if (length(issues) == 0 && length(suggestions) == 0) {
    message("\n  All checks passed.")
  }

  invisible(list(
    valid = length(issues) == 0,
    project_dir = normalizePath(project_dir, mustWork = FALSE),
    config = config,
    tables = tables,
    issues = issues,
    suggestions = suggestions
  ))
}

#' @rdname lpr_check
#' @export
checkPackRat <- lpr_check

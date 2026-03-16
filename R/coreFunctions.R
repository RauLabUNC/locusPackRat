#' Initialize locusPackRat Project
#'
#' Creates a .locusPackRat directory structure and initializes project with gene or region data
#'
#' @param data data.frame or data.table with gene/region information
#' @param mode Character: "gene" or "region" mode
#' @param species Character: "human" or "mouse"
#' @param genome Character: "hg19", "hg38", "mm10", or "mm39"
#' @param project_dir Character: Path to project directory (default: current directory)
#' @param keep_pseudo Logical: If FALSE, removes pseudogenes (Gm* and *Rik genes)
#'   from region list. Does nothing in gene mode.
#' @param force Logical: Overwrite existing .locusPackRat directory if it exists
#' @param overlap_mode Character: "any" (default) keeps genes with any overlap,
#'                      "complete" keeps only genes entirely within regions
#' @return Invisible TRUE on success
#'
#' @note Genome coordinates and orthology mappings are loaded from cached data files
#'       in tests/test_data/ (or inst/extdata/ when installed as a package).
#'       Currently supports mm39 (mouse) and hg38 (human) genomes.
#'
#' @examples
#' \dontrun{
#' # Gene mode: Initialize with a gene list
#' genes <- data.frame(
#'   gene_symbol = c("Myc", "Tp53", "Egfr"),
#'   log2FC = c(1.5, -2.1, 0.8)
#' )
#' initPackRat(genes, mode = "gene", species = "mouse", genome = "mm39",
#'             project_dir = "my_analysis")
#'
#' # Region mode: Initialize with genomic intervals
#' regions <- data.frame(
#'   chr = c("1", "3"),
#'   start = c(1000000, 5000000),
#'   end = c(2000000, 6000000)
#' )
#' initPackRat(regions, mode = "region", species = "mouse", genome = "mm39",
#'             project_dir = "qtl_analysis")
#' }
#'
#' @importFrom utils head packageVersion tail
#' @importFrom data.table setnames fread fwrite setDT := copy
#' @importFrom dplyr bind_rows
#' @importFrom jsonlite write_json read_json
#'
#' @export
initPackRat <- function(data,
                       mode = c("gene", "region"),
                       species = c("human", "mouse"),
                       genome = c("hg38", "hg19", "mm39", "mm10"),
                       project_dir = ".",
                       keep_pseudo = FALSE,
                       overlap_mode = c("any", "complete"),
                       force = FALSE) {

  # Validate inputs
  mode <- match.arg(mode)
  species <- match.arg(species)
  genome <- match.arg(genome)
  overlap_mode <- match.arg(overlap_mode)

  # Check species/genome compatibility
  if ((species == "human" && genome %in% c("mm10", "mm39")) ||
      (species == "mouse" && genome %in% c("hg19", "hg38"))) {
    stop("Incompatible species and genome specified")
  }

  # Set up directory structure
  packrat_dir <- file.path(project_dir, ".locusPackRat")

  if (dir.exists(packrat_dir) && !force) {
    stop(".locusPackRat directory already exists. Use force=TRUE to overwrite.")
  }

  message("Initializing locusPackRat project...")
  message(sprintf("Note: All data added to this project must use the %s genome build. ",
                  genome),
          "Mixing genome builds will cause incorrect coordinate matching. ",
          "Use rtracklayer::liftOver() to convert coordinates if needed.")

  # Create directory structure
  dirs <- c(
    packrat_dir,
    file.path(packrat_dir, "input"),
    file.path(packrat_dir, "supplementary"),
    file.path(packrat_dir, "cache"),
    file.path(packrat_dir, "output")
  )

  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
    }
  }

  # Convert input to data.table
  setDT(data)

  # Process based on mode
  if (mode == "gene") {
    message("Processing gene list...")
    processed_data <- .processGeneInput(data, species, genome)
    output_file <- file.path(packrat_dir, "input", "genes.csv")
  } else {
    message("Processing region list...")
    processed_data <- .processRegionInput(data, species, genome, keep_pseudo, overlap_mode)
    output_file <- file.path(packrat_dir, "input", "regions.csv")
  }

  # Generate orthology table
  message("Generating orthology information...")
  orthology_data <- .generateOrthologyTable(processed_data, species, genome)

  # Save processed data
  fwrite(processed_data, output_file)
  message(sprintf("Saved %s data to %s", mode, output_file))

  # Save orthology data
  if (!is.null(orthology_data) && nrow(orthology_data) > 0) {
    ortho_file <- file.path(packrat_dir, "input", "orthology.csv")
    fwrite(orthology_data, ortho_file)
    message(sprintf("Saved orthology data to %s", ortho_file))
  }

  # Create config file
  # Get package version safely (may not be installed as package yet)
  pkg_version <- tryCatch(
    packageVersion("locusPackRat"),
    error = function(e) "0.1.0-dev"
  )

  config <- list(
    mode = mode,
    species = species,
    genome = genome,
    overlap_mode = if (mode == "region") overlap_mode else NULL,
    initialization_date = Sys.Date(),
    package_version = as.character(pkg_version),
    input_columns = names(processed_data),
    n_entries = nrow(processed_data)
  )

  config_file <- file.path(packrat_dir, "config.json")
  jsonlite::write_json(config, config_file, pretty = TRUE, auto_unbox = TRUE)
  message(sprintf("Created config file: %s", config_file))

  message("\nlocusPackRat project initialized successfully!")
  message(sprintf("Mode: %s | Species: %s | Genome: %s", mode, species, genome))
  message(sprintf("Processed %d %s", nrow(processed_data), ifelse(mode == "gene", "genes", "regions")))

  invisible(TRUE)
} 

#' Add Supplementary Table to locusPackRat Project
#'
#' Adds supplementary data tables linked to genes or regions in an existing project.
#' Can handle gene-based data, genomic regions, or point coordinates (e.g., SNPs).
#'
#' @param data data.frame or data.table with supplementary information
#' @param table_name Character: Name for the output table (without .csv extension)
#' @param link_type Character: "gene", "region", or "point".
#'                  "point" treats data as genomic positions (SNPs) and maps them to overlapping regions.
#' @param link_by Character: Column name in data to use for linking.
#' @param abbreviation Character: Optional short prefix (e.g., "mm", "otd") used to
#'                      distinguish columns from this table when exporting to Excel.
#'                      If another table has the same column name, this prefix will
#'                      be prepended (e.g., "mm_description"). Must be unique per project.
#' @param project_dir Character: Path to project directory containing .locusPackRat
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' # Link by gene symbol
#' expression_data <- data.frame(
#'   gene_symbol = c("Myc", "Tp53"),
#'   tissue_expr = c(100, 50)
#' )
#' addRatTable(expression_data, table_name = "expression",
#'             link_type = "gene", link_by = "gene_symbol",
#'             project_dir = "my_analysis")
#'
#' # Link by genomic region
#' peak_data <- data.frame(
#'   chr = c("1", "1"), start = c(1000, 2000), end = c(1500, 2500),
#'   score = c(10, 20)
#' )
#' addRatTable(peak_data, table_name = "peaks", abbreviation = "pk",
#'             link_type = "region", link_by = "chr,start,end",
#'             project_dir = "my_analysis")
#' }
#'
#' @importFrom data.table fread fwrite setDT merge.data.table := setnames
#' @importFrom jsonlite read_json write_json
#'
#' @export
addRatTable <- function(data,
                        table_name,
                        link_type = c("gene", "region", "point"),
                        link_by = NULL,
                        abbreviation = NULL,
                        project_dir = ".") {

  # Validate inputs
  link_type <- match.arg(link_type)
  # Preserve the original type for the config file
  original_link_type <- link_type
  
  if (missing(table_name) || is.null(table_name)) {
    stop("table_name must be provided")
  }

  # Validate table_name characters (alphanumeric, underscore, hyphen, space only)
  if (!grepl("^[a-zA-Z0-9_ -]+$", table_name)) {
    stop("table_name must contain only letters, numbers, underscores, hyphens, or spaces")
  }

  # Check for existing project
  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) {
    stop("No .locusPackRat directory found. Run initPackRat() first.")
  }

  # Read config
  config_file <- file.path(packrat_dir, "config.json")
  if (!file.exists(config_file)) {
    stop("Config file not found. Project may be corrupted.")
  }
  config <- jsonlite::read_json(config_file)

  # Validate abbreviation uniqueness
  if (!is.null(abbreviation)) {
    existing_abbrevs <- sapply(config$supplementary_tables, function(t) t$abbreviation)
    existing_abbrevs <- existing_abbrevs[!sapply(existing_abbrevs, is.null)]
    if (abbreviation %in% existing_abbrevs) {
      stop(sprintf("Abbreviation '%s' already used by another table.", abbreviation))
    }
  }

  message(sprintf("Adding supplementary table to %s %s project...",
                  config$species, config$genome))

  # Convert input to data.table
  setDT(data)

  # --- HANDLE POINT DATA ---
  # If "point", normalize to "region" format (start == end)
  if (link_type == "point") {
    message("Normalizing point data to genomic intervals...")
    
    # 1. Standardize Chromosome column
    if (!"chr" %in% names(data)) {
       # Try standard alternatives (e.g., chrom, chromosome, seqname)
       chr_col <- grep("^chr|^chrom|^seqname", names(data), value=TRUE, ignore.case=TRUE)[1]
       if (!is.na(chr_col)) {
         setnames(data, chr_col, "chr")
         message(sprintf("  Mapping '%s' to 'chr'", chr_col))
       } else {
         stop("Link type 'point' requires a chromosome column (e.g., 'chr', 'chrom')")
       }
    }
    
    # 2. Standardize Position column
    if (!"pos" %in% names(data)) {
       # Try standard alternatives (e.g., pos, bp, location, start)
       pos_col <- grep("^pos|^bp|^loc|^start", names(data), value=TRUE, ignore.case=TRUE)[1]
       if (!is.na(pos_col)) {
         setnames(data, pos_col, "pos")
         message(sprintf("  Mapping '%s' to 'pos'", pos_col))
       } else {
         stop("Link type 'point' requires a position column (e.g., 'pos', 'bp')")
       }
    }
    
    # 3. Create start/end for overlap joining (Point is an interval of width 0 or 1)
    data[, start := as.numeric(pos)]
    data[, end := as.numeric(pos)]
    
    # Switch processing mode to "region" for the rest of the function
    link_type <- "region"
  }

  # Load existing data based on link type and config mode
  if (link_type == "gene") {
    # Check if we have genes to link to
    genes_file <- file.path(packrat_dir, "input", "genes.csv")
    if (config$mode == "region") {
      # For region mode, we need to extract genes from regions first
      regions_file <- file.path(packrat_dir, "input", "regions.csv")
      if (!file.exists(regions_file)) {
        stop("Regions file not found, cannot link gene data")
      }
      existing_data <- fread(regions_file)
      # Get unique genes from regions
      if ("genes" %in% names(existing_data)) {
        # Assuming genes are stored as comma-separated in a column
        all_genes <- unlist(strsplit(existing_data$genes, ","))
        all_genes <- unique(trimws(all_genes))
        existing_data <- data.table(gene_symbol = all_genes)
      } else {
        stop("No gene information found in regions data")
      }
    } else {
      # Gene mode - direct load
      if (!file.exists(genes_file)) {
        stop("Genes file not found")
      }
      existing_data <- fread(genes_file)
    }

    # Perform linking
    message(sprintf("Linking data by %s...", ifelse(is.null(link_by), "auto-detection", link_by)))
    linked_data <- .linkGeneData(data, existing_data, link_by)

  } else {
    linked_data <- data
  }

  # Save supplementary table
  output_file <- file.path(packrat_dir, "supplementary", paste0(table_name, ".csv"))
  fwrite(linked_data, output_file)
  message(sprintf("Saved supplementary table to %s", output_file))
  message(sprintf("Linked %d of %d input rows", nrow(linked_data), nrow(data)))

  # Update config with table info
  if (is.null(config$supplementary_tables)) {
    config$supplementary_tables <- list()
  }

  config$supplementary_tables[[table_name]] <- list(
    link_type = original_link_type, # Store what the user actually requested
    link_by = ifelse(is.null(link_by), "auto", link_by),
    abbreviation = abbreviation,
    date_added = Sys.Date(),
    n_rows = nrow(linked_data),
    n_cols = ncol(linked_data),
    columns = names(linked_data)
  )

  jsonlite::write_json(config, config_file, pretty = TRUE, auto_unbox = TRUE)
  message("Updated config file")

  invisible(TRUE)
}

#' Remove a Supplementary Table from locusPackRat Project
#'
#' Removes a supplementary table and its metadata from a project
#'
#' @param table_name Character: Name of the supplementary table to remove
#' @param project_dir Character: Path to project directory (default: ".")
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' # Remove a supplementary table
#' removeRatTable("old_data", project_dir = "my_analysis")
#' }
#'
#' @importFrom jsonlite read_json write_json
#'
#' @export
removeRatTable <- function(table_name, project_dir = ".") {
  # Validate project exists

  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) {
    stop("No .locusPackRat directory found. Run initPackRat() first.")
  }

  # Read config
  config_file <- file.path(packrat_dir, "config.json")
  if (!file.exists(config_file)) {
    stop("Config file not found. Project may be corrupted.")
  }
  config <- jsonlite::read_json(config_file)

  # Check table exists in config
  if (is.null(config$supplementary_tables[[table_name]])) {
    stop(sprintf("Table '%s' not found in project. Available tables: %s",
                 table_name,
                 paste(names(config$supplementary_tables), collapse = ", ")))
  }

  # Delete CSV file
  csv_path <- file.path(packrat_dir, "supplementary", paste0(table_name, ".csv"))
  if (file.exists(csv_path)) {
    file.remove(csv_path)
    message(sprintf("Removed file: %s", csv_path))
  } else {
    warning(sprintf("CSV file not found: %s", csv_path))
  }

  # Remove entry from config
  config$supplementary_tables[[table_name]] <- NULL
  jsonlite::write_json(config, config_file, pretty = TRUE, auto_unbox = TRUE)
  message(sprintf("Removed '%s' from config", table_name))

  invisible(TRUE)
}

#' Build Shareable Packet from locusPackRat Project
#'
#' Zips the contents of \code{.locusPackRat/output/} with an auto-generated
#' \code{README.md} summarizing the project, producing a single shareable archive.
#' Optionally includes raw supplementary CSV files from \code{.locusPackRat/supplementary/}.
#'
#' @param project_dir Character: Path to project directory containing .locusPackRat
#' @param output_path Character: Directory where the zip will be written (default: \code{project_dir})
#' @param filename Character: Name for the zip file (default: \code{"{basename}_packet.zip"})
#' @param include_readme Logical: Whether to include an auto-generated README.md (default: TRUE)
#' @param include_supplementary Logical or character vector: \code{FALSE} (default) includes
#'   only output files; \code{TRUE} includes all supplementary CSVs; a character vector
#'   (e.g., \code{c("ot_diseases", "mouse_phenotypes")}) includes only the named tables.
#'   Supplementary files are placed in a \code{supplementary/} subdirectory in the zip.
#' @param overwrite Logical: Whether to overwrite an existing zip file (default: FALSE)
#'
#' @return Invisible path to the created zip file
#'
#' @examples
#' \dontrun{
#' # Build a packet from the current project
#' buildPacket(project_dir = "my_analysis")
#'
#' # Custom output location
#' buildPacket(project_dir = "my_analysis", output_path = "~/shared", overwrite = TRUE)
#'
#' # Include all supplementary CSVs
#' buildPacket(project_dir = "my_analysis", include_supplementary = TRUE)
#'
#' # Include only specific supplementary tables
#' buildPacket(project_dir = "my_analysis",
#'             include_supplementary = c("ot_diseases", "mouse_phenotypes"))
#' }
#'
#' @importFrom jsonlite read_json
#' @importFrom utils zip
#'
#' @export
buildPacket <- function(project_dir = ".",
                        output_path = NULL,
                        filename = NULL,
                        include_readme = TRUE,
                        include_supplementary = FALSE,
                        overwrite = FALSE) {
  # Validate project
  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) {
    stop("No .locusPackRat directory found. Run initPackRat() first.")
  }

  config_file <- file.path(packrat_dir, "config.json")
  if (!file.exists(config_file)) {
    stop("Config file not found. Project may be corrupted.")
  }
  config <- jsonlite::read_json(config_file)

  # Check output directory has files
  output_dir <- file.path(packrat_dir, "output")
  if (!dir.exists(output_dir)) {
    stop("No output directory found. Run makeGeneSheet() first to generate output.")
  }
  output_files <- list.files(output_dir, full.names = TRUE)
  if (length(output_files) == 0) {
    stop("Output directory is empty. Run makeGeneSheet() first to generate output.")
  }

  # Resolve output path and filename
  if (is.null(output_path)) output_path <- project_dir
  if (!dir.exists(output_path)) {
    stop(sprintf("Output path does not exist: %s", output_path))
  }
  if (is.null(filename)) {
    filename <- paste0(basename(normalizePath(project_dir)), "_packet.zip")
  }
  zip_path <- file.path(normalizePath(output_path), filename)

  # Check overwrite

  if (file.exists(zip_path) && !overwrite) {
    stop(sprintf("File already exists: %s\nSet overwrite = TRUE to replace it.", zip_path))
  }

  # Create staging directory
  staging_dir <- file.path(tempdir(), paste0("packrat_packet_", format(Sys.time(), "%Y%m%d%H%M%S")))
  dir.create(staging_dir, recursive = TRUE)
  on.exit(unlink(staging_dir, recursive = TRUE), add = TRUE)

  # Copy output files to staging
  file.copy(output_files, staging_dir)

  # Handle supplementary files
  supp_copied <- NULL
  if (!identical(include_supplementary, FALSE)) {
    supp_dir <- file.path(packrat_dir, "supplementary")
    if (!dir.exists(supp_dir)) {
      warning("No supplementary directory found -- skipping supplementary files.")
    } else {
      all_csvs <- list.files(supp_dir, pattern = "\\.csv$", full.names = TRUE)
      if (is.character(include_supplementary)) {
        # Include only named tables
        requested <- paste0(include_supplementary, ".csv")
        matched <- all_csvs[basename(all_csvs) %in% requested]
        missing <- setdiff(requested, basename(matched))
        if (length(missing) > 0) {
          warning(sprintf("Supplementary table(s) not found: %s",
                          paste(sub("\\.csv$", "", missing), collapse = ", ")))
        }
        all_csvs <- matched
      }
      if (length(all_csvs) > 0) {
        staging_supp <- file.path(staging_dir, "supplementary")
        dir.create(staging_supp, recursive = TRUE)
        file.copy(all_csvs, staging_supp)
        supp_copied <- all_csvs
      }
    }
  }

  # Generate README
  if (include_readme) {
    readme_lines <- .generatePacketReadme(config, output_files, project_dir,
                                          supp_files = supp_copied)
    writeLines(readme_lines, file.path(staging_dir, "README.md"))
  }

  # Create zip from staging directory
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(staging_dir)

  files_to_zip <- list.files(".", full.names = FALSE, recursive = TRUE)

  # Create zip: prefer the 'zip' R package, fall back to system zip
  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zipr(zip_path, files = files_to_zip)
  } else {
    zip_cmd <- Sys.getenv("R_ZIPCMD")
    if (!nzchar(zip_cmd)) {
      zip_cmd <- Sys.which("zip")
      if (!nzchar(zip_cmd)) {
        stop(
          "No zip tool found. Either:\n",
          "  1. Install the 'zip' R package: install.packages('zip')\n",
          "  2. Install system zip (e.g., 'sudo apt install zip' on Linux)\n",
          "  3. On Windows, install Rtools: https://cran.r-project.org/bin/windows/Rtools/"
        )
      }
    }
    utils::zip(zip_path, files = files_to_zip, zip = zip_cmd)
  }

  zip_size <- .formatFileSize(file.info(zip_path)$size)
  message(sprintf("Packet created: %s (%s)", zip_path, zip_size))
  message(sprintf("Contains %d file(s)", length(files_to_zip)))

  invisible(zip_path)
}

#' List Supplementary Tables in locusPackRat Project
#'
#' Returns information about available supplementary tables in a project
#'
#' @param project_dir Character: Path to project directory containing .locusPackRat
#' @param full_info Boolean: If true, prints all column names for each supplementary table
#' @param criteria_info Boolean: If true, prints detailed column summaries useful for
#'   building upset plot criteria. Shows data types, ranges (numeric), and sample values (character).
#'
#' @return data.table with table information (name, link_type, link_by, n_rows, n_cols)
#'
#' @examples
#' \dontrun{
#' # List all supplementary tables
#' listPackRatTables(project_dir = "my_analysis")
#'
#' # List with full column information
#' listPackRatTables(project_dir = "my_analysis", full_info = TRUE)
#'
#' # List with criteria-building info (types, ranges, sample values)
#' listPackRatTables(project_dir = "my_analysis", criteria_info = TRUE)
#' }
#'
#' @importFrom data.table data.table fread
#' @importFrom dplyr bind_rows
#' @importFrom jsonlite read_json
#'
#' @export
listPackRatTables <- function(project_dir = ".", full_info = FALSE, criteria_info = FALSE) {
  # Check for existing project
  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) {
    stop("No .locusPackRat directory found. Run initPackRat() first.")
  }

  # Read config
  config_file <- file.path(packrat_dir, "config.json")
  if (!file.exists(config_file)) {
    warning("Config file not found. Project may be corrupted.")
    config <- list()
  } else {
    config <- jsonlite::read_json(config_file)
  }

  # Check supplementary directory
  supp_dir <- file.path(packrat_dir, "supplementary")
  if (!dir.exists(supp_dir)) {
    message("No supplementary tables found")
    return(data.table())
  }

  # Get list of tables
  supp_files <- list.files(supp_dir, pattern = "\\.csv$", full.names = FALSE)
  if (length(supp_files) == 0) {
    message("No supplementary tables found")
    return(data.table())
  }

  # Build table information
  table_info <- data.table()
  for (sf in supp_files) {
    table_name <- sub("\\.csv$", "", sf)

    # Get info from config if available
    if (!is.null(config$supplementary_tables) &&
        table_name %in% names(config$supplementary_tables)) {
      tbl_config <- config$supplementary_tables[[table_name]]
      info <- data.table(
        table_name = table_name,
        table_abbr = ifelse(is.null(tbl_config$abbreviation) | length(tbl_config$abbreviation)==0, 
          NA_character_, tbl_config$abbreviation),
        link_type = ifelse(is.null(tbl_config$link_type), NA_character_, tbl_config$link_type),
        link_by = ifelse(is.null(tbl_config$link_by), NA_character_, tbl_config$link_by),
        n_rows = ifelse(is.null(tbl_config$n_rows), NA_integer_, tbl_config$n_rows),
        n_cols = ifelse(is.null(tbl_config$n_cols), NA_integer_, tbl_config$n_cols),
        date_added = as.character(ifelse(is.null(tbl_config$date_added), NA_character_, tbl_config$date_added))
      )
    } else {
      # Try to get basic info from file
      file_path <- file.path(supp_dir, sf)
      tryCatch({
        temp_data <- fread(file_path, nrows = 1)
        info <- data.table(
          table_name = table_name,
          table_abbr = NA_character_,
          link_type = NA_character_,
          link_by = NA_character_,
          n_rows = as.integer(system(sprintf("wc -l < '%s'", file_path), intern = TRUE)) - 1L,
          n_cols = ncol(temp_data),
          date_added = NA_character_
        )
      }, error = function(e) {
        info <- data.table(
          table_name = table_name,
          table_abbr = NA_character_,
          link_type = NA_character_,
          link_by = NA_character_,
          n_rows = NA_integer_,
          n_cols = NA_integer_,
          date_added = NA_character_
        )
      })
    }
  #  message(sprintf("on %s table",sf))
    table_info <- bind_rows(table_info, info)
  }

  message(sprintf("Found %d supplementary table(s):", nrow(table_info)))
  for (i in seq_len(nrow(table_info))) {
    message(sprintf("  - %s: %d rows with %d cols, linked by '%s'",
                   table_info$table_name[i],
                   table_info$n_rows[i],
                   table_info$n_cols[i],
                   table_info$link_by[i]))
  }
  if (full_info) {
    message(sprintf("Printing full column names for %d supplementary table(s)",nrow(table_info)))
    for(i in 1:length(config$supplementary_tables)){
      name=names(config$supplementary_tables[i][1])
      columns=paste(unlist(config$supplementary_tables[i][[1]]$columns),collapse=" ; ")
      message(sprintf("Columns in %s:",name))
      message(sprintf(columns))
    }
    message(sprintf("Completed"))
  }

  # Criteria info mode: print detailed column summaries for building criteria
  if (criteria_info) {
    message("\n=== Column Details for Criteria Building ===\n")

    for (sf in supp_files) {
      table_name <- sub("\\.csv$", "", sf)
      file_path <- file.path(supp_dir, sf)

      # Get table metadata
      tbl_row <- table_info[table_info$table_name == table_name, ]
      link_info <- if (!is.na(tbl_row$link_by[1])) {
        sprintf(", linked by %s", tbl_row$link_by[1])
      } else {
        ""
      }

      message(sprintf("Table: %s (%d rows%s)",
                      table_name, tbl_row$n_rows[1], link_info))
      message("  Available columns for criteria:")

      # Load full table for column analysis
      tryCatch({
        table_data <- data.table::fread(file_path)

        for (col_name in names(table_data)) {
          col_data <- table_data[[col_name]]

          if (is.numeric(col_data)) {
            # Numeric: show min, max, median
            col_min <- signif(min(col_data, na.rm = TRUE), 4)
            col_max <- signif(max(col_data, na.rm = TRUE), 4)
            col_med <- signif(stats::median(col_data, na.rm = TRUE), 4)
            n_na <- sum(is.na(col_data))
            na_info <- if (n_na > 0) sprintf(", %d NA", n_na) else ""
            message(sprintf("    - %s (numeric): min=%s, max=%s, median=%s%s",
                            col_name, col_min, col_max, col_med, na_info))

          } else if (is.logical(col_data)) {
            # Logical: show TRUE/FALSE counts
            n_true <- sum(col_data, na.rm = TRUE)
            n_false <- sum(!col_data, na.rm = TRUE)
            n_na <- sum(is.na(col_data))
            na_info <- if (n_na > 0) sprintf(", %d NA", n_na) else ""
            message(sprintf("    - %s (logical): %d TRUE, %d FALSE%s",
                            col_name, n_true, n_false, na_info))

          } else {
            # Character: show unique count and sample values
            unique_vals <- unique(col_data[!is.na(col_data) & col_data != ""])
            n_unique <- length(unique_vals)
            n_na <- sum(is.na(col_data) | col_data == "")

            if (n_unique <= 6) {
              # Few unique values - show all
              vals_str <- paste(sprintf("\"%s\"", unique_vals), collapse = ", ")
              message(sprintf("    - %s (character): %d unique values: %s",
                              col_name, n_unique, vals_str))
            } else {
              # Many unique values - show first 3 and last 3
              sorted_vals <- sort(unique_vals)
              first3 <- paste(sprintf("\"%s\"", head(sorted_vals, 3)), collapse = ", ")
              last3 <- paste(sprintf("\"%s\"", tail(sorted_vals, 3)), collapse = ", ")
              message(sprintf("    - %s (character): %d unique values",
                              col_name, n_unique))
              message(sprintf("        First 3: %s", first3))
              message(sprintf("        Last 3: %s", last3))
            }
          }
        }
        message("")  # Blank line between tables

      }, error = function(e) {
        message(sprintf("    Error reading table: %s", e$message))
      })
    }
  }

  return(table_info)
}

#' Collapse Column Values for Merging
#' @noRd
.collapse_column <- function(x) {
  # Drop NAs and empty strings (if character)
  if (is.character(x)) {
    x <- x[x != ""]
  }
  x <- unique(x[!is.na(x)])

  # If no values left, return NA of the *same type* as the column
  if (length(x) == 0L) {
    return(x[NA_integer_])   # NA with same typeof(x)
  }

  # If only one value, just return it (keeps original type)
  if (length(x) == 1L) {
    return(x)
  }

  # If character and multiple values, paste them
  if (is.character(x)) {
    return(paste(x, collapse = "; "))
  }

  # For non-character with multiple distinct values, just keep the first
  # (keeps type numeric/logical/etc, avoids mixing with character)
  x[1L]
}

#' Generate Gene Sheet from locusPackRat Project
#'
#' Creates formatted gene sheets (CSV or Excel) from project data with optional filtering and highlighting
#'
#' @param filter_expr Character: Expression to filter genes (uses data.table syntax, e.g., "chr == 1 & start > 1000000")
#' @param highlight_genes Character vector: Gene symbols to highlight
#' @param highlight_expr Character: Expression to determine which genes to highlight
#' @param highlight_color Character: Background color for highlighted genes (hex color)
#' @param summary_level Character: "gene" for one row per gene, "locus" for one row per locus/region
#' @param output_file Character: Output file path (defaults to .locusPackRat/output/gene_sheet.csv or .xlsx)
#' @param format Character: "csv" for plain CSV, "excel" for formatted Excel workbook
#' @param split_by Character: For Excel, how to split into tabs - "none", "chromosome", "criteria"
#' @param split_criteria Named list: For split_by="criteria", list of tab_name = filter_expression
#' @param include_supplementary Logical or character vector: TRUE to include all supplementary tables,
#'                               FALSE for none, or character vector of specific table names to include
#' @param collapse_multiples Logical: when merging supplementary tables,
#'   collapse multiple matches per key into a single row by aggregating
#'   values (e.g. concatenating text). If FALSE, keep all matches as
#'   separate rows.
#' @param prefix_mode Character: How to handle column name prefixing from supplementary tables.
#'                    "collision" (default) only prefixes columns that would overwrite existing ones.
#'                    "abbreviated" prefixes all columns from tables that have an abbreviation set.
#'                    "always" prefixes all columns, using table_name as fallback if no abbreviation.
#' @param project_dir Character: Path to project directory containing .locusPackRat
#' @param exclude_tables Character vector: Names of supplementary tables to exclude from output
#' @param format_numbers Logical: Auto-detect and format numeric columns (default TRUE)
#' @param format_pvalues Logical: Apply scientific notation to p-value columns (default TRUE)
#' @param format_coordinates Logical: Apply thousands separator to coordinate columns (default TRUE)
#' @param include_summary Logical: Add a summary sheet with workbook metadata (default TRUE)
#' @param header_color Character: Hex color for header background (default "#4472C4" - Excel blue)
#' @param min_col_width Numeric: Minimum column width in characters (default 10)
#' @param merge_repeated Logical: Merge consecutive cells with same gene_symbol in Phenotypes_Detail sheet (default FALSE)
#'
#' @return data.table of the gene sheet (invisible)
#'
#' @examples
#' \dontrun{
#' # Basic CSV export
#' makeGeneSheet(format = "csv", output_file = "results.csv",
#'               project_dir = "my_analysis")
#'
#' # Filtered Excel export with multiple sheets
#' makeGeneSheet(
#'   format = "excel",
#'   output_file = "analysis.xlsx",
#'   split_by = "criteria",
#'   split_criteria = list(
#'     "All_Genes" = "TRUE",
#'     "Significant" = "padj < 0.05",
#'     "Upregulated" = "log2FC > 1"
#'   ),
#'   highlight_genes = c("Myc", "Tp53"),
#'   project_dir = "my_analysis"
#' )
#' }
#'
#' @importFrom data.table fread fwrite setDT merge.data.table := .SD
#' @importFrom jsonlite read_json
#' @importFrom openxlsx createWorkbook addWorksheet writeDataTable createStyle addStyle setColWidths freezePane saveWorkbook
#'
#' @export
makeGeneSheet <- function(filter_expr = NULL,
                          highlight_genes = NULL,
                          highlight_expr = NULL,
                          highlight_color = "#FFFF99",
                          summary_level = c("gene", "locus"),
                          output_file = NULL,
                          format = c("excel", "csv"),
                          split_by = c("none", "chromosome", "criteria"),
                          split_criteria = NULL,
                          include_supplementary = TRUE,
                          exclude_tables = NULL,
                          collapse_multiples = TRUE,
                          prefix_mode = c("collision", "abbreviated", "always"),
                          project_dir = ".",
                          format_numbers = TRUE,
                          format_pvalues = TRUE,
                          format_coordinates = TRUE,
                          include_summary = TRUE,
                          header_color = "#4472C4",
                          min_col_width = 10,
                          merge_repeated = FALSE) {

  summary_level <- match.arg(summary_level)
  format <- match.arg(format)
  split_by <- match.arg(split_by)
  prefix_mode <- match.arg(prefix_mode)

  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) stop(sprintf("No .locusPackRat directory found at %s", packrat_dir))

  config <- jsonlite::read_json(file.path(packrat_dir, "config.json"))
  message(sprintf("Generating %s sheet from %s...", config$mode, packrat_dir))

  # 1. LOAD BASE DATA
  if (config$mode == "gene") {
    base_data <- fread(file.path(packrat_dir, "input", "genes.csv"))
    if(!"pk_id" %in% names(base_data)) base_data[, pk_id := .I]
  } else {
    region_data <- fread(file.path(packrat_dir, "input", "regions.csv"))
    base_data <- .expandRegionsToGenes(region_data, config)
    if(!"pk_id" %in% names(base_data)) base_data[, pk_id := .I]
  }

  # 2. LOAD ORTHOLOGY
  ortho_file <- file.path(packrat_dir, "input", "orthology.csv")
  if (file.exists(ortho_file)) {
    ortho_data <- fread(ortho_file)
    if ("ortho_species" %in% names(ortho_data)) ortho_data[, ortho_species := NULL]

    if (config$species == "mouse") {
      mx <- "gene_symbol"; my <- "mouse_gene_symbol"; keep <- c("human_gene_symbol", "human_ensembl_id")
    } else {
      mx <- "gene_symbol"; my <- "human_gene_symbol"; keep <- c("mouse_gene_symbol", "mouse_ensembl_id")
    }

    if (mx %in% names(base_data) && my %in% names(ortho_data)) {
      base_data <- merge(base_data, ortho_data[, c(my, keep), with=FALSE], by.x = mx, by.y = my, all.x = TRUE)
    }
  }

  # 3. LINK SUPPLEMENTARY TABLES
  if (!isFALSE(include_supplementary)) {
    supp_dir <- file.path(packrat_dir, "supplementary")
    
    if (dir.exists(supp_dir)) {
      all_supp <- list.files(supp_dir, pattern = "\\.csv$", full.names = FALSE)
      
      # Debug: Check what files are found
      if(length(all_supp) == 0) {
        warning(sprintf("No CSV files found in supplementary directory: %s", supp_dir))
      } else {
        message(sprintf("Found %d supplementary files: %s", length(all_supp), paste(all_supp, collapse=", ")))
      }

      if (is.character(include_supplementary)) {
        tables_to_include <- intersect(paste0(include_supplementary, ".csv"), all_supp)
      } else {
        tables_to_include <- all_supp
      }

      # Filter out excluded tables
      if (!is.null(exclude_tables)) {
        tables_to_include <- setdiff(tables_to_include, paste0(exclude_tables, ".csv"))
        if (length(exclude_tables) > 0) {
          message(sprintf("Excluding tables: %s", paste(exclude_tables, collapse = ", ")))
        }
      }

      for (sf in tables_to_include) {
        table_name <- sub("\\.csv$", "", sf)
        supp_data <- fread(file.path(supp_dir, sf))
        
        # --- A. DETECT KEYS ---
        merge_cols <- NULL
        # Check Config
        if (!is.null(config$supplementary_tables[[table_name]]$link_by) && config$supplementary_tables[[table_name]]$link_by != "auto") {
          merge_cols <- trimws(strsplit(config$supplementary_tables[[table_name]]$link_by, ",")[[1]])
        }
        # Fallback Auto-detect
        if (is.null(merge_cols)) {
          common <- intersect(names(base_data), names(supp_data))
          if ("region_id" %in% common) merge_cols <- "region_id"
          else if (all(c("chr", "start", "end") %in% common)) merge_cols <- c("chr", "start", "end")
          else if ("gene_symbol" %in% common) merge_cols <- "gene_symbol"
        }

        # --- COLUMN PREFIXING ---
        # Get abbreviation from config (or use table_name as fallback for "always" mode)
        abbrev <- config$supplementary_tables[[table_name]]$abbreviation
        val_cols <- setdiff(names(supp_data), merge_cols)

        # Determine which columns to prefix based on mode
        cols_to_prefix <- character(0)
        prefix <- NULL

        if (prefix_mode == "always") {
          # Use abbreviation if set, otherwise table_name
          prefix <- if (!is.null(abbrev)) abbrev else table_name
          cols_to_prefix <- val_cols
        } else if (prefix_mode == "abbreviated" && !is.null(abbrev)) {
          # Prefix all columns from tables that have abbreviation
          prefix <- abbrev
          cols_to_prefix <- val_cols
        } else if (prefix_mode == "collision" && !is.null(abbrev)) {
          # Only prefix columns that would collide
          prefix <- abbrev
          cols_to_prefix <- intersect(names(base_data), val_cols)
        }

        if (length(cols_to_prefix) > 0 && !is.null(prefix)) {
          new_names <- paste0(prefix, "_", cols_to_prefix)
          setnames(supp_data, cols_to_prefix, new_names)
          message(sprintf("    Prefixed %d columns with '%s_'", length(cols_to_prefix), prefix))
        }

        # --- B. EXECUTE MERGE LOGIC ---
        if (!is.null(merge_cols)) {
          # Capture columns we EXPECT to add (for safety filling later)
          expected_cols <- setdiff(names(supp_data), merge_cols)

          # CASE 1: COORDINATE OVERLAP (Interval Join)
          if (all(c("chr", "start", "end") %in% merge_cols)) {
             message(sprintf("  Linking %s using Interval Overlap...", table_name))
             
             # Normalize Coordinates
             supp_data[, chr := gsub("^chr", "", as.character(chr), ignore.case=TRUE)]
             supp_data[, start := as.integer(start)][, end := as.integer(end)]
             base_data[, chr := gsub("^chr", "", as.character(chr), ignore.case=TRUE)]
             base_data[, start := as.integer(start)][, end := as.integer(end)]

             # Non-Equi Join
             # Note: This join might produce 0 rows if no overlaps
             overlap_join <- supp_data[base_data, 
                                       on = .(chr, start <= end, end >= start), 
                                       nomatch = NULL] 
             
             val_cols <- setdiff(names(supp_data), c("chr", "start", "end"))
             
             if (length(val_cols) > 0) {
               # Aggregation
              collapsed_vals <- overlap_join[
                ,
                lapply(.SD, .collapse_column),
                by = pk_id,
                .SDcols = val_cols
              ]
               
               # Merge
               base_data <- merge(base_data, collapsed_vals, by = "pk_id", all.x = TRUE)
             }

          } else {
            # CASE 2: EXACT MERGE
            message(sprintf("  Linking %s (Exact Match: %s)", table_name, paste(merge_cols, collapse=",")))
            
            if (collapse_multiples) {
              val_cols <- setdiff(names(supp_data), merge_cols)
              if (length(val_cols) > 0) {
                supp_data <- supp_data[
                  ,
                  lapply(.SD, .collapse_column),
                  by = merge_cols,
                  .SDcols = val_cols
                ]
              } else {
                supp_data <- unique(supp_data)
              }
            }
            base_data <- merge(base_data, supp_data, by = merge_cols, all.x = TRUE, suffixes = c("", paste0(".", table_name)))
          }
          
          # SAFETY CHECK: Ensure expected columns exist
          # If merge failed or yielded 0 matches, columns might be missing.
          # We force add them as NA to prevent 'Object not found' errors in filtering.
          missing_cols <- setdiff(expected_cols, names(base_data))
          if(length(missing_cols) > 0) {
             message(sprintf("    Warning: %s matched 0 rows. Adding empty columns: %s", table_name, paste(missing_cols, collapse=", ")))
             base_data[, (missing_cols) := NA]
          }

          # Clean suffixes
          bad_cols <- grep(paste0("\\.", table_name, "$"), names(base_data), value = TRUE)
          if(length(bad_cols) > 0) base_data[, (bad_cols) := NULL]

        } else {
          warning(sprintf("  Skipping %s: No link columns found.", table_name))
        }
      } 
    }
  }

  if("pk_id" %in% names(base_data)) base_data[, pk_id := NULL]

  # 4. FILTERING
  if (!is.null(filter_expr)) {
    tryCatch({
      base_data <- base_data[eval(parse(text = filter_expr))]
    }, error = function(e) stop("Filter Expression Error: ", e$message))
  }

  # 5. SUMMARIZE & SAVE
  output_data <- if (summary_level == "locus") .summarizeAtLocusLevel(base_data, config) else base_data

  output_dir <- file.path(packrat_dir, "output")
  if (is.null(output_file)) output_file <- paste0("gene_sheet", ifelse(format == "excel", ".xlsx", ".csv"))
  final_path <- file.path(output_dir, basename(output_file))

  if (format == "csv") {
    fwrite(output_data, final_path)
    message(sprintf("Saved CSV to %s", final_path))
  } else {
    wb <- .createFormattedExcel(
      data = output_data,
      highlight_genes = highlight_genes,
      highlight_expr = highlight_expr,
      highlight_color = highlight_color,
      split_by = split_by,
      split_criteria = split_criteria,
      format_numbers = format_numbers,
      format_pvalues = format_pvalues,
      format_coordinates = format_coordinates,
      include_summary = include_summary,
      header_color = header_color,
      min_col_width = min_col_width,
      filter_expr = filter_expr,
      merge_repeated = merge_repeated
    )
    openxlsx::saveWorkbook(wb, final_path, overwrite = TRUE)
    message(sprintf("Saved Excel to %s", final_path))
  }

  invisible(output_data)
}

# ========== Helper Functions ==========

#' Generate README.md Content for a Packet
#' @param config List: parsed config.json
#' @param output_files Character: full paths to output files
#' @param project_dir Character: path to project directory
#' @param supp_files Character or NULL: full paths to included supplementary CSVs
#' @return Character vector of markdown lines
#' @noRd
.generatePacketReadme <- function(config, output_files, project_dir, supp_files = NULL) {
  safe <- function(x) if (is.null(x)) "unknown" else as.character(x)

  lines <- c(
    sprintf("# %s -- locusPackRat Packet", basename(normalizePath(project_dir))),
    "",
    "## Project Info",
    "",
    "| Field | Value |",
    "| ----- | ----- |",
    sprintf("| Mode | %s |", safe(config$mode)),
    sprintf("| Species | %s |", safe(config$species)),
    sprintf("| Genome | %s |", safe(config$genome)),
    sprintf("| Entries | %s |", safe(config$n_entries)),
    sprintf("| Initialized | %s |", safe(config$initialization_date)),
    sprintf("| Package version | %s |", safe(config$package_version)),
    ""
  )

  # Supplementary data sources
  tbl_info <- tryCatch(
    suppressMessages(listPackRatTables(project_dir)),
    error = function(e) NULL
  )
  if (!is.null(tbl_info) && nrow(tbl_info) > 0) {
    lines <- c(lines,
      "## Supplementary Data Sources",
      "",
      "| Table | Abbreviation | Link Type | Rows | Cols | Date Added |",
      "| ----- | ------------ | --------- | ---- | ---- | ---------- |"
    )
    for (i in seq_len(nrow(tbl_info))) {
      r <- tbl_info[i, ]
      lines <- c(lines, sprintf("| %s | %s | %s | %s | %s | %s |",
        r$table_name,
        ifelse(is.na(r$table_abbr), "-", r$table_abbr),
        ifelse(is.na(r$link_type), "-", r$link_type),
        ifelse(is.na(r$n_rows), "-", as.character(r$n_rows)),
        ifelse(is.na(r$n_cols), "-", as.character(r$n_cols)),
        ifelse(is.na(r$date_added), "-", r$date_added)
      ))
    }
    lines <- c(lines, "")
  }

  # Included supplementary files (bundled in the zip)
  if (!is.null(supp_files) && length(supp_files) > 0) {
    sinfo <- file.info(supp_files)
    lines <- c(lines,
      "## Included Supplementary Files",
      "",
      "| File | Size |",
      "| ---- | ---- |"
    )
    for (i in seq_along(supp_files)) {
      lines <- c(lines, sprintf("| supplementary/%s | %s |",
        basename(supp_files[i]),
        .formatFileSize(sinfo$size[i])
      ))
    }
    lines <- c(lines, "")
  }

  # Included files
  finfo <- file.info(output_files)
  lines <- c(lines,
    "## Included Files",
    "",
    "| File | Size |",
    "| ---- | ---- |"
  )
  for (i in seq_along(output_files)) {
    lines <- c(lines, sprintf("| %s | %s |",
      basename(output_files[i]),
      .formatFileSize(finfo$size[i])
    ))
  }

  # Footer
  pkg_version <- tryCatch(
    as.character(packageVersion("locusPackRat")),
    error = function(e) "dev"
  )
  lines <- c(lines, "",
    "---",
    sprintf("*Generated by locusPackRat %s on %s*", pkg_version, Sys.time())
  )

  lines
}

#' Format Bytes to Human-Readable Size
#' @param bytes Numeric: file size in bytes
#' @return Character string (e.g., "1.2 MB")
#' @noRd
.formatFileSize <- function(bytes) {
  if (is.na(bytes) || bytes == 0) return("0 B")
  units <- c("B", "KB", "MB", "GB")
  exp <- min(floor(log(bytes, 1024)), length(units) - 1)
  val <- bytes / (1024^exp)
  if (exp == 0) return(sprintf("%d B", as.integer(val)))
  sprintf("%.1f %s", val, units[exp + 1])
}

#' Process Gene Input Data
#' @noRd
.processGeneInput <- function(data, species, genome) {
  # Make a defensive copy to avoid modifying the caller's data
  data <- copy(data)

  # Detect ID type based on column names
  gene_cols <- c("gene_symbol", "gene_name", "symbol", "ensembl_id", "gene_id", "id")
  id_col <- intersect(names(data), gene_cols)[1]

  if (is.na(id_col)) {
    stop("Could not find gene identifier column. Expected one of: ",
         paste(gene_cols, collapse = ", "))
  }

  # Standardize column names
  setnames(data, id_col, "input_id", skip_absent = TRUE)

  # Load genome coordinate data from cached files
  coords_file <- NULL
  if (species == "mouse") {
    if (genome == "mm39") {
      coords_file <- system.file("extdata", "mouse_coords_mm39.csv", package = "locusPackRat")
      if (coords_file == "") {
        # Try test data location if package not installed
        coords_file <- "tests/test_data/mouse_coords_mm39.csv"
      }
    }
    if (genome == "mm10") {
      coords_file <- system.file("extdata", "mouse_coords_mm10.csv", package = "locusPackRat")
      if (coords_file == "") {
        # Try test data location if package not installed
        coords_file <- "tests/test_data/mouse_coords_mm10.csv"
      }
    }
  } else if (species == "human") {
    if (genome == "hg38") {
      coords_file <- system.file("extdata", "human_coords_hg38.csv", package = "locusPackRat")
      if (coords_file == "") {
        # Try test data location if package not installed
        coords_file <- "tests/test_data/human_coords_hg38.csv"
      }
    }
    if (genome == "hg19") {
      coords_file <- system.file("extdata", "human_coords_hg19.csv", package = "locusPackRat")
      if (coords_file == "") {
        # Try test data location if package not installed
        coords_file <- "tests/test_data/human_coords_hg19.csv"
      }
    }
  }

  # Initialize result with input gene symbols
  result <- data.table(
    gene_symbol = data$input_id,
    ensembl_id = NA_character_,
    chr = NA_character_,
    start = NA_integer_,
    end = NA_integer_,
    strand = NA_character_
  )

  # Load coordinates if file exists
  if (!is.null(coords_file) && file.exists(coords_file)) {
    coords_data <- fread(coords_file)

    # Match on gene_symbol or ensembl_id
    id_type <- ifelse(grepl("^ENS", data$input_id[1]), "ensembl_id", "gene_symbol")

    if (id_type == "gene_symbol") {
      # Merge by gene symbol
      result <- merge(result[, .(gene_symbol)], coords_data,
                     by = "gene_symbol", all.x = TRUE)
    } else {
      # Merge by ensembl ID
      result[, ensembl_id := data$input_id]
      result <- merge(result[, .(ensembl_id)], coords_data,
                     by = "ensembl_id", all.x = TRUE)
    }

    # Ensure numeric coordinates
    if ("start" %in% names(result)) result[, start := as.integer(start)]
    if ("end" %in% names(result)) result[, end := as.integer(end)]

    # Report match statistics
    n_total <- nrow(result)
    n_matched <- sum(!is.na(result$chr))
    n_unmatched <- n_total - n_matched
    if (n_unmatched > 0) {
      unmatched_genes <- result[is.na(chr), gene_symbol]
      if (length(unmatched_genes) > 5) {
        unmatched_preview <- paste(c(head(unmatched_genes, 5), "..."), collapse = ", ")
      } else {
        unmatched_preview <- paste(unmatched_genes, collapse = ", ")
      }
      message(sprintf("  Matched %d/%d genes to coordinates (%d not found: %s)",
                      n_matched, n_total, n_unmatched, unmatched_preview))
    } else {
      message(sprintf("  Matched all %d genes to coordinates", n_total))
    }
  } else {
    message("  Note: Coordinate data file not found, using NA placeholders")
  }

  # Preserve any additional columns from input data (but not the original id column)
  extra_cols <- setdiff(names(data), c(id_col, "input_id"))
  if (length(extra_cols) > 0) {
    # Merge back the extra columns
    data_with_extras <- copy(data)
    data_with_extras[, gene_symbol := input_id]
    result <- merge(result, data_with_extras[, c("gene_symbol", extra_cols), with = FALSE],
                   by = "gene_symbol", all.x = TRUE)
  }

  return(result)
}

#' Process Region Input Data
#' @noRd
.processRegionInput <- function(data, species, genome, keep_pseudo, overlap_mode = "any") {
  # Make a defensive copy to avoid modifying the caller's data
  data <- copy(data)

  # Check for required columns
  required <- c("chr", "start", "end")
  missing <- setdiff(required, names(data))

  if (length(missing) > 0) {
    # Try alternative names
    alt_names <- list(
      chr = c("chrom", "chromosome", "seqname"),
      start = c("start_pos", "chromStart", "begin"),
      end = c("end_pos", "chromEnd", "stop")
    )

    for (req in missing) {
      alts <- alt_names[[req]]
      found <- intersect(names(data), alts)
      if (length(found) > 0) {
        setnames(data, found[1], req)
      }
    }
  }
  # Final check
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns for regions: ", paste(missing, collapse = ", "))
  }

  # Ensure numeric coordinates and character chromosome
  data[, start := as.numeric(start)]
  data[, end := as.numeric(end)]
  data[, chr := as.character(chr)]

  # Add region ID if not present
  if (!"region_id" %in% names(data)) {
    data[, region_id := paste0("region_", .I)]
  }
  # Get genes in regions
  coord_type <- paste0(species, "_coords_", genome, ".csv")
  coords_file <- system.file("extdata", coord_type, package = "locusPackRat")
  coords <- fread(coords_file)
  if (!keep_pseudo) {
    # Filter out Gm* genes and genes ending in Rik (both are pseudogenes)
    pseudo_pattern <- "^Gm[0-9]|Rik$"
    pseudo_idx <- grep(pseudo_pattern, coords$gene_symbol)
    if (length(pseudo_idx) > 0) {
      coords <- coords[-pseudo_idx]
    }
  }

  # Ensure numeric coordinates
  coords[, start := as.numeric(start)]
  coords[, end := as.numeric(end)]
  coords[, chr := as.character(chr)]

  # Look for overlaps - capture gene coordinates to determine overlap type
  overlaps <- coords[data,
                      .(region_id = i.region_id,
                        region_start = i.start,
                        region_end = i.end,
                        gene_symbol = x.gene_symbol,
                        gene_start = x.start,
                        gene_end = x.end),
                      on = .(chr, start <= end, end >= start),
                      nomatch = NULL]

  if (nrow(overlaps) > 0) {
      # Calculate overlap type: complete if gene is entirely within region
      overlaps[, overlap_type := ifelse(
        gene_start >= region_start & gene_end <= region_end,
        "complete",
        "partial"
      )]

      # Filter if overlap_mode == "complete"
      if (overlap_mode == "complete") {
        overlaps <- overlaps[overlap_type == "complete"]
      }

      if (nrow(overlaps) > 0) {
        # Collapse genes and overlap_types into comma-separated strings (in same order)
        genes_aggregated <- overlaps[, .(
          genes = paste(unique(gene_symbol), collapse = ", "),
          overlap_types = paste(overlap_type[!duplicated(gene_symbol)], collapse = ", ")
        ), by = region_id]

        # Merge back to the original data
        data <- merge(data, genes_aggregated, by = "region_id", all.x = TRUE)
      } else {
        data[, genes := NA_character_]
        data[, overlap_types := NA_character_]
      }

  } else {
      # If no overlaps found, create empty columns
      data[, genes := NA_character_]
      data[, overlap_types := NA_character_]
  }
  return(data)
}

#' Generate Orthology Table
#' @noRd
.generateOrthologyTable <- function(data, species, genome) {
  # Load orthology data from cached file
  ortho_file <- system.file("extdata", "orthology.csv", package = "locusPackRat")
  if (ortho_file == "") {
    # Try test data location if package not installed
    ortho_file <- "tests/test_data/orthology_mapping.csv"
  }

  if (!file.exists(ortho_file)) {
    message("  Note: Orthology data file not found, using NA placeholders")

    # Fallback to placeholder data
    if ("gene_symbol" %in% names(data)) {
      genes <- unique(data$gene_symbol)
      genes <- genes[!is.na(genes)]

      if (length(genes) > 0) {
        ortho <- data.table(
          source_gene_symbol = genes,
          source_ensembl_id = NA_character_,
          ortho_gene_symbol = NA_character_,
          ortho_ensembl_id = NA_character_,
          ortho_species = ifelse(species == "human", "mouse", "human")
        )

        # Add species-specific column names
        if (species == "human") {
          setnames(ortho,
                  c("source_gene_symbol", "source_ensembl_id", "ortho_gene_symbol", "ortho_ensembl_id"),
                  c("human_gene_symbol", "human_ensembl_id", "mouse_gene_symbol", "mouse_ensembl_id"))
        } else {
          setnames(ortho,
                  c("source_gene_symbol", "source_ensembl_id", "ortho_gene_symbol", "ortho_ensembl_id"),
                  c("mouse_gene_symbol", "mouse_ensembl_id", "human_gene_symbol", "human_ensembl_id"))
        }

        return(ortho)
      }
    }
    return(NULL)
  }

  # Load orthology mapping data
  ortho_mapping <- fread(ortho_file)

  # Remove duplicates if any
  ortho_mapping <- unique(ortho_mapping)

  # Get unique genes from input data
  if ("gene_symbol" %in% names(data)) {
    genes <- unique(data$gene_symbol)
    genes <- genes[!is.na(genes)]
  } else if ("genes" %in% names(data)) {
    # Region mode: extract genes from comma-separated column
    genes_raw <- data$genes[!is.na(data$genes) & data$genes != ""]
    genes <- unique(unlist(strsplit(genes_raw, ", ")))
    genes <- trimws(genes)
    genes <- genes[genes != ""]
  } else if ("ensembl_id" %in% names(data)) {
    ensembl_ids <- unique(data$ensembl_id)
    ensembl_ids <- ensembl_ids[!is.na(ensembl_ids)]
    # We'll match on Ensembl IDs instead
    genes <- NULL
  } else {
    return(NULL)
  }

  # Filter orthology data for our genes
  if (species == "mouse") {
    # Filter for mouse genes
    if (!is.null(genes)) {
      result <- ortho_mapping[mouse_gene_symbol %in% genes]
    } else {
      result <- ortho_mapping[mouse_ensembl_id %in% ensembl_ids]
    }

    # Ensure we have all input genes in result (even if no ortholog)
    if (!is.null(genes)) {
      missing_genes <- setdiff(genes, result$mouse_gene_symbol)
      if (length(missing_genes) > 0) {
        missing_rows <- data.table(
          mouse_gene_symbol = missing_genes,
          mouse_ensembl_id = NA_character_,
          human_gene_symbol = NA_character_,
          human_ensembl_id = NA_character_
        )
        result <- bind_rows(result, missing_rows)
      }
    }
  } else {
    # Filter for human genes
    if (!is.null(genes)) {
      result <- ortho_mapping[human_gene_symbol %in% genes]
    } else {
      result <- ortho_mapping[human_ensembl_id %in% ensembl_ids]
    }

    # Ensure we have all input genes in result (even if no ortholog)
    if (!is.null(genes)) {
      missing_genes <- setdiff(genes, result$human_gene_symbol)
      if (length(missing_genes) > 0) {
        missing_rows <- data.table(
          human_gene_symbol = missing_genes,
          human_ensembl_id = NA_character_,
          mouse_gene_symbol = NA_character_,
          mouse_ensembl_id = NA_character_
        )
        result <- bind_rows(result, missing_rows)
      }
    }
  }

  # Remove the ortho_species column if it exists (it's redundant)
  if ("ortho_species" %in% names(result)) {
    result[, ortho_species := NULL]
  }

  return(result)
}

#' Link Gene Data
#' @noRd
.linkGeneData <- function(new_data, existing_genes, link_by = NULL) {
  if (is.null(link_by)) {
    # Auto-detect linking column
    possible_links <- c("gene_symbol", "ensembl_id", "gene_name", "gene_id")
    link_cols <- intersect(names(new_data), possible_links)

    if (length(link_cols) == 0) {
      stop("Could not auto-detect linking column. Please specify link_by parameter.")
    }
    link_by <- link_cols[1]
  }

  # Check if link column exists in both datasets
  if (!link_by %in% names(new_data)) {
    stop(sprintf("Link column '%s' not found in input data", link_by))
  }

  # Find matching column in existing data
  existing_link <- intersect(names(existing_genes), c(link_by, "gene_symbol", "ensembl_id"))
  if (length(existing_link) == 0) {
    stop("Could not find matching column in existing gene data")
  }

  # Perform merge to find matching rows
  # But only keep the linking column and new data columns
  matched <- merge(existing_genes[, ..existing_link], new_data,
                 by.x = existing_link[1], by.y = link_by,
                 all = FALSE)

  # Rename the linking column back to its original name if needed
  if (existing_link[1] != link_by) {
    setnames(matched, existing_link[1], link_by)
  }

  return(matched)
}

#' Link Region Data
#' @importFrom data.table foverlaps setkey
#' @noRd
.linkRegionData <- function(new_data, existing_regions, genome) {
  # Check for coordinate columns in new data
  coord_cols <- c("chr", "start", "end")
  if (!all(coord_cols %in% names(new_data))) {
    stop("Region linking requires chr, start, end columns in input data")
  }

  # Convert to numeric
  new_data[, start := as.numeric(start)]
  new_data[, end := as.numeric(end)]
  existing_regions[, start := as.numeric(start)]
  existing_regions[, end := as.numeric(end)]

  # Force chromosome to character to avoid integer/character join errors
  # (fread often guesses integer for simple chromosomes like "1", "2")
  new_data[, chr := as.character(chr)]
  existing_regions[, chr := as.character(chr)]

  # Find overlaps
  setkey(new_data, chr, start, end)
  setkey(existing_regions, chr, start, end)

  # Use data.table's foverlaps for efficient overlap detection
  overlaps <- data.table::foverlaps(new_data, existing_regions, nomatch = NULL)

  return(overlaps)
}

#' Expand Regions to Genes
#' @noRd
.expandRegionsToGenes <- function(region_data, config) {
  if (!"genes" %in% names(region_data) || all(is.na(region_data$genes)) || all(region_data$genes == "")) {
    warning("`genes` column is missing or empty, returning empty gene table")
    return(data.table())
  }

  region_data[, genes := as.character(genes)]
  # Split comma-separated genes
  genes_list <- strsplit(region_data$genes, ", ")

  # Split overlap_types if present (same order as genes)
  has_overlap_types <- "overlap_types" %in% names(region_data)
  if (has_overlap_types) {
    region_data[, overlap_types := as.character(overlap_types)]
    overlap_types_list <- strsplit(region_data$overlap_types, ", ")
  }

  # Create expanded table (keep regions even if they have no genes)
  n_genes <- lengths(genes_list)
  n_genes[n_genes == 0] <- 1  # Ensure at least 1 row per region
  result <- region_data[rep(seq_len(nrow(region_data)), n_genes)]

  # Assign gene symbols, using NA for regions with no genes
  gene_symbols <- lapply(genes_list, function(x) if(length(x) == 0) NA_character_ else x)
  result[, gene_symbol := unlist(gene_symbols)]
  result[, gene_symbol := trimws(gene_symbol)]

  # Add overlap_type for each gene (NA for regions with no genes)
  if (has_overlap_types) {
    overlap_types_expanded <- lapply(overlap_types_list, function(x) if(length(x) == 0) NA_character_ else x)
    result[, overlap_type := unlist(overlap_types_expanded)]
    result[, overlap_type := trimws(overlap_type)]
    # Remove the aggregated column
    result[, overlap_types := NULL]
  }

  # Add region info and rename region coordinates
  result[, from_region := paste0(chr, ":", start, "-", end)]
  setnames(result, c("chr", "start", "end"), c("region_chr", "region_start", "region_end"))

  # Load gene-specific coordinates from genome annotation
  species <- config$species
  genome <- config$genome
  coords_file <- NULL

  if (species == "mouse") {
    if (genome == "mm39") {
      coords_file <- system.file("extdata", "mouse_coords_mm39.csv", package = "locusPackRat")
      if (coords_file == "") coords_file <- "tests/test_data/mouse_coords_mm39.csv"
    } else if (genome == "mm10") {
      coords_file <- system.file("extdata", "mouse_coords_mm10.csv", package = "locusPackRat")
      if (coords_file == "") coords_file <- "tests/test_data/mouse_coords_mm10.csv"
    }
  } else if (species == "human") {
    if (genome == "hg38") {
      coords_file <- system.file("extdata", "human_coords_hg38.csv", package = "locusPackRat")
      if (coords_file == "") coords_file <- "tests/test_data/human_coords_hg38.csv"
    } else if (genome == "hg19") {
      coords_file <- system.file("extdata", "human_coords_hg19.csv", package = "locusPackRat")
      if (coords_file == "") coords_file <- "tests/test_data/human_coords_hg19.csv"
    }
  }

  # Merge with gene coordinates
  if (!is.null(coords_file) && file.exists(coords_file)) {
    coords_data <- fread(coords_file)
    # Keep only the coordinate columns we need
    coords_cols <- intersect(c("gene_symbol", "chr", "start", "end", "strand", "ensembl_id"), names(coords_data))
    coords_data <- coords_data[, ..coords_cols]

    result <- merge(result, coords_data, by = "gene_symbol", all.x = TRUE)

    # Ensure numeric coordinates
    if ("start" %in% names(result)) result[, start := as.integer(start)]
    if ("end" %in% names(result)) result[, end := as.integer(end)]

    message(sprintf("  Loaded gene coordinates for %d/%d genes", sum(!is.na(result$start)), nrow(result)))
  } else {
    message("  Note: Coordinate file not found, gene-level overlap will use region bounds")
    # Fall back to region coordinates
    result[, chr := region_chr]
    result[, start := region_start]
    result[, end := region_end]
  }

  return(result)
}

#' Summarize at Locus Level
#' @noRd
.summarizeAtLocusLevel <- function(data, config) {
  if (config$mode == "region") {
    # Group by region
    if ("region_id" %in% names(data)) {
      summary <- data[, .(
        n_genes = .N,
        genes = paste(unique(gene_symbol), collapse = ", "),
        chr = chr[1],
        start = min(start, na.rm = TRUE),
        end = max(end, na.rm = TRUE)
      ), by = region_id]
    } else {
      # Group by coordinates
      summary <- data[, .(
        n_genes = .N,
        genes = paste(unique(gene_symbol), collapse = ", ")
      ), by = .(chr, start, end)]
    }
  } else {
    # For gene mode, group by chromosome
    summary <- data[, .(
      n_genes = .N,
      genes = paste(unique(gene_symbol), collapse = ", "),
      start = min(start, na.rm = TRUE),
      end = max(end, na.rm = TRUE)
    ), by = chr]
  }

  return(summary)
}

#' Sort Data by region_id then gene_symbol
#' @noRd
.sortByRegionAndGene <- function(dt) {
  sort_cols <- c()
  if ("region_id" %in% names(dt)) sort_cols <- c(sort_cols, "region_id")
  if ("gene_symbol" %in% names(dt)) sort_cols <- c(sort_cols, "gene_symbol")
  if (length(sort_cols) > 0) {
    data.table::setorderv(dt, sort_cols, na.last = TRUE)
  }
  return(dt)
}

#' Detect Column Format Based on Name and Data
#' @noRd
.detectColumnFormat <- function(col_name, col_data) {
  name_lower <- tolower(col_name)

  # P-value patterns
  if (grepl("pval|p\\.value|p_value|fdr|qvalue|padj|adj\\.p", name_lower)) {
    return("scientific")
  }
  # Fold change patterns
  if (grepl("log2|lfc|fold|fc$|_fc_", name_lower)) {
    return("decimal2")
  }
  # Coordinate patterns
  if (grepl("^(start|end)$|_start$|_end$|position|^pos$|_pos$", name_lower)) {
    return("thousands")
  }
  # Percentage patterns
  if (grepl("percent|pct|proportion", name_lower)) {
    return("percent")
  }

  return("general")
}

#' Calculate Smart Column Width
#' @noRd
.calculateColumnWidth <- function(col_name, col_data, min_width = 10, max_width = 50) {
  # Get header width
  header_width <- nchar(col_name)

  # Sample data width (first 100 non-NA values)
  if (is.character(col_data) || is.factor(col_data)) {
    sample_data <- head(col_data[!is.na(col_data)], 100)
    if (length(sample_data) > 0) {
      data_width <- max(nchar(as.character(sample_data)), na.rm = TRUE)
    } else {
      data_width <- 0
    }
  } else if (is.numeric(col_data)) {
    # Numbers need less width
    data_width <- 12
  } else {
    data_width <- 10
  }

  # Return constrained width
  width <- max(header_width, data_width) + 2  # Add padding
  return(min(max(width, min_width), max_width))
}

#' Create Summary Sheet
#' @noRd
.createSummarySheet <- function(wb, sheets_info, filter_expr = NULL,
                                phenotypes_data = NULL, diseases_data = NULL) {
  openxlsx::addWorksheet(wb, "Summary", gridLines = FALSE)

  # Get package version
 pkg_version <- tryCatch(
    as.character(utils::packageVersion("locusPackRat")),
    error = function(e) "dev"
  )

  # Build summary content
  summary_lines <- c(
    "locusPackRat Output Summary",
    paste(rep("=", 30), collapse = ""),
    "",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("Package version:", pkg_version),
    "",
    "Sheets in this workbook:",
    paste(rep("-", 25), collapse = "")
  )

  # Add sheet info
  for (sheet_name in names(sheets_info)) {
    n_rows <- sheets_info[[sheet_name]]
    summary_lines <- c(summary_lines, sprintf("  %s: %d rows", sheet_name, n_rows))
  }

  summary_lines <- c(summary_lines, "", "Formatting Legend:", paste(rep("-", 25), collapse = ""))
  summary_lines <- c(summary_lines,
    "  Yellow highlight: genes matching highlight criteria",
    "  Striped rows: alternating for readability",
    "  Scientific notation: p-values and FDR columns",
    "  Decimal format: fold change columns",
    "  Thousands separator: genomic coordinates"
  )

  if (!is.null(filter_expr)) {
    summary_lines <- c(summary_lines, "", "Filter applied:", paste(" ", filter_expr))
  }

  # Add top 15 phenotypes section
  if (!is.null(phenotypes_data) && nrow(phenotypes_data) > 0) {
    pheno_col <- intersect(c("phenotype", "mm_phenotype"), names(phenotypes_data))[1]
    if (!is.na(pheno_col)) {
      pheno_counts <- sort(table(phenotypes_data[[pheno_col]]), decreasing = TRUE)
      top_pheno <- head(pheno_counts, 15)
      summary_lines <- c(summary_lines, "", "Top 15 Phenotypes (by frequency):",
                         paste(rep("-", 35), collapse = ""))
      for (i in seq_along(top_pheno)) {
        summary_lines <- c(summary_lines,
                           sprintf("  %d. %s (%d)", i, names(top_pheno)[i], top_pheno[i]))
      }
    }
  }

  # Add top 15 diseases section
  if (!is.null(diseases_data) && nrow(diseases_data) > 0) {
    disease_col <- intersect(c("disease_name", "otd_disease_name"), names(diseases_data))[1]
    if (!is.na(disease_col)) {
      disease_counts <- sort(table(diseases_data[[disease_col]]), decreasing = TRUE)
      top_diseases <- head(disease_counts, 15)
      summary_lines <- c(summary_lines, "", "Top 15 Diseases (by frequency):",
                         paste(rep("-", 35), collapse = ""))
      for (i in seq_along(top_diseases)) {
        summary_lines <- c(summary_lines,
                           sprintf("  %d. %s (%d)", i, names(top_diseases)[i], top_diseases[i]))
      }
    }
  }

  # Write summary content
  summary_df <- data.frame(Summary = summary_lines, stringsAsFactors = FALSE)
  openxlsx::writeData(wb, sheet = "Summary", x = summary_df, colNames = FALSE)

  # Style the title
  title_style <- openxlsx::createStyle(
    fontSize = 14,
    textDecoration = "bold",
    fontColour = "#1F4E79"
  )
  openxlsx::addStyle(wb, sheet = "Summary", style = title_style, rows = 1, cols = 1)

  # Set column width for readability
  openxlsx::setColWidths(wb, sheet = "Summary", cols = 1, widths = 60)
}

#' Create Phenotypes Detail Sheet from MouseMine Data
#' @noRd
.createPhenotypesDetailSheet <- function(wb, data, header_style, merge_repeated = FALSE) {
  # Check for mouseMine phenotype columns
  phenotype_cols <- c("mm_phenotype", "phenotype")
  pheno_col <- intersect(phenotype_cols, names(data))[1]

  if (is.na(pheno_col) || !pheno_col %in% names(data)) {
    return(wb)
  }

  # Only proceed if we have semicolon-separated values
  has_multi <- any(grepl(";", data[[pheno_col]], fixed = TRUE), na.rm = TRUE)
  if (!has_multi) {
    return(wb)
  }

  message("  Creating Phenotypes_Detail sheet...")

  # Identify related columns to explode together
  # Include from_region if present (region mode)
  possible_cols <- c("gene_symbol", "from_region", "mm_mp_id", "mp_id",
                     "mm_phenotype", "phenotype", "mm_pubmed_id", "pubmed_id",
                     "mm_description", "description")
  cols_present <- intersect(possible_cols, names(data))

  if (length(cols_present) < 2) {
    return(wb)
  }

  # Extract relevant columns
  detail_data <- data[, cols_present, with = FALSE]
  detail_data <- unique(detail_data)

  # Explode semicolon-separated values
  result_list <- list()

  for (i in seq_len(nrow(detail_data))) {
    row <- detail_data[i]

    # Skip if phenotype is NA or empty
    if (is.na(row[[pheno_col]]) || row[[pheno_col]] == "") {
      next
    }

    # Split phenotype column (primary)
    pheno_vals <- trimws(unlist(strsplit(row[[pheno_col]], ";", fixed = TRUE)))
    pheno_vals <- pheno_vals[pheno_vals != ""]

    n_vals <- length(pheno_vals)
    if (n_vals == 0) next

    # Create expanded rows
    expanded <- data.table(gene_symbol = rep(row$gene_symbol, n_vals))

    # Add from_region if present (region mode)
    if ("from_region" %in% cols_present && !is.na(row$from_region)) {
      expanded[["from_region"]] <- rep(row$from_region, n_vals)
    }

    # Add phenotype
    expanded[[pheno_col]] <- pheno_vals

    # Handle mp_id column (may have fewer values)
    mp_col <- intersect(c("mm_mp_id", "mp_id"), cols_present)[1]
    if (!is.na(mp_col)) {
      mp_vals <- if (!is.na(row[[mp_col]])) {
        trimws(unlist(strsplit(row[[mp_col]], ";", fixed = TRUE)))
      } else {
        NA_character_
      }
      # Recycle to match length
      expanded[[mp_col]] <- rep_len(mp_vals, n_vals)
    }

    # Handle pubmed_id column
    pub_col <- intersect(c("mm_pubmed_id", "pubmed_id"), cols_present)[1]
    if (!is.na(pub_col)) {
      pub_vals <- if (!is.na(row[[pub_col]])) {
        trimws(unlist(strsplit(as.character(row[[pub_col]]), ";", fixed = TRUE)))
      } else {
        NA_character_
      }
      expanded[[pub_col]] <- rep_len(pub_vals, n_vals)
    }

    # Handle description column
    desc_col <- intersect(c("mm_description", "description"), cols_present)[1]
    if (!is.na(desc_col)) {
      desc_vals <- if (!is.na(row[[desc_col]])) {
        trimws(unlist(strsplit(row[[desc_col]], ";", fixed = TRUE)))
      } else {
        NA_character_
      }
      expanded[[desc_col]] <- rep_len(desc_vals, n_vals)
    }

    result_list[[i]] <- expanded
  }

  if (length(result_list) == 0) {
    return(wb)
  }

  detail_result <- rbindlist(result_list, fill = TRUE)
  detail_result <- unique(detail_result)

  # Filter: remove rows where only gene_symbol has a value
  non_gene_cols <- setdiff(names(detail_result), "gene_symbol")
  if (length(non_gene_cols) > 0) {
    has_data <- rowSums(!is.na(detail_result[, non_gene_cols, with = FALSE]) &
                        detail_result[, non_gene_cols, with = FALSE] != "") > 0
    detail_result <- detail_result[has_data]
  }

  if (nrow(detail_result) == 0) {
    return(wb)
  }

  # Sort by region_id then gene_symbol for consistent ordering
  detail_result <- .sortByRegionAndGene(detail_result)

  # Add worksheet (use writeData instead of writeDataTable if merging)
  openxlsx::addWorksheet(wb, "Phenotypes_Detail")

  # Style for gene_symbol column: bold, centered vert+horiz
  gene_symbol_style <- openxlsx::createStyle(
    textDecoration = "bold",
    halign = "center",
    valign = "center"
  )

  # Border styles for group boundaries
  top_border_style <- openxlsx::createStyle(border = "top", borderStyle = "thin")
  bottom_border_style <- openxlsx::createStyle(border = "bottom", borderStyle = "thin")

  if (merge_repeated) {
    # Write without table format for merging
    openxlsx::writeData(wb, sheet = "Phenotypes_Detail", x = detail_result,
                        headerStyle = header_style)

    # Find gene_symbol column index
    gene_col_idx <- which(names(detail_result) == "gene_symbol")

    if (length(gene_col_idx) > 0) {
      # Merge consecutive cells with same gene_symbol
      merge_style <- openxlsx::createStyle(
        textDecoration = "bold",
        valign = "center",
        halign = "center"
      )
      genes <- detail_result$gene_symbol
      start_row <- 2
      current_gene <- genes[1]

      for (row_idx in 2:length(genes)) {
        if (genes[row_idx] != current_gene || row_idx == length(genes)) {
          end_row <- if (genes[row_idx] != current_gene) row_idx else row_idx + 1
          if (end_row - start_row > 1) {
            openxlsx::mergeCells(wb, sheet = "Phenotypes_Detail",
                                 cols = gene_col_idx, rows = start_row:end_row)
            openxlsx::addStyle(wb, sheet = "Phenotypes_Detail", style = merge_style,
                               rows = start_row, cols = gene_col_idx)
          }
          # Add borders around gene groups
          openxlsx::addStyle(wb, sheet = "Phenotypes_Detail", style = top_border_style,
                             rows = start_row, cols = 1:ncol(detail_result), stack = TRUE)
          openxlsx::addStyle(wb, sheet = "Phenotypes_Detail", style = bottom_border_style,
                             rows = end_row, cols = 1:ncol(detail_result), stack = TRUE)
          start_row <- row_idx + 1
          current_gene <- genes[row_idx]
        }
      }
    }

    # Add filter manually
    openxlsx::addFilter(wb, sheet = "Phenotypes_Detail", rows = 1, cols = 1:ncol(detail_result))
  } else {
    openxlsx::writeDataTable(wb, sheet = "Phenotypes_Detail", x = detail_result,
                             tableName = "Phenotypes_Detail", withFilter = TRUE,
                             tableStyle = "TableStyleLight9")

    # Find gene_symbol column index and apply styling + borders
    gene_col_idx <- which(names(detail_result) == "gene_symbol")

    if (length(gene_col_idx) > 0 && nrow(detail_result) > 0) {
      # Apply bold centered style to all gene_symbol cells
      openxlsx::addStyle(wb, sheet = "Phenotypes_Detail", style = gene_symbol_style,
                         rows = 2:(nrow(detail_result) + 1), cols = gene_col_idx, stack = TRUE)

      # Add borders around gene groups
      genes <- detail_result$gene_symbol
      group_start <- 2  # First data row (1 is header)

      for (row_idx in seq_along(genes)) {
        is_last <- (row_idx == length(genes))
        is_boundary <- is_last || (genes[row_idx] != genes[row_idx + 1])

        if (is_boundary) {
          data_row <- row_idx + 1  # +1 for header row
          # Top border on first row of group
          openxlsx::addStyle(wb, sheet = "Phenotypes_Detail", style = top_border_style,
                             rows = group_start, cols = 1:ncol(detail_result), stack = TRUE)
          # Bottom border on last row of group
          openxlsx::addStyle(wb, sheet = "Phenotypes_Detail", style = bottom_border_style,
                             rows = data_row, cols = 1:ncol(detail_result), stack = TRUE)
          group_start <- data_row + 1
        }
      }
    }
  }

  openxlsx::addStyle(wb, sheet = "Phenotypes_Detail", style = header_style,
                     rows = 1, cols = 1:ncol(detail_result), stack = TRUE)
  openxlsx::setColWidths(wb, sheet = "Phenotypes_Detail", cols = 1:ncol(detail_result), widths = "auto")
  openxlsx::freezePane(wb, sheet = "Phenotypes_Detail", firstActiveRow = 2, firstActiveCol = 2)

  return(wb)
}

#' Create Diseases Detail Sheet from OpenTargets Data
#' @noRd
.createDiseasesDetailSheet <- function(wb, data, header_style) {
  # Check for OpenTargets disease columns (with or without prefix)
  disease_cols <- c("otd_disease_name", "disease_name")
  disease_col <- intersect(disease_cols, names(data))[1]

  if (is.na(disease_col) || !disease_col %in% names(data)) {
    return(wb)
  }

  # Only proceed if we have semicolon-separated values
  has_multi <- any(grepl(";", data[[disease_col]], fixed = TRUE), na.rm = TRUE)
  if (!has_multi) {
    return(wb)
  }

  message("  Creating Diseases_Detail sheet...")

  # Identify related columns to explode together
  # Include from_region if present (region mode)
  possible_cols <- c("gene_symbol", "from_region", "otd_disease_id", "disease_id",
                     "otd_disease_name", "disease_name", "otd_score", "score")
  cols_present <- intersect(possible_cols, names(data))

  if (length(cols_present) < 2) {
    return(wb)
  }

  # Extract relevant columns
  detail_data <- data[, cols_present, with = FALSE]
  detail_data <- unique(detail_data)

  # Explode semicolon-separated values
  result_list <- list()

  for (i in seq_len(nrow(detail_data))) {
    row <- detail_data[i]

    # Skip if disease_name is NA or empty
    if (is.na(row[[disease_col]]) || row[[disease_col]] == "") {
      next
    }

    # Split disease_name column (primary)
    disease_vals <- trimws(unlist(strsplit(row[[disease_col]], ";", fixed = TRUE)))
    disease_vals <- disease_vals[disease_vals != ""]

    n_vals <- length(disease_vals)
    if (n_vals == 0) next

    # Create expanded rows
    expanded <- data.table(gene_symbol = rep(row$gene_symbol, n_vals))

    # Add from_region if present (region mode)
    if ("from_region" %in% cols_present && !is.na(row$from_region)) {
      expanded[["from_region"]] <- rep(row$from_region, n_vals)
    }

    # Add disease_name
    expanded[[disease_col]] <- disease_vals

    # Handle disease_id column (may have fewer values)
    id_col <- intersect(c("otd_disease_id", "disease_id"), cols_present)[1]
    if (!is.na(id_col)) {
      id_vals <- if (!is.na(row[[id_col]])) {
        trimws(unlist(strsplit(row[[id_col]], ";", fixed = TRUE)))
      } else {
        NA_character_
      }
      # Recycle to match length
      expanded[[id_col]] <- rep_len(id_vals, n_vals)
    }

    # Handle score column
    score_col <- intersect(c("otd_score", "score"), cols_present)[1]
    if (!is.na(score_col)) {
      score_vals <- if (!is.na(row[[score_col]])) {
        trimws(unlist(strsplit(as.character(row[[score_col]]), ";", fixed = TRUE)))
      } else {
        NA_character_
      }
      # Try to convert to numeric
      score_numeric <- suppressWarnings(as.numeric(score_vals))
      expanded[[score_col]] <- rep_len(score_numeric, n_vals)
    }

    result_list[[i]] <- expanded
  }

  if (length(result_list) == 0) {
    return(wb)
  }

  detail_result <- rbindlist(result_list, fill = TRUE)
  detail_result <- unique(detail_result)

  # Filter: remove rows where only gene_symbol has a value
  non_gene_cols <- setdiff(names(detail_result), "gene_symbol")
  if (length(non_gene_cols) > 0) {
    has_data <- rowSums(!is.na(detail_result[, non_gene_cols, with = FALSE]) &
                        detail_result[, non_gene_cols, with = FALSE] != "") > 0
    detail_result <- detail_result[has_data]
  }

  if (nrow(detail_result) == 0) {
    return(wb)
  }

  # Sort by region_id then gene_symbol for consistent ordering
  detail_result <- .sortByRegionAndGene(detail_result)

  # Add worksheet
  openxlsx::addWorksheet(wb, "Diseases_Detail")

  # Style for gene_symbol column: bold, centered vert+horiz
  gene_symbol_style <- openxlsx::createStyle(
    textDecoration = "bold",
    halign = "center",
    valign = "center"
  )

  # Border styles for group boundaries
  top_border_style <- openxlsx::createStyle(border = "top", borderStyle = "thin")
  bottom_border_style <- openxlsx::createStyle(border = "bottom", borderStyle = "thin")

  openxlsx::writeDataTable(wb, sheet = "Diseases_Detail", x = detail_result,
                           tableName = "Diseases_Detail", withFilter = TRUE,
                           tableStyle = "TableStyleLight9")

  # Find gene_symbol column index and apply styling + borders
  gene_col_idx <- which(names(detail_result) == "gene_symbol")

  if (length(gene_col_idx) > 0 && nrow(detail_result) > 0) {
    # Apply bold centered style to all gene_symbol cells
    openxlsx::addStyle(wb, sheet = "Diseases_Detail", style = gene_symbol_style,
                       rows = 2:(nrow(detail_result) + 1), cols = gene_col_idx, stack = TRUE)

    # Add borders around gene groups
    genes <- detail_result$gene_symbol
    group_start <- 2  # First data row (1 is header)

    for (row_idx in seq_along(genes)) {
      is_last <- (row_idx == length(genes))
      is_boundary <- is_last || (genes[row_idx] != genes[row_idx + 1])

      if (is_boundary) {
        data_row <- row_idx + 1  # +1 for header row
        # Top border on first row of group
        openxlsx::addStyle(wb, sheet = "Diseases_Detail", style = top_border_style,
                           rows = group_start, cols = 1:ncol(detail_result), stack = TRUE)
        # Bottom border on last row of group
        openxlsx::addStyle(wb, sheet = "Diseases_Detail", style = bottom_border_style,
                           rows = data_row, cols = 1:ncol(detail_result), stack = TRUE)
        group_start <- data_row + 1
      }
    }
  }

  openxlsx::addStyle(wb, sheet = "Diseases_Detail", style = header_style,
                     rows = 1, cols = 1:ncol(detail_result), stack = TRUE)
  openxlsx::setColWidths(wb, sheet = "Diseases_Detail", cols = 1:ncol(detail_result), widths = "auto")
  openxlsx::freezePane(wb, sheet = "Diseases_Detail", firstActiveRow = 2, firstActiveCol = 2)

  return(wb)
}

#' Create Regions Summary Sheet
#' @noRd
.createRegionsSheet <- function(wb, data, header_style) {
  # Check if this is region mode data
  region_cols <- c("region_id", "region_chr", "region_start", "region_end", "genes", "overlap_types")
  present_cols <- intersect(region_cols, names(data))

  if (length(present_cols) < 3) {
    return(list(wb = wb, data = data))
  }

  message("  Creating Regions sheet...")

  # Extract unique region data
  # Include any user columns that might be at region level (e.g., peak_id, max_lod)
  key_cols <- c("region_id")
  key_col <- intersect(key_cols, names(data))[1]

  if (is.na(key_col)) {
    return(list(wb = wb, data = data))
  }

  # Get all columns that have same value per region
  region_level_cols <- c()
  for (col in names(data)) {
    if (col == "gene_symbol") next
    # Check if column is constant within each region
    test <- data[, .(n_unique = data.table::uniqueN(get(col))), by = key_col]
    if (all(test$n_unique <= 1, na.rm = TRUE)) {
      region_level_cols <- c(region_level_cols, col)
    }
  }

  # Create regions summary
  if (length(region_level_cols) > 0) {
    regions_data <- unique(data[, region_level_cols, with = FALSE])

    # Sort alphabetically by region_id for consistent ordering
    if ("region_id" %in% names(regions_data)) {
      data.table::setorderv(regions_data, "region_id", na.last = TRUE)
    }

    # Remove columns where all values are NA
    all_na_cols <- sapply(regions_data, function(x) all(is.na(x)))
    if (any(all_na_cols)) {
      regions_data <- regions_data[, !all_na_cols, with = FALSE]
    }

    # Add worksheet
    openxlsx::addWorksheet(wb, "Regions")
    openxlsx::writeDataTable(wb, sheet = "Regions", x = regions_data,
                             tableName = "Regions", withFilter = TRUE,
                             tableStyle = "TableStyleLight9")
    openxlsx::addStyle(wb, sheet = "Regions", style = header_style,
                       rows = 1, cols = 1:ncol(regions_data), stack = TRUE)

    # Auto-wrap columns with long content (> 50 chars)
    wrap_style <- openxlsx::createStyle(wrapText = TRUE, valign = "top")
    for (col_idx in seq_along(names(regions_data))) {
      col_data <- regions_data[[col_idx]]
      if (is.character(col_data)) {
        max_len <- max(nchar(col_data), na.rm = TRUE)
        if (!is.na(max_len) && max_len > 50) {
          # Apply wrap to all data rows in this column
          openxlsx::addStyle(wb, sheet = "Regions", style = wrap_style,
                             rows = 2:(nrow(regions_data) + 1), cols = col_idx, stack = TRUE)
          # Set column width to 80 for wrapped columns
          openxlsx::setColWidths(wb, sheet = "Regions", cols = col_idx, widths = 80)
        } else {
          openxlsx::setColWidths(wb, sheet = "Regions", cols = col_idx, widths = "auto")
        }
      } else {
        openxlsx::setColWidths(wb, sheet = "Regions", cols = col_idx, widths = "auto")
      }
    }

    openxlsx::freezePane(wb, sheet = "Regions", firstActiveRow = 2, firstActiveCol = 2)
  }

  # Remove redundant region columns from main data
  cols_to_remove <- present_cols
  data_cleaned <- data[, !names(data) %in% cols_to_remove, with = FALSE]

  return(list(wb = wb, data = data_cleaned))
}

#' Create Formatted Excel Workbook
#' @noRd
.createFormattedExcel <- function(data, highlight_genes, highlight_expr,
                                  highlight_color, split_by, split_criteria,
                                  format_numbers = TRUE, format_pvalues = TRUE,
                                  format_coordinates = TRUE, include_summary = TRUE,
                                  header_color = "#4472C4", min_col_width = 10,
                                  filter_expr = NULL, merge_repeated = FALSE) {

  wb <- openxlsx::createWorkbook()

  # Create styles
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fontColour = "#FFFFFF",
    fgFill = header_color,
    halign = "center",
    valign = "center",
    border = "bottom",
    borderStyle = "medium",
    wrapText = TRUE
  )
  highlight_style <- openxlsx::createStyle(fgFill = highlight_color)
  stripe_style <- openxlsx::createStyle(fgFill = "#F0F0F0")

  # Number format styles
  scientific_style <- openxlsx::createStyle(numFmt = "0.00E+00")
  decimal2_style <- openxlsx::createStyle(numFmt = "0.00")
  thousands_style <- openxlsx::createStyle(numFmt = "#,##0")
  percent_style <- openxlsx::createStyle(numFmt = "0.0%")

  # Store original data for detail sheets before cleaning
  original_data <- copy(data)

  # Handle region mode: create Regions sheet and clean data
  regions_result <- .createRegionsSheet(wb, data, header_style)
  wb <- regions_result$wb
  data <- regions_result$data

  # Determine sheets
  sheets <- list()

  if (split_by == "none") {
    sheets[["All Genes"]] <- data
  } else if (split_by == "chromosome") {
    if ("chr" %in% names(data)) {
      sheets <- split(data, data$chr)
      names(sheets) <- paste0("Chr_", names(sheets))
    } else {
      sheets[["All Genes"]] <- data
    }
  } else if (split_by == "criteria" && !is.null(split_criteria)) {
    sheets[["All_Data"]] <- data
    for (tab_name in names(split_criteria)) {
      f_expr <- split_criteria[[tab_name]]
      tryCatch({
        filtered <- data[eval(parse(text = f_expr))]
        sheets[[tab_name]] <- filtered
      }, error = function(e) {
        warning(sprintf("Skipping tab '%s': %s", tab_name, e$message))
      })
    }
  }

  # Track sheet row counts for summary
  sheets_info <- list()

  for (sheet_name in names(sheets)) {
    sheet_data <- sheets[[sheet_name]]
    if (nrow(sheet_data) == 0) next

    # Remove columns where all values are NA
    all_na_cols <- sapply(sheet_data, function(x) all(is.na(x)))
    if (any(all_na_cols)) {
      sheet_data <- sheet_data[, !all_na_cols, with = FALSE]
    }

    # Sort by region_id then gene_symbol for consistent ordering
    sheet_data <- .sortByRegionAndGene(sheet_data)

    sheets_info[[sheet_name]] <- nrow(sheet_data)

    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeDataTable(
      wb,
      sheet = sheet_name,
      x = sheet_data,
      tableName = gsub("[^A-Za-z0-9]", "_", sheet_name),
      withFilter = TRUE,
      tableStyle = "TableStyleLight9"
    )

    # Apply header style
    openxlsx::addStyle(wb, sheet = sheet_name, style = header_style,
                       rows = 1, cols = 1:ncol(sheet_data), stack = TRUE)

    # Apply striped rows (need at least 2 data rows for stripes to make sense)
    if (nrow(sheet_data) >= 2) {
      rows <- seq(3, nrow(sheet_data) + 1, by = 2)
      if (length(rows) > 0) {
        openxlsx::addStyle(wb, sheet = sheet_name, style = stripe_style,
                           rows = rows, cols = 1:ncol(sheet_data), gridExpand = TRUE, stack = TRUE)
      }
    }

    # Apply number formatting per column
    if (format_numbers) {
      for (col_idx in seq_along(names(sheet_data))) {
        col_name <- names(sheet_data)[col_idx]
        col_data <- sheet_data[[col_name]]

        # Only format numeric columns
        if (!is.numeric(col_data)) next

        fmt <- .detectColumnFormat(col_name, col_data)

        style_to_apply <- switch(fmt,
          "scientific" = if (format_pvalues) scientific_style else NULL,
          "decimal2" = decimal2_style,
          "thousands" = if (format_coordinates) thousands_style else NULL,
          "percent" = percent_style,
          NULL
        )

        if (!is.null(style_to_apply)) {
          data_rows <- 2:(nrow(sheet_data) + 1)
          openxlsx::addStyle(wb, sheet = sheet_name, style = style_to_apply,
                             rows = data_rows, cols = col_idx, stack = TRUE)
        }
      }
    }

    # Apply highlighting
    rows_to_hl <- c()
    if (!is.null(highlight_genes) && "gene_symbol" %in% names(sheet_data)) {
      rows_to_hl <- which(sheet_data$gene_symbol %in% highlight_genes) + 1
    }
    if (!is.null(highlight_expr)) {
      tryCatch({
        expr_rows <- which(sheet_data[, eval(parse(text = highlight_expr))]) + 1
        rows_to_hl <- unique(c(rows_to_hl, expr_rows))
      }, error = function(e) NULL)
    }
    if (length(rows_to_hl) > 0) {
      openxlsx::addStyle(wb, sheet = sheet_name, style = highlight_style,
                         rows = rows_to_hl, cols = 1:ncol(sheet_data), gridExpand = TRUE, stack = TRUE)
    }

    # Set smart column widths
    col_widths <- sapply(seq_along(names(sheet_data)), function(i) {
      .calculateColumnWidth(names(sheet_data)[i], sheet_data[[i]], min_width = min_col_width)
    })
    openxlsx::setColWidths(wb, sheet = sheet_name, cols = 1:ncol(sheet_data), widths = col_widths)

    # Freeze header row
    openxlsx::freezePane(wb, sheet = sheet_name, firstActiveRow = 2, firstActiveCol = 2)
  }

  # Add phenotypes detail sheet if mouseMine data present
  wb <- .createPhenotypesDetailSheet(wb, original_data, header_style, merge_repeated)

  # Add diseases detail sheet if OpenTargets data present
  wb <- .createDiseasesDetailSheet(wb, original_data, header_style)

  # Update sheets_info to include Regions, Phenotypes_Detail, and Diseases_Detail if created
  all_sheets <- openxlsx::sheets(wb)
  if ("Regions" %in% all_sheets) {
    # Get row count from Regions sheet
    regions_nrow <- nrow(openxlsx::readWorkbook(wb, sheet = "Regions"))
    sheets_info <- c(list(Regions = regions_nrow), sheets_info)
  }
  phenotypes_data <- NULL
  diseases_data <- NULL
  if ("Phenotypes_Detail" %in% all_sheets) {
    # Get data from Phenotypes_Detail sheet
    phenotypes_data <- data.table::as.data.table(openxlsx::readWorkbook(wb, sheet = "Phenotypes_Detail"))
    sheets_info[["Phenotypes_Detail"]] <- nrow(phenotypes_data)
  }
  if ("Diseases_Detail" %in% all_sheets) {
    # Get data from Diseases_Detail sheet
    diseases_data <- data.table::as.data.table(openxlsx::readWorkbook(wb, sheet = "Diseases_Detail"))
    sheets_info[["Diseases_Detail"]] <- nrow(diseases_data)
  }

  # Add summary sheet (at the beginning)
  if (include_summary && length(sheets_info) > 0) {
    .createSummarySheet(wb, sheets_info, filter_expr, phenotypes_data, diseases_data)
    # Move summary sheet to first position
    openxlsx::worksheetOrder(wb) <- c(length(openxlsx::sheets(wb)), 1:(length(openxlsx::sheets(wb)) - 1))
  }

  return(wb)
}

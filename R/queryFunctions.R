# Suppress R CMD check NOTEs for data.table NSE variables
utils::globalVariables(c(
  "disease_id", "disease_name", "therapeutic_areas", "l2g_score",
  "disease", "human_ensembl_id"
))

#' @noRd
.getProjectGenes <- function(packrat_dir) {
  config <- jsonlite::read_json(file.path(packrat_dir, "config.json"))
  
  if (config$mode == "gene") {
    dt <- data.table::fread(file.path(packrat_dir, "input", "genes.csv"))
    return(unique(dt$gene_symbol))
  } else {
    # In region mode, genes are nested in the regions file
    dt <- data.table::fread(file.path(packrat_dir, "input", "regions.csv"))
    # Expand comma-separated genes
    genes <- unlist(strsplit(dt$genes, ", "))
    return(unique(trimws(genes)))
  }
}

#' Get Project Regions for QTL Query
#' @noRd
.getProjectRegions <- function(packrat_dir, config) {
  if (config$mode == "gene") {
    # Gene mode: look up coordinates from genes.csv
    genes_dt <- data.table::fread(file.path(packrat_dir, "input", "genes.csv"))
    coord_file <- system.file("extdata",
                              paste0(config$species, "_coords_", config$genome, ".csv"),
                              package = "locusPackRat")
    if (coord_file == "") {
      coord_file <- paste0("inst/extdata/", config$species, "_coords_", config$genome, ".csv")
    }
    if (!file.exists(coord_file)) {
      stop("Coordinate file not found: ", coord_file)
    }
    coords <- data.table::fread(coord_file)
    result <- merge(genes_dt[, .(gene_symbol)],
                    coords[, .(gene_symbol, chr, start, end)],
                    by = "gene_symbol", all.x = TRUE)
    result <- result[!is.na(chr)]
    # Add region_id based on gene symbol for gene mode
    result[, region_id := gene_symbol]
    return(result)
  } else {
    # Region mode: read regions directly
    regions_dt <- data.table::fread(file.path(packrat_dir, "input", "regions.csv"))
    return(regions_dt[, .(region_id, chr, start, end)])
  }
}

#' Query MouseMine for Phenotypes
#'
#' Queries the MouseMine database for gene phenotype annotations using the
#' MouseMine REST API directly (no InterMineR dependency required).
#'
#' @param project_dir Character path to the project directory
#' @param limit Integer (optional). If provided, limits the query to the first N genes. Useful for testing.
#' @param chunk_size Integer. Number of genes to query per API call. Default 200.
#'
#' @return Invisible data.table of results
#'
#' @examples
#' \dontrun{
#' # Query MouseMine for all genes in project
#' queryMouseMine(project_dir = "my_analysis")
#'
#' # Test with limited genes
#' queryMouseMine(project_dir = "my_analysis", limit = 10)
#' }
#'
#' @export
queryMouseMine <- function(project_dir = ".", limit = NULL, chunk_size = 200) {

  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) stop("Project not found.")

  # 1. Get Genes

  message("Retrieving gene list from project...")
  gene_syms <- .getProjectGenes(packrat_dir)

  # Limit for testing
  if (!is.null(limit)) {
    message(sprintf("TEST MODE: Limiting query to first %d genes", limit))
    gene_syms <- head(gene_syms, limit)
  }

  # 2. MouseMine REST API endpoint
  base_url <- "http://www.mousemine.org/mousemine/service/template/results"

  # 3. Batch Processing
  chunks <- split(gene_syms, ceiling(seq_along(gene_syms) / chunk_size))
  results_list <- list()

  message(sprintf("Querying %d genes in %d batches...", length(gene_syms), length(chunks)))

  for (i in seq_along(chunks)) {
    syms <- chunks[[i]]
    message(sprintf("  Processing batch %d/%d...", i, length(chunks)))

    tryCatch({
      response <- httr::GET(
        base_url,
        query = list(
          name = "_Feature_Phenotype",
          constraint1 = "OntologyAnnotation.subject",
          op1 = "LOOKUP",
          value1 = paste(syms, collapse = ","),
          format = "json"
        )
      )

      if (httr::status_code(response) == 200) {
        content <- httr::content(response, "text", encoding = "UTF-8")
        parsed <- jsonlite::fromJSON(content, flatten = TRUE)

        # Extract results from the JSON structure
        if (!is.null(parsed$results) && length(parsed$results) > 0) {
          # Convert to data.table with column names from columnHeaders
          res_dt <- data.table::as.data.table(parsed$results)
          if (!is.null(parsed$columnHeaders)) {
            col_names <- parsed$columnHeaders
            if (length(col_names) == ncol(res_dt)) {
              names(res_dt) <- col_names
            }
          }
          if (nrow(res_dt) > 0) results_list[[i]] <- res_dt
        }
      } else {
        warning(sprintf("Batch %d failed with status %d", i, httr::status_code(response)))
      }
    }, error = function(e) warning("Batch failed: ", e$message))

    Sys.sleep(0.2) # Rate limiting
  }

  final_dt <- data.table::rbindlist(results_list, fill = TRUE)

  # 4. Clean & Save
  if (nrow(final_dt) > 0) {
    # Rename columns to be more user-friendly
    # API returns human-readable names like "Ontology Annotation > Subject . Symbol"
    col_map <- c(
      "Ontology Annotation > Subject . Symbol" = "gene_symbol",
      "Ontology Annotation > Subject > Primary Identifier" = "mgi_id",
      "Ontology Annotation > Term Name" = "phenotype",
      "Ontology Annotation > Ontology Term . Identifier" = "mp_id",
      "Ontology Annotation > Evidence > Publications > PubMed ID" = "pubmed_id",
      "Ontology Annotation > Evidence > Comments > Description" = "description"
    )
    for (old_name in names(col_map)) {
      if (old_name %in% names(final_dt)) {
        data.table::setnames(final_dt, old_name, col_map[[old_name]])
      }
    }
    # Drop less useful columns
    drop_cols <- c("Ontology Annotation > Subject > Type",
                   "Ontology Annotation > Evidence > Comments > Type")
    final_dt[, (intersect(drop_cols, names(final_dt))) := NULL]

    addRatTable(
      data = final_dt,
      table_name = "mouse_phenotypes",
      link_type = "gene",
      link_by = "gene_symbol",
      abbreviation = "mm",
      project_dir = project_dir
    )
    message("Success: 'mouse_phenotypes' added to project.")
  } else {
    warning("No phenotypes found for these genes.")
  }

  invisible(final_dt)
}

#' Query Open Targets Platform
#'
#' Fetches comprehensive gene annotation data from the Open Targets GraphQL API
#' including disease associations, genetic constraints, drug tractability, expression
#' profiles, protein interactions, and more.
#'
#' @param project_dir Character path to the project directory
#' @param data_types Character vector specifying which data types to retrieve.
#'   Available options: \code{"diseases"}, \code{"constraints"}, \code{"tractability"},
#'   \code{"expression"}, \code{"mouse_phenotypes"}, \code{"homologues"},
#'   \code{"interactions"}, \code{"known_drugs"}, \code{"depmap"}, \code{"pharmacogenomics"}.
#'   Default retrieves diseases, constraints, and tractability for backward compatibility.
#' @param disease_limit Integer. Maximum number of disease associations per gene.
#'   Default 500. Set to NULL for no limit.
#' @param limit Integer (optional). Limit query to first N genes for testing.
#' @param chunk_size Integer. Number of genes per API call (default 50).
#' @param verbose Logical. Print progress messages. Default TRUE.
#' @param max_retries Integer. Maximum retry attempts for failed requests. Default 3.
#'
#' @return Invisible list of data.tables for each requested data type
#'
#' @details
#' The function queries the Open Targets Platform API for comprehensive gene annotations.
#' For mouse projects, genes are first mapped to human orthologs before querying.
#'
#' **Output tables by data type:**
#' \itemize{
#'   \item \code{diseases} -> \code{ot_diseases} (otd): Disease associations with scores
#'   \item \code{constraints} -> \code{ot_constraints} (otc): GnomAD genetic constraint metrics
#'   \item \code{tractability} -> \code{ot_tractability} (ott): Drug tractability assessments
#'   \item \code{expression} -> \code{ot_expression} (ote): RNA/protein expression by tissue
#'   \item \code{mouse_phenotypes} -> \code{ot_mouse_phenotypes} (otmp): Mouse phenotypes from OT
#'   \item \code{homologues} -> \code{ot_homologues} (oth): Cross-species orthologs
#'   \item \code{interactions} -> \code{ot_interactions} (oti): Protein-protein interactions
#'   \item \code{known_drugs} -> \code{ot_known_drugs} (otkd): Clinical drug precedence
#'   \item \code{depmap} -> \code{ot_depmap} (otdm): DepMap CRISPR essentiality data
#'   \item \code{pharmacogenomics} -> \code{ot_pharmacogenomics} (otpg): Drug response variants
#' }
#'
#' @examples
#' \dontrun{
#' # Query Open Targets for all genes in project (default: diseases, constraints, tractability)
#' queryOpenTargets(project_dir = "my_analysis")
#'
#' # Query with additional data types
#' queryOpenTargets(project_dir = "my_analysis",
#'                  data_types = c("diseases", "expression", "interactions"))
#'
#' # Query all available data types
#' queryOpenTargets(project_dir = "my_analysis",
#'                  data_types = c("diseases", "constraints", "tractability",
#'                                 "expression", "mouse_phenotypes", "homologues",
#'                                 "interactions", "known_drugs", "depmap",
#'                                 "pharmacogenomics"))
#'
#' # Test with limited genes and more disease associations
#' queryOpenTargets(project_dir = "my_analysis", limit = 5, disease_limit = 1000)
#' }
#'
#' @importFrom stats na.omit
#' @export
queryOpenTargets <- function(project_dir = ".",
                             data_types = c("diseases", "constraints", "tractability"),
                             disease_limit = 500,
                             limit = NULL,
                             chunk_size = 50,
                             verbose = TRUE,
                             max_retries = 3) {

  # 1. Setup and Checks
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for this function.")
  }

  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) stop("Project not found.")

  # Validate data_types
  valid_types <- c("diseases", "constraints", "tractability", "expression",
                   "mouse_phenotypes", "homologues", "interactions",
                   "known_drugs", "depmap", "pharmacogenomics")
  invalid <- setdiff(data_types, valid_types)
  if (length(invalid) > 0) {
    stop("Invalid data_types: ", paste(invalid, collapse = ", "),
         "\nValid options: ", paste(valid_types, collapse = ", "))
  }

  # 2. Read Config
  config_path <- file.path(packrat_dir, "config.json")
  if (!file.exists(config_path)) stop("config.json not found in project.")
  config <- jsonlite::read_json(config_path)

  if (verbose) message(sprintf("Project Config: %s - %s", config$species, config$genome))

  # 3. Get Genes and Resolve Human IDs
  gene_syms <- .getProjectGenes(packrat_dir)
  id_map <- .resolveHumanIds(gene_syms, config)

  query_ids <- unique(na.omit(id_map$human_ensembl_id))

  if (length(query_ids) == 0) {
    warning("No valid Human Ensembl IDs found for these genes. Skipping query.")
    return(NULL)
  }

  if (!is.null(limit)) {
    if (verbose) message(sprintf("TEST MODE: Limiting to first %d genes", limit))
    query_ids <- head(query_ids, limit)
  }

  # 4. Build GraphQL Query based on requested data types
  build_query <- function(ids, types, dis_limit) {
    blocks <- vapply(ids, function(id) {
      alias <- paste0("g_", gsub("-", "_", id))

      # Build field list based on requested types
      fields <- c("id", "approvedSymbol")

      if ("constraints" %in% types) {
        fields <- c(fields, "geneticConstraint { constraintType exp obs score oe oeLower oeUpper }")
      }
      if ("tractability" %in% types) {
        fields <- c(fields, "tractability { label modality value }")
      }
      if ("diseases" %in% types) {
        dis_size <- if (is.null(dis_limit)) 10000 else dis_limit
        fields <- c(fields, sprintf(
          "associatedDiseases(page: { size: %d, index: 0 }) { rows { disease { id name therapeuticAreas { id name } } score datatypeScores { id score } } }",
          dis_size
        ))
      }
      if ("expression" %in% types) {
        fields <- c(fields, "expressions { tissue { id label anatomicalSystems organs } rna { level value unit zscore } protein { level reliability cellType { name level reliability } } }")
      }
      if ("mouse_phenotypes" %in% types) {
        fields <- c(fields, "mousePhenotypes { modelPhenotypeLabel modelPhenotypeId modelPhenotypeClasses { id label } targetFromSourceId targetInModel targetInModelMgiId biologicalModels { id allelicComposition geneticBackground } }")
      }
      if ("homologues" %in% types) {
        fields <- c(fields, "homologues { speciesId speciesName homologyType queryPercentageIdentity targetPercentageIdentity targetGeneId targetGeneSymbol isHighConfidence }")
      }
      if ("interactions" %in% types) {
        fields <- c(fields, "interactions { count rows { intA intABiologicalRole intB intBBiologicalRole score sourceDatabase } }")
      }
      if ("known_drugs" %in% types) {
        fields <- c(fields, "knownDrugs { count rows { drugId drug { id name drugType maximumClinicalTrialPhase isApproved } mechanismOfAction disease { id name } phase status } }")
      }
      if ("depmap" %in% types) {
        fields <- c(fields, "depMapEssentiality { screens { depmapId cellLineName diseaseFromSource mutation geneEffect expression } }")
      }
      if ("pharmacogenomics" %in% types) {
        fields <- c(fields, "pharmacogenomics { variantRsId variantFunctionalConsequence { id label } genotypeId drugs { drugId drugFromSource } phenotypeText genotypeAnnotationText phenotypeFromSourceId evidenceLevel literature datasourceId studyId }")
      }

      sprintf('%s: target(ensemblId: "%s") {\n  %s\n}', alias, id, paste(fields, collapse = "\n  "))
    }, FUN.VALUE = character(1))

    paste0("query multiTarget {\n", paste(blocks, collapse = "\n"), "\n}")
  }

  # 5. Batch Processing with retry logic
  # Reduce chunk size for heavy data types to avoid API timeouts
  if (any(c("mouse_phenotypes", "known_drugs") %in% data_types)) {
    chunk_size <- min(chunk_size, 5)
  }
  batches <- split(query_ids, ceiling(seq_along(query_ids) / chunk_size))
  endpoint <- "https://api.platform.opentargets.org/api/v4/graphql"

  # Initialize result lists for each data type
  results <- list(
    diseases = list(), constraints = list(), tractability = list(),
    expression = list(), mouse_phenotypes = list(), homologues = list(),
    interactions = list(), known_drugs = list(), depmap = list(),
    pharmacogenomics = list()
  )

  if (verbose) message(sprintf("Querying Open Targets for %d genes in %d batches...", length(query_ids), length(batches)))
  if (verbose) message(sprintf("Data types: %s", paste(data_types, collapse = ", ")))

  for (i in seq_along(batches)) {
    batch_ids <- batches[[i]]
    if (verbose) message(sprintf("  Processing batch %d/%d...", i, length(batches)))

    # Retry logic
    success <- FALSE
    for (attempt in seq_len(max_retries)) {
      tryCatch({
        response <- httr::POST(
          endpoint,
          body = list(query = build_query(batch_ids, data_types, disease_limit)),
          encode = "json",
          httr::timeout(120)
        )
        httr::stop_for_status(response)

        content <- httr::content(response, "text", encoding = "UTF-8")
        data_json <- jsonlite::fromJSON(content, flatten = FALSE)$data

        for (alias in names(data_json)) {
          g_dat <- data_json[[alias]]
          if (is.null(g_dat)) next

          h_id <- g_dat$id

          # Parse each requested data type
          results <- .parseOtResponse(g_dat, h_id, data_types, results)
        }

        success <- TRUE
        break
      }, error = function(e) {
        if (attempt < max_retries) {
          if (verbose) message(sprintf("    Attempt %d failed, retrying...", attempt))
          Sys.sleep(2^attempt)  # Exponential backoff
        } else {
          warning("Batch ", i, " failed after ", max_retries, " attempts: ", e$message)
        }
      })
    }

    Sys.sleep(0.5)
  }

  # 6. Save Results
  abbrev_map <- c(
    diseases = "otd", constraints = "otc", tractability = "ott",
    expression = "ote", mouse_phenotypes = "otmp", homologues = "oth",
    interactions = "oti", known_drugs = "otkd", depmap = "otdm",
    pharmacogenomics = "otpg"
  )

  saved_tables <- list()
  for (dtype in data_types) {
    if (length(results[[dtype]]) > 0) {
      dt <- data.table::rbindlist(results[[dtype]], fill = TRUE, use.names = TRUE)

      # Merge back original symbols (e.g. Mouse Symbols)
      dt <- merge(dt, id_map, by = "human_ensembl_id", all.x = TRUE, allow.cartesian = TRUE)

      # Clean up nested columns if they persist
      list_cols <- names(dt)[sapply(dt, is.list)]
      for (lcol in list_cols) {
        dt[, (lcol) := NULL]
      }

      # Reorder to put symbol first
      if ("gene_symbol" %in% names(dt)) {
        col_order <- c("gene_symbol", setdiff(names(dt), "gene_symbol"))
        data.table::setcolorder(dt, col_order)
      }

      table_name <- paste0("ot_", dtype)
      addRatTable(
        data = dt,
        table_name = table_name,
        link_type = "gene",
        link_by = "gene_symbol",
        abbreviation = abbrev_map[[dtype]],
        project_dir = project_dir
      )

      saved_tables[[dtype]] <- dt
      if (verbose) message(sprintf("  Saved %s with %d rows", table_name, nrow(dt)))
    }
  }

  invisible(saved_tables)
}


#' Resolve Gene Symbols to Human Ensembl IDs
#' @noRd
.resolveHumanIds <- function(gene_syms, config) {
  id_map <- data.table::data.table(gene_symbol = gene_syms)

  if (config$species == "human") {
    coords_file <- system.file("extdata", "human_coords_hg38.csv", package = "locusPackRat")
    if (coords_file == "") coords_file <- "inst/extdata/human_coords_hg38.csv"
    if (!file.exists(coords_file)) stop("Internal 'human_coords_hg38.csv' not found.")

    ref_coords <- data.table::fread(coords_file)
    id_map <- merge(id_map, ref_coords[, .(gene_symbol, human_ensembl_id = ensembl_id)],
                    by = "gene_symbol", all.x = TRUE)

  } else if (config$species == "mouse") {
    coords_file <- system.file("extdata", "mouse_coords_mm39.csv", package = "locusPackRat")
    if (coords_file == "") coords_file <- "inst/extdata/mouse_coords_mm39.csv"
    if (!file.exists(coords_file)) stop("Internal 'mouse_coords_mm39.csv' not found.")

    ref_coords <- data.table::fread(coords_file)

    ortho_file <- system.file("extdata", "orthology.csv", package = "locusPackRat")
    if (ortho_file == "") ortho_file <- "inst/extdata/orthology.csv"
    if (!file.exists(ortho_file)) stop("Internal 'orthology.csv' not found.")

    ortho_dat <- data.table::fread(ortho_file)

    step1 <- merge(id_map, ref_coords[, .(gene_symbol, mouse_ensembl_id = ensembl_id)],
                   by = "gene_symbol", all.x = TRUE)
    final_map <- merge(step1, ortho_dat[, .(mouse_ensembl_id, human_ensembl_id)],
                       by = "mouse_ensembl_id", all.x = TRUE)

    id_map <- final_map[, .(gene_symbol, human_ensembl_id)]
  } else {
    stop("Unsupported species in config: ", config$species)
  }

  return(id_map)
}


#' Parse Open Targets API Response
#' @noRd
.parseOtResponse <- function(g_dat, h_id, data_types, results) {

  # --- Diseases ---
  if ("diseases" %in% data_types &&
      !is.null(g_dat$associatedDiseases$rows) &&
      length(g_dat$associatedDiseases$rows) > 0) {
    d_dt <- data.table::as.data.table(g_dat$associatedDiseases$rows)

    # Flatten nested disease columns
    if ("disease.id" %in% names(d_dt)) {
      data.table::setnames(d_dt, old = c("disease.id", "disease.name"),
                           new = c("disease_id", "disease_name"), skip_absent = TRUE)
    }
    # Extract therapeutic areas if present
    if ("disease" %in% names(d_dt) && is.list(d_dt$disease)) {
      tryCatch({
        d_dt[, disease_id := sapply(disease, function(x) if (is.list(x)) x$id else NA_character_)]
        d_dt[, disease_name := sapply(disease, function(x) if (is.list(x)) x$name else NA_character_)]
        d_dt[, therapeutic_areas := sapply(disease, function(x) {
          if (is.list(x) && !is.null(x$therapeuticAreas)) {
            paste(sapply(x$therapeuticAreas, function(ta) ta$name), collapse = "; ")
          } else NA_character_
        })]
      }, error = function(e) NULL)
    }

    d_dt[, human_ensembl_id := h_id]
    results$diseases[[length(results$diseases) + 1]] <- d_dt
  }

  # --- Constraints ---
  if ("constraints" %in% data_types &&
      !is.null(g_dat$geneticConstraint) &&
      length(g_dat$geneticConstraint) > 0) {
    c_dt <- data.table::as.data.table(g_dat$geneticConstraint)
    c_dt[, human_ensembl_id := h_id]
    results$constraints[[length(results$constraints) + 1]] <- c_dt
  }

  # --- Tractability ---
  if ("tractability" %in% data_types &&
      !is.null(g_dat$tractability) &&
      length(g_dat$tractability) > 0) {
    t_dt <- data.table::as.data.table(g_dat$tractability)
    t_dt[, human_ensembl_id := h_id]
    results$tractability[[length(results$tractability) + 1]] <- t_dt
  }

  # --- Expression ---
  if ("expression" %in% data_types &&
      !is.null(g_dat$expressions) &&
      length(g_dat$expressions) > 0) {
    # Guard: fromJSON may auto-simplify to data.frame; convert to list-of-rows
    expr_input <- g_dat$expressions
    if (is.data.frame(expr_input)) {
      expr_input <- lapply(seq_len(nrow(expr_input)), function(i) {
        lapply(expr_input[i, ], function(x) if (is.list(x)) x[[1]] else x)
      })
    }
    expr_list <- lapply(expr_input, function(ex) {
      row <- data.table::data.table(
        tissue_id = ex$tissue$id %||% NA_character_,
        tissue_label = ex$tissue$label %||% NA_character_,
        anatomical_systems = paste(ex$tissue$anatomicalSystems %||% "", collapse = "; "),
        organs = paste(ex$tissue$organs %||% "", collapse = "; "),
        rna_level = ex$rna$level %||% NA_integer_,
        rna_value = ex$rna$value %||% NA_real_,
        rna_unit = ex$rna$unit %||% NA_character_,
        rna_zscore = ex$rna$zscore %||% NA_real_,
        protein_level = ex$protein$level %||% NA_integer_,
        protein_reliability = ex$protein$reliability %||% NA_character_
      )
      # Add cell type info if available
      if (!is.null(ex$protein$cellType) && length(ex$protein$cellType) > 0) {
        row$protein_cell_types <- paste(sapply(ex$protein$cellType, function(ct) ct$name), collapse = "; ")
      } else {
        row$protein_cell_types <- NA_character_
      }
      row
    })
    e_dt <- data.table::rbindlist(expr_list, fill = TRUE)
    e_dt[, human_ensembl_id := h_id]
    results$expression[[length(results$expression) + 1]] <- e_dt
  }

  # --- Mouse Phenotypes ---
  if ("mouse_phenotypes" %in% data_types &&
      !is.null(g_dat$mousePhenotypes) &&
      length(g_dat$mousePhenotypes) > 0) {
    # Guard: fromJSON may auto-simplify to data.frame; convert to list-of-rows
    mp_input <- g_dat$mousePhenotypes
    if (is.data.frame(mp_input)) {
      mp_input <- lapply(seq_len(nrow(mp_input)), function(i) {
        lapply(mp_input[i, ], function(x) if (is.list(x)) x[[1]] else x)
      })
    }
    mp_list <- lapply(mp_input, function(mp) {
      row <- data.table::data.table(
        phenotype_label = mp$modelPhenotypeLabel %||% NA_character_,
        phenotype_id = mp$modelPhenotypeId %||% NA_character_,
        target_from_source_id = mp$targetFromSourceId %||% NA_character_,
        target_in_model = mp$targetInModel %||% NA_character_,
        target_mgi_id = mp$targetInModelMgiId %||% NA_character_
      )
      # Add phenotype classes
      if (!is.null(mp$modelPhenotypeClasses) && length(mp$modelPhenotypeClasses) > 0) {
        row$phenotype_classes <- paste(sapply(mp$modelPhenotypeClasses, function(pc) pc$label), collapse = "; ")
      } else {
        row$phenotype_classes <- NA_character_
      }
      # Add biological models
      if (!is.null(mp$biologicalModels) && length(mp$biologicalModels) > 0) {
        row$allelic_composition <- paste(sapply(mp$biologicalModels, function(bm) bm$allelicComposition), collapse = "; ")
        row$genetic_background <- paste(sapply(mp$biologicalModels, function(bm) bm$geneticBackground), collapse = "; ")
      } else {
        row$allelic_composition <- NA_character_
        row$genetic_background <- NA_character_
      }
      row
    })
    mp_dt <- data.table::rbindlist(mp_list, fill = TRUE)
    mp_dt[, human_ensembl_id := h_id]
    results$mouse_phenotypes[[length(results$mouse_phenotypes) + 1]] <- mp_dt
  }

  # --- Homologues ---
  if ("homologues" %in% data_types &&
      !is.null(g_dat$homologues) &&
      length(g_dat$homologues) > 0) {
    h_dt <- data.table::as.data.table(g_dat$homologues)
    data.table::setnames(h_dt,
                         old = c("speciesId", "speciesName", "homologyType",
                                 "queryPercentageIdentity", "targetPercentageIdentity",
                                 "targetGeneId", "targetGeneSymbol", "isHighConfidence"),
                         new = c("species_id", "species_name", "homology_type",
                                 "query_pct_identity", "target_pct_identity",
                                 "target_gene_id", "target_gene_symbol", "is_high_confidence"),
                         skip_absent = TRUE)
    h_dt[, human_ensembl_id := h_id]
    results$homologues[[length(results$homologues) + 1]] <- h_dt
  }

  # --- Interactions ---
  if ("interactions" %in% data_types &&
      !is.null(g_dat$interactions$rows) &&
      length(g_dat$interactions$rows) > 0) {
    i_dt <- data.table::as.data.table(g_dat$interactions$rows)
    data.table::setnames(i_dt,
                         old = c("intA", "intABiologicalRole", "intB", "intBBiologicalRole",
                                 "score", "sourceDatabase"),
                         new = c("interactor_a", "role_a", "interactor_b", "role_b",
                                 "interaction_score", "source_database"),
                         skip_absent = TRUE)
    i_dt[, human_ensembl_id := h_id]
    results$interactions[[length(results$interactions) + 1]] <- i_dt
  }

  # --- Known Drugs ---
  if ("known_drugs" %in% data_types &&
      !is.null(g_dat$knownDrugs$rows) &&
      length(g_dat$knownDrugs$rows) > 0) {
    # Guard: fromJSON may auto-simplify to data.frame; convert to list-of-rows
    kd_input <- g_dat$knownDrugs$rows
    if (is.data.frame(kd_input)) {
      kd_input <- lapply(seq_len(nrow(kd_input)), function(i) {
        lapply(kd_input[i, ], function(x) if (is.list(x)) x[[1]] else x)
      })
    }
    kd_list <- lapply(kd_input, function(kd) {
      data.table::data.table(
        drug_id = kd$drugId %||% kd$drug$id %||% NA_character_,
        drug_name = kd$drug$name %||% NA_character_,
        drug_type = kd$drug$drugType %||% NA_character_,
        max_clinical_phase = kd$drug$maximumClinicalTrialPhase %||% NA_integer_,
        is_approved = kd$drug$isApproved %||% NA,
        mechanism_of_action = kd$mechanismOfAction %||% NA_character_,
        disease_id = kd$disease$id %||% NA_character_,
        disease_name = kd$disease$name %||% NA_character_,
        trial_phase = kd$phase %||% NA_integer_,
        trial_status = kd$status %||% NA_character_
      )
    })
    kd_dt <- data.table::rbindlist(kd_list, fill = TRUE)
    kd_dt[, human_ensembl_id := h_id]
    results$known_drugs[[length(results$known_drugs) + 1]] <- kd_dt
  }

  # --- DepMap ---
  if ("depmap" %in% data_types &&
      !is.null(g_dat$depMapEssentiality$screens) &&
      length(g_dat$depMapEssentiality$screens) > 0) {
    # Guard: fromJSON may auto-simplify to data.frame with wide-format columns;
    # extract only the expected columns to prevent column explosion
    dm_input <- g_dat$depMapEssentiality$screens
    if (is.data.frame(dm_input)) {
      expected_cols <- c("depmapId", "cellLineName", "diseaseFromSource",
                         "mutation", "geneEffect", "expression")
      avail_cols <- intersect(expected_cols, names(dm_input))
      dm_dt <- data.table::as.data.table(dm_input[, avail_cols, drop = FALSE])
    } else {
      dm_list <- lapply(dm_input, function(sc) {
        data.table::data.table(
          depmapId = sc$depmapId %||% NA_character_,
          cellLineName = sc$cellLineName %||% NA_character_,
          diseaseFromSource = sc$diseaseFromSource %||% NA_character_,
          mutation = sc$mutation %||% NA_character_,
          geneEffect = sc$geneEffect %||% NA_real_,
          expression = sc$expression %||% NA_real_
        )
      })
      dm_dt <- data.table::rbindlist(dm_list, fill = TRUE)
    }
    data.table::setnames(dm_dt,
                         old = c("depmapId", "cellLineName", "diseaseFromSource",
                                 "mutation", "geneEffect", "expression"),
                         new = c("depmap_id", "cell_line", "disease",
                                 "mutation", "gene_effect", "expression"),
                         skip_absent = TRUE)
    dm_dt[, human_ensembl_id := h_id]
    results$depmap[[length(results$depmap) + 1]] <- dm_dt
  }

  # --- Pharmacogenomics ---
  if ("pharmacogenomics" %in% data_types &&
      !is.null(g_dat$pharmacogenomics) &&
      length(g_dat$pharmacogenomics) > 0) {
    # Guard: fromJSON may auto-simplify to data.frame; convert to list-of-rows
    pg_input <- g_dat$pharmacogenomics
    if (is.data.frame(pg_input)) {
      pg_input <- lapply(seq_len(nrow(pg_input)), function(i) {
        lapply(pg_input[i, ], function(x) if (is.list(x)) x[[1]] else x)
      })
    }
    pg_list <- lapply(pg_input, function(pg) {
      # Extract variant consequence from SequenceOntologyTerm object
      variant_consequence <- if (is.list(pg$variantFunctionalConsequence)) {
        pg$variantFunctionalConsequence$label %||% pg$variantFunctionalConsequence$id %||% NA_character_
      } else {
        pg$variantFunctionalConsequence %||% NA_character_
      }
      # Extract drugs from nested list
      drug_ids <- NA_character_
      drug_names <- NA_character_
      if (!is.null(pg$drugs) && length(pg$drugs) > 0) {
        if (is.data.frame(pg$drugs)) {
          drug_ids <- paste(pg$drugs$drugId[!is.na(pg$drugs$drugId)], collapse = "; ")
          drug_names <- paste(pg$drugs$drugFromSource[!is.na(pg$drugs$drugFromSource)], collapse = "; ")
        } else if (is.list(pg$drugs)) {
          drug_ids <- paste(sapply(pg$drugs, function(d) d$drugId %||% NA_character_), collapse = "; ")
          drug_names <- paste(sapply(pg$drugs, function(d) d$drugFromSource %||% NA_character_), collapse = "; ")
        }
      }
      data.table::data.table(
        variant_rsid = pg$variantRsId %||% NA_character_,
        variant_consequence = variant_consequence,
        genotype_id = pg$genotypeId %||% NA_character_,
        drug_id = drug_ids,
        drug_name = drug_names,
        phenotype_text = pg$phenotypeText %||% NA_character_,
        genotype_annotation = pg$genotypeAnnotationText %||% NA_character_,
        phenotype_source_id = pg$phenotypeFromSourceId %||% NA_character_,
        evidence_level = pg$evidenceLevel %||% NA_character_,
        datasource_id = pg$datasourceId %||% NA_character_,
        study_id = pg$studyId %||% NA_character_
      )
    })
    pg_dt <- data.table::rbindlist(pg_list, fill = TRUE)
    pg_dt[, human_ensembl_id := h_id]
    results$pharmacogenomics[[length(results$pharmacogenomics) + 1]] <- pg_dt
  }

  return(results)
}


#' Null-coalescing operator
#' @noRd
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x


#' Query Open Targets for QTL Data
#'
#' Queries the Open Targets Genetics API for quantitative trait locus (QTL) data
#' including eQTL, pQTL, scQTL, and other molecular QTL types from credible sets.
#' Queries are performed by genomic region (for both gene and region mode projects).
#'
#' @param project_dir Character path to the project directory
#' @param study_types Character vector specifying which QTL study types to retrieve.
#'   Available options: \code{"eqtl"} (bulk expression), \code{"pqtl"} (bulk protein),
#'   \code{"sceqtl"} (single-cell expression), \code{"scpqtl"} (single-cell protein),
#'   \code{"sqtl"} (splicing), \code{"scsqtl"} (single-cell splicing),
#'   \code{"tuqtl"} (transcript usage), \code{"sctuqtl"} (single-cell transcript usage).
#'   Default includes bulk and single-cell eQTL/pQTL.
#' @param include_colocalisation Logical. Include GWAS colocalisation data. Default FALSE.
#' @param include_l2g Logical. Include Locus-to-Gene (L2G) prediction scores. Default TRUE.
#' @param min_l2g_score Numeric. Minimum L2G score threshold (0-1). Default 0.5.
#' @param limit Integer (optional). Limit query to first N regions for testing.
#' @param page_size Integer. Number of credible sets per API request. Default 500.
#' @param verbose Logical. Print progress messages. Default TRUE.
#'
#' @return Invisible list of data.tables with QTL results
#'
#' @details
#' This function queries the Open Targets Genetics credibleSets endpoint to retrieve
#' fine-mapped QTL associations. It provides:
#'
#' \itemize{
#'   \item \strong{eQTL data}: Expression QTL from GTEx, eQTLGen, and other sources
#'   \item \strong{pQTL data}: Protein QTL from deCODE, INTERVAL, etc.
#'   \item \strong{scQTL data}: Single-cell QTL addressing reviewer concern R2.2
#'   \item \strong{L2G scores}: Machine learning predictions linking variants to causal genes
#'   \item \strong{Colocalisation}: Evidence that QTL and GWAS signals share the same causal variant
#' }
#'
#' **Output tables:**
#' \itemize{
#'   \item \code{ot_qtl_studies} (otqs): Study metadata including tissue, cell type, sample size
#'   \item \code{ot_qtl_credible_sets} (otqc): Fine-mapped loci with lead variants and statistics
#'   \item \code{ot_qtl_coloc} (otql): Colocalisation with GWAS traits (if requested)
#' }
#'
#' @examples
#' \dontrun{
#' # Query for all bulk and single-cell eQTL data
#' queryOpenTargetsQTL(project_dir = "my_analysis",
#'                     study_types = c("eqtl", "sceqtl"))
#'
#' # Include colocalisation with GWAS
#' queryOpenTargetsQTL(project_dir = "my_analysis",
#'                     include_colocalisation = TRUE)
#'
#' # Higher L2G threshold for stringent gene assignment
#' queryOpenTargetsQTL(project_dir = "my_analysis",
#'                     min_l2g_score = 0.7)
#'
#' # Query all available QTL types
#' queryOpenTargetsQTL(project_dir = "my_analysis",
#'                     study_types = c("eqtl", "pqtl", "sceqtl", "scpqtl",
#'                                     "sqtl", "scsqtl"))
#' }
#'
#' @importFrom stats na.omit
#' @export
queryOpenTargetsQTL <- function(project_dir = ".",
                                study_types = c("eqtl", "pqtl", "sceqtl", "scpqtl"),
                                include_colocalisation = FALSE,
                                include_l2g = TRUE,
                                min_l2g_score = 0.5,
                                limit = NULL,
                                page_size = 500,
                                verbose = TRUE) {

  # 1. Setup and Checks
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for this function.")
  }

  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) stop("Project not found.")

  # Validate study_types
  valid_types <- c("eqtl", "pqtl", "sceqtl", "scpqtl", "sqtl", "scsqtl", "tuqtl", "sctuqtl")
  invalid <- setdiff(tolower(study_types), valid_types)
  if (length(invalid) > 0) {
    stop("Invalid study_types: ", paste(invalid, collapse = ", "),
         "\nValid options: ", paste(valid_types, collapse = ", "))
  }
  study_types <- tolower(study_types)

  # 2. Read Config
  config_path <- file.path(packrat_dir, "config.json")
  if (!file.exists(config_path)) stop("config.json not found in project.")
  config <- jsonlite::read_json(config_path)

  if (verbose) message(sprintf("Project Config: %s - %s (mode: %s)",
                               config$species, config$genome, config$mode))

  # 3. Get regions from project (works for both gene and region mode)
  regions_dt <- .getProjectRegions(packrat_dir, config)

  if (nrow(regions_dt) == 0) {
    warning("No valid regions found for querying.")
    return(NULL)
  }

  # Build region strings in API format: "chr:start-end"
  # Note: Open Targets expects "chr" prefix in region queries
  regions_dt[, chr_clean := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
  region_strings <- sprintf("%s:%d-%d", regions_dt$chr_clean, regions_dt$start, regions_dt$end)

  # Store mapping for later linking results back to project
  region_mapping <- data.table::data.table(
    region_string = region_strings,
    region_id = regions_dt$region_id,
    chr = regions_dt$chr,
    start = regions_dt$start,
    end = regions_dt$end
  )

  if (!is.null(limit)) {
    if (verbose) message(sprintf("TEST MODE: Limiting to first %d regions", limit))
    region_strings <- head(region_strings, limit)
    region_mapping <- head(region_mapping, limit)
  }

  # 4. Build GraphQL query for credible sets (region-based)
  build_qtl_query <- function(region_strs, types, incl_l2g, incl_coloc, pg_size) {
    # Study types as unquoted lowercase enums
    type_filter <- paste(tolower(types), collapse = ", ")

    # Regions as quoted strings
    regions_str <- paste(sprintf('"%s"', region_strs), collapse = ", ")

    # L2G sub-query
    l2g_fields <- if (incl_l2g) {
      "l2GPredictions(page: { size: 10, index: 0 }) { rows { studyLocusId score target { id approvedSymbol } } }"
    } else ""

    # Colocalisation sub-query
    coloc_fields <- if (incl_coloc) {
      "colocalisation(page: { size: 50, index: 0 }) {
        rows {
          rightStudyLocusId leftStudyLocusId rightStudyType
          colocalisationMethod h3 h4
          otherStudyLocus { studyId study { traitFromSource } }
        }
      }"
    } else ""

    sprintf('
      query credibleSetsForRegions {
        credibleSets(
          regions: [%s]
          studyTypes: [%s]
          page: { size: %d, index: 0 }
        ) {
          count
          rows {
            studyLocusId
            studyId
            studyType
            chromosome
            position
            region
            locusStart
            locusEnd
            qtlGeneId
            isTransQtl
            sampleSize
            pValueMantissa
            pValueExponent
            beta
            standardError
            effectAlleleFrequencyFromSource
            finemappingMethod
            confidence
            credibleSetIndex
            variant {
              id
              chromosome
              position
              referenceAllele
              alternateAllele
              rsIds
            }
            study {
              id
              studyType
              traitFromSource
              projectId
              hasSumstats
              nSamples
              biosample { biosampleId biosampleName }
              target { id approvedSymbol }
            }
            %s
            %s
          }
        }
      }
    ', regions_str, type_filter, pg_size, l2g_fields, coloc_fields)
  }

  # 5. Query by batched regions
  endpoint <- "https://api.platform.opentargets.org/api/v4/graphql"

  results_studies <- list()
  results_credsets <- list()
  results_coloc <- list()

  if (verbose) {
    message(sprintf("Querying Open Targets QTL for %d regions...", length(region_strings)))
    message(sprintf("Study types: %s", paste(toupper(study_types), collapse = ", ")))
  }

  # Batch regions (API may have limits)
  batch_size <- 20
  batches <- split(region_strings, ceiling(seq_along(region_strings) / batch_size))

  for (i in seq_along(batches)) {
    batch_regions <- batches[[i]]
    if (verbose) message(sprintf("  Processing batch %d/%d (%d regions)...",
                                  i, length(batches), length(batch_regions)))

    tryCatch({
      query <- build_qtl_query(batch_regions, study_types, include_l2g,
                               include_colocalisation, page_size)

      response <- httr::POST(
        endpoint,
        body = list(query = query),
        encode = "json",
        httr::timeout(120)
      )

      # Check for errors
      if (httr::status_code(response) != 200) {
        content <- httr::content(response, "text", encoding = "UTF-8")
        if (verbose) warning(sprintf("Batch %d failed with status %d: %s",
                                     i, httr::status_code(response), content))
        next
      }

      content <- httr::content(response, "text", encoding = "UTF-8")
      data_json <- jsonlite::fromJSON(content, flatten = FALSE)

      # Check for GraphQL errors
      if (!is.null(data_json$errors)) {
        if (verbose) warning(sprintf("Batch %d GraphQL error: %s",
                                     i, data_json$errors[[1]]$message))
        next
      }

      cs_data <- data_json$data$credibleSets
      if (is.null(cs_data) || is.null(cs_data$rows) || length(cs_data$rows) == 0) {
        next
      }

      rows <- cs_data$rows

      # Parse each credible set row
      for (j in seq_len(if (is.data.frame(rows)) nrow(rows) else length(rows))) {
        row <- if (is.data.frame(rows)) as.list(rows[j, ]) else rows[[j]]

        # Extract study info
        study <- row$study
        study_dt <- data.table::data.table(
          study_id = row$studyId %||% NA_character_,
          study_type = row$studyType %||% NA_character_,
          trait = study$traitFromSource %||% NA_character_,
          project_id = study$projectId %||% NA_character_,
          biosample_id = study$biosample$biosampleId %||% NA_character_,
          biosample_name = study$biosample$biosampleName %||% NA_character_,
          n_samples = row$sampleSize %||% NA_integer_,
          has_sumstats = study$hasSumstats %||% NA,
          target_gene_id = study$target$id %||% NA_character_,
          target_gene_symbol = study$target$approvedSymbol %||% NA_character_
        )
        results_studies[[length(results_studies) + 1]] <- study_dt

        # Extract variant info
        variant <- row$variant

        # Extract credible set info
        credset_dt <- data.table::data.table(
          study_locus_id = row$studyLocusId %||% NA_character_,
          study_id = row$studyId %||% NA_character_,
          study_type = row$studyType %||% NA_character_,
          qtl_gene_id = row$qtlGeneId %||% NA_character_,
          is_trans_qtl = row$isTransQtl %||% NA,
          chromosome = row$chromosome %||% NA_character_,
          position = row$position %||% NA_integer_,
          region = row$region %||% NA_character_,
          locus_start = row$locusStart %||% NA_integer_,
          locus_end = row$locusEnd %||% NA_integer_,
          lead_variant_id = variant$id %||% NA_character_,
          ref_allele = variant$referenceAllele %||% NA_character_,
          alt_allele = variant$alternateAllele %||% NA_character_,
          rsid = if (!is.null(variant$rsIds) && length(variant$rsIds) > 0) variant$rsIds[1] else NA_character_,
          pvalue_mantissa = row$pValueMantissa %||% NA_real_,
          pvalue_exponent = row$pValueExponent %||% NA_integer_,
          beta = row$beta %||% NA_real_,
          se = row$standardError %||% NA_real_,
          eaf = row$effectAlleleFrequencyFromSource %||% NA_real_,
          finemapping_method = row$finemappingMethod %||% NA_character_,
          confidence = row$confidence %||% NA_character_,
          credible_set_index = row$credibleSetIndex %||% NA_integer_,
          biosample_id = study$biosample$biosampleId %||% NA_character_,
          biosample_name = study$biosample$biosampleName %||% NA_character_
        )

        # Add L2G scores if available
        if (include_l2g && !is.null(row$l2GPredictions$rows) && length(row$l2GPredictions$rows) > 0) {
          # Get top L2G prediction
          l2g_rows <- row$l2GPredictions$rows
          if (is.data.frame(l2g_rows) && nrow(l2g_rows) > 0) {
            top_l2g <- l2g_rows[1, ]
            credset_dt$l2g_gene_id <- top_l2g$target$id %||% NA_character_
            credset_dt$l2g_gene_symbol <- top_l2g$target$approvedSymbol %||% NA_character_
            credset_dt$l2g_score <- top_l2g$score %||% NA_real_
          } else if (is.list(l2g_rows) && length(l2g_rows) > 0) {
            top_l2g <- l2g_rows[[1]]
            credset_dt$l2g_gene_id <- top_l2g$target$id %||% NA_character_
            credset_dt$l2g_gene_symbol <- top_l2g$target$approvedSymbol %||% NA_character_
            credset_dt$l2g_score <- top_l2g$score %||% NA_real_
          }
        }

        results_credsets[[length(results_credsets) + 1]] <- credset_dt

        # Extract colocalisation if requested
        if (include_colocalisation && !is.null(row$colocalisation$rows)) {
          coloc_rows <- row$colocalisation$rows
          coloc_list <- if (is.data.frame(coloc_rows)) {
            lapply(seq_len(nrow(coloc_rows)), function(k) as.list(coloc_rows[k, ]))
          } else {
            coloc_rows
          }

          for (coloc in coloc_list) {
            other <- coloc$otherStudyLocus
            coloc_dt <- data.table::data.table(
              study_locus_id = row$studyLocusId %||% NA_character_,
              qtl_study_id = row$studyId %||% NA_character_,
              right_study_locus_id = coloc$rightStudyLocusId %||% NA_character_,
              right_study_type = coloc$rightStudyType %||% NA_character_,
              gwas_study_id = other$studyId %||% NA_character_,
              gwas_trait = other$study$traitFromSource %||% NA_character_,
              coloc_method = coloc$colocalisationMethod %||% NA_character_,
              h3 = coloc$h3 %||% NA_real_,
              h4 = coloc$h4 %||% NA_real_
            )
            results_coloc[[length(results_coloc) + 1]] <- coloc_dt
          }
        }
      }

    }, error = function(e) {
      if (verbose) warning(sprintf("Batch %d failed: %s", i, e$message))
    })

    Sys.sleep(0.3)  # Rate limiting
  }

  # 6. Link results back to project and save

  # Build comprehensive id_map from both project genes AND QTL response gene IDs
  gene_syms <- .getProjectGenes(packrat_dir)
  id_map <- .resolveHumanIds(gene_syms, config)

  # Collect all ensembl IDs from QTL results that may not be in the project gene list
  extra_ids <- character(0)
  if (length(results_studies) > 0) {
    studies_tmp <- data.table::rbindlist(results_studies, fill = TRUE, use.names = TRUE)
    if ("target_gene_id" %in% names(studies_tmp)) {
      extra_ids <- c(extra_ids, stats::na.omit(unique(studies_tmp$target_gene_id)))
    }
  }
  if (length(results_credsets) > 0) {
    credsets_tmp <- data.table::rbindlist(results_credsets, fill = TRUE, use.names = TRUE)
    if ("qtl_gene_id" %in% names(credsets_tmp)) {
      extra_ids <- c(extra_ids, stats::na.omit(unique(credsets_tmp$qtl_gene_id)))
    }
    if ("l2g_gene_id" %in% names(credsets_tmp)) {
      extra_ids <- c(extra_ids, stats::na.omit(unique(credsets_tmp$l2g_gene_id)))
    }
  }
  extra_ids <- setdiff(unique(extra_ids), id_map$human_ensembl_id)

  # Look up extra IDs in the human coordinate reference file
  if (length(extra_ids) > 0) {
    coords_file <- system.file("extdata", "human_coords_hg38.csv", package = "locusPackRat")
    if (coords_file == "") coords_file <- "inst/extdata/human_coords_hg38.csv"
    if (file.exists(coords_file)) {
      ref_coords <- data.table::fread(coords_file)
      extra_map <- ref_coords[ref_coords$ensembl_id %in% extra_ids,
                              .(gene_symbol, human_ensembl_id = ensembl_id)]
      if (nrow(extra_map) > 0) {
        id_map <- data.table::rbindlist(list(id_map, extra_map), fill = TRUE, use.names = TRUE)
        id_map <- unique(id_map, by = "human_ensembl_id")
      }
    }
  }

  saved_tables <- list()

  # Studies table
  if (length(results_studies) > 0) {
    studies_dt <- data.table::rbindlist(results_studies, fill = TRUE, use.names = TRUE)
    studies_dt <- unique(studies_dt)

    # Link by target gene if available
    if ("target_gene_id" %in% names(studies_dt)) {
      studies_dt <- merge(studies_dt,
                          id_map[, .(human_ensembl_id, gene_symbol)],
                          by.x = "target_gene_id", by.y = "human_ensembl_id",
                          all.x = TRUE)
    }

    if ("gene_symbol" %in% names(studies_dt)) {
      col_order <- c("gene_symbol", setdiff(names(studies_dt), "gene_symbol"))
      data.table::setcolorder(studies_dt, col_order)
    }

    addRatTable(
      data = studies_dt,
      table_name = "ot_qtl_studies",
      link_type = "gene",
      link_by = "gene_symbol",
      abbreviation = "otqs",
      project_dir = project_dir
    )
    saved_tables$studies <- studies_dt
    if (verbose) message(sprintf("  Saved ot_qtl_studies with %d rows", nrow(studies_dt)))
  }

  # Credible sets table
  if (length(results_credsets) > 0) {
    credsets_dt <- data.table::rbindlist(results_credsets, fill = TRUE, use.names = TRUE)

    # Link by qtlGeneId (the gene being studied in molQTL)
    if ("qtl_gene_id" %in% names(credsets_dt)) {
      credsets_dt <- merge(credsets_dt,
                           id_map[, .(human_ensembl_id, gene_symbol)],
                           by.x = "qtl_gene_id", by.y = "human_ensembl_id",
                           all.x = TRUE)
    }

    # Filter by L2G score if specified
    if (include_l2g && "l2g_score" %in% names(credsets_dt)) {
      before_filter <- nrow(credsets_dt)
      credsets_dt <- credsets_dt[is.na(l2g_score) | l2g_score >= min_l2g_score]
      if (verbose && before_filter > nrow(credsets_dt)) {
        message(sprintf("  Filtered %d rows by L2G score threshold (%.2f)",
                        before_filter - nrow(credsets_dt), min_l2g_score))
      }
    }

    if ("gene_symbol" %in% names(credsets_dt)) {
      col_order <- c("gene_symbol", setdiff(names(credsets_dt), "gene_symbol"))
      data.table::setcolorder(credsets_dt, col_order)
    }

    addRatTable(
      data = credsets_dt,
      table_name = "ot_qtl_credible_sets",
      link_type = "gene",
      link_by = "gene_symbol",
      abbreviation = "otqc",
      project_dir = project_dir
    )
    saved_tables$credible_sets <- credsets_dt
    if (verbose) message(sprintf("  Saved ot_qtl_credible_sets with %d rows", nrow(credsets_dt)))
  }

  # Colocalisation table
  if (include_colocalisation && length(results_coloc) > 0) {
    coloc_dt <- data.table::rbindlist(results_coloc, fill = TRUE, use.names = TRUE)

    # Link coloc results to credible sets to get gene info
    if (nrow(coloc_dt) > 0 && exists("credsets_dt") && "study_locus_id" %in% names(credsets_dt)) {
      gene_link <- unique(credsets_dt[, .(study_locus_id, gene_symbol)])
      coloc_dt <- merge(coloc_dt, gene_link, by = "study_locus_id", all.x = TRUE)
    }

    if ("gene_symbol" %in% names(coloc_dt)) {
      col_order <- c("gene_symbol", setdiff(names(coloc_dt), "gene_symbol"))
      data.table::setcolorder(coloc_dt, col_order)
    }

    addRatTable(
      data = coloc_dt,
      table_name = "ot_qtl_coloc",
      link_type = "gene",
      link_by = "gene_symbol",
      abbreviation = "otql",
      project_dir = project_dir
    )
    saved_tables$colocalisation <- coloc_dt
    if (verbose) message(sprintf("  Saved ot_qtl_coloc with %d rows", nrow(coloc_dt)))
  }

  if (length(saved_tables) == 0) {
    warning("No QTL data found for the specified regions and study types.")
  }

  invisible(saved_tables)
}
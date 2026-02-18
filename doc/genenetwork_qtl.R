## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## ----libraries----------------------------------------------------------------
# library(locusPackRat)
# library(data.table)
# library(httr)
# library(jsonlite)

## ----api_base-----------------------------------------------------------------
# # Base URL for GeneNetwork2 API
# gn2_base <- "https://genenetwork.org/api/v_pre1"

## ----list_bxd_datasets--------------------------------------------------------
# # Fetch all datasets for the BXD group
# bxd_datasets <- fromJSON(content(
#   GET(paste0(gn2_base, "/datasets/BXD")),
#   as = "text", encoding = "UTF-8"
# ))
# 
# # Filter for expression datasets (mRNA / RNA-seq)
# expression_datasets <- bxd_datasets[grepl("mRNA|RNA-seq|expression",
#                                           bxd_datasets$tissue, ignore.case = TRUE), ]
# 
# # Preview available expression datasets
# display_cols <- intersect(names(expression_datasets), c("id", "name", "tissue", "description"))
# head(expression_datasets[, display_cols, drop = FALSE])

## ----query_gene_traits--------------------------------------------------------
# # Example: query Cisd2 expression traits in BXD heart dataset
# # Dataset ID will vary — use the listing above to find the correct one
# dataset_name <- "HC_M2_0606_P"  # Example: heart dataset
# 
# # Search for traits matching a gene symbol
# trait_search_url <- paste0(gn2_base, "/traits/", dataset_name)
# heart_traits <- fromJSON(content(
#   GET(trait_search_url),
#   as = "text", encoding = "UTF-8"
# ))
# 
# # Filter for a specific gene of interest
# gene_of_interest <- "Cisd2"
# cisd2_traits <- heart_traits[grepl(gene_of_interest,
#                                     heart_traits$symbol, ignore.case = TRUE), ]
# 
# if (nrow(cisd2_traits) > 0) {
#   cat("Found", nrow(cisd2_traits), "probes/traits for", gene_of_interest, "\n")
#   print(cisd2_traits[, c("id", "symbol", "description")])
# }

## ----get_trait_values---------------------------------------------------------
# # Get expression values for a specific trait across BxD strains
# trait_id <- cisd2_traits$id[1]  # Use the first matching probe
# 
# trait_url <- paste0(gn2_base, "/trait/", dataset_name, "/", trait_id)
# trait_data <- fromJSON(content(
#   GET(trait_url),
#   as = "text", encoding = "UTF-8"
# ))
# 
# # Convert to data.table for locusPackRat integration
# if (!is.null(trait_data$sample_data)) {
#   strain_expression <- data.table(
#     strain = names(trait_data$sample_data),
#     expression = unlist(trait_data$sample_data)
#   )
#   cat("Expression data for", nrow(strain_expression), "BxD strains\n")
#   head(strain_expression)
# }

## ----run_gemma----------------------------------------------------------------
# # Request GEMMA mapping for a trait
# # Note: this endpoint may take some time to return results
# gemma_url <- paste0(gn2_base, "/mapping?trait_id=", trait_id,
#                     "&db=", dataset_name,
#                     "&method=gemma")
# 
# gemma_results <- fromJSON(content(
#   GET(gemma_url, timeout(120)),
#   as = "text", encoding = "UTF-8"
# ))
# 
# # Convert to data.table
# if (length(gemma_results) > 0) {
#   eqtl_scan <- data.table(
#     chr = gemma_results$chr,
#     pos = gemma_results$Mb * 1e6,  # Convert Mb to bp
#     lod = gemma_results$lod_score,
#     gene_symbol = gene_of_interest,
#     source = "GeneNetwork2_BxD"
#   )
# 
#   # Filter for significant peaks (e.g., LOD > 3)
#   significant_eqtl <- eqtl_scan[lod > 3]
#   cat("Found", nrow(significant_eqtl), "significant eQTL positions\n")
# }

## ----list_do_datasets---------------------------------------------------------
# # Fetch datasets for the DO group
# do_datasets <- fromJSON(content(
#   GET(paste0(gn2_base, "/datasets/DO")),
#   as = "text", encoding = "UTF-8"
# ))
# 
# # Preview
# display_cols <- intersect(names(do_datasets), c("id", "name", "tissue", "description"))
# head(do_datasets[, display_cols, drop = FALSE])

## ----query_do_phenotypes------------------------------------------------------
# # Example: query a DO phenotype dataset
# do_dataset_name <- do_datasets$name[1]  # Use first available
# 
# do_traits <- fromJSON(content(
#   GET(paste0(gn2_base, "/traits/", do_dataset_name)),
#   as = "text", encoding = "UTF-8"
# ))
# 
# # Preview available traits
# if (nrow(do_traits) > 0) {
#   head(do_traits[, c("id", "symbol", "description")])
# }

## ----integrate_gene_level-----------------------------------------------------
# # Summarize eQTL results per gene (e.g., maximum cis-eQTL LOD)
# # In practice, you would loop over multiple candidate genes
# eqtl_summary <- data.table(
#   gene_symbol = c("Cisd2", "Pdlim5", "Manba"),
#   bxd_max_lod = c(5.2, 3.1, 4.8),
#   bxd_peak_chr = c("3", "3", "3"),
#   bxd_peak_pos = c(135000000, 142000000, 138000000),
#   bxd_cis_eqtl = c(TRUE, FALSE, TRUE),
#   gn2_dataset = rep("HC_M2_0606_P", 3)
# )
# 
# # Add to locusPackRat project
# addRatTable(
#   data = eqtl_summary,
#   table_name = "bxd_eqtl",
#   abbreviation = "bxd",
#   link_type = "gene",
#   link_by = "gene_symbol",
#   project_dir = "my_qtl_project"
# )

## ----integrate_point_level----------------------------------------------------
# # Use the full scan data from GEMMA
# # Ensure columns match the expected format: chr, pos (or start/end)
# scan_for_project <- eqtl_scan[, .(
#   chr = chr,
#   start = as.integer(pos),
#   end = as.integer(pos),
#   lod = lod,
#   gene_symbol = gene_symbol,
#   source = source
# )]
# 
# addRatTable(
#   data = scan_for_project,
#   table_name = "bxd_eqtl_scan",
#   abbreviation = "bxds",
#   link_type = "region",
#   link_by = "chr,start,end",
#   project_dir = "my_qtl_project"
# )

## ----gnapi_example------------------------------------------------------------
# # Install GNapi
# # remotes::install_github("kbroman/GNapi")
# 
# library(GNapi)
# 
# # List available species
# species <- list_species()
# 
# # List BxD datasets
# bxd_datasets <- list_datasets("BXD")
# 
# # Get trait data
# trait_data <- get_trait("HC_M2_0606_P", "your_trait_id")

## ----session_info-------------------------------------------------------------
# sessionInfo()


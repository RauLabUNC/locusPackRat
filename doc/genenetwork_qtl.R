## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE
)

## ----libraries----------------------------------------------------------------
library(locusPackRat)
library(data.table)
library(httr)
library(jsonlite)

## ----api_setup----------------------------------------------------------------
# Base URL for GeneNetwork2 API
gn2_base <- "https://genenetwork.org/api/v_pre1"

# Helper: fetch from API with cached fallback
gn_fetch_or_cache <- function(url, cache_file, timeout_sec = 30) {
  result <- tryCatch({
    resp <- httr::GET(url, httr::timeout(timeout_sec))
    httr::stop_for_status(resp)
    jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"))
  }, error = function(e) NULL)

  if (!is.null(result) && length(result) > 0) {
    message("GeneNetwork API: live data retrieved")
    return(result)
  }

  cache_path <- system.file("extdata", "gn_cache", cache_file,
                            package = "locusPackRat")
  if (nzchar(cache_path) && file.exists(cache_path)) {
    message("GeneNetwork API: using cached data (API unavailable)")
    return(jsonlite::fromJSON(cache_path))
  }
  stop("GeneNetwork API unavailable and no cached data found for: ", cache_file)
}

## ----list_bxd_datasets--------------------------------------------------------
# Fetch all datasets for the BXD group
bxd_datasets <- gn_fetch_or_cache(
  paste0(gn2_base, "/datasets/BXD"),
  "bxd_datasets.json"
)

# The FullName column contains tissue and platform info
# Filter for expression datasets related to heart/cardiac tissue
heart_datasets <- bxd_datasets[grepl("Heart|Cardiac|heart|cardiac",
                                     bxd_datasets$FullName), ]

cat("Total BXD datasets:", nrow(bxd_datasets), "\n")
cat("Heart-related datasets:", nrow(heart_datasets), "\n\n")

# Preview heart datasets
display_cols <- intersect(names(heart_datasets),
                          c("Id", "FullName", "Short_Abbreviation"))
head(heart_datasets[, display_cols, drop = FALSE])

## ----query_gene_traits--------------------------------------------------------
# Query Cisd2 expression traits in a BXD heart dataset
dataset_name <- "HC_M2_0606_P"
gene_of_interest <- "Cisd2"

# Fetch traits (pre-filtered to Cisd2 in cached version)
cisd2_traits <- gn_fetch_or_cache(
  paste0(gn2_base, "/traits/", dataset_name),
  "cisd2_traits.json"
)

# The live endpoint returns ALL traits; filter to our gene
# The cached version is already pre-filtered
if (nrow(cisd2_traits) > 20) {
  # Live data -- filter by Symbol column
  cisd2_traits <- cisd2_traits[grepl(gene_of_interest, cisd2_traits$Symbol,
                                     ignore.case = TRUE), ]
}

cat("Found", nrow(cisd2_traits), "probes/traits for", gene_of_interest, "\n\n")
print(cisd2_traits[, c("Name", "Symbol", "Description", "Chr", "Mb")])

## ----get_trait_values---------------------------------------------------------
# Get detail for the first (strongest) Cisd2 probe
probe_name <- cisd2_traits$Name[1]

trait_detail <- gn_fetch_or_cache(
  paste0(gn2_base, "/trait/", dataset_name, "/", probe_name),
  "cisd2_trait_detail.json"
)

cat("Probe:", trait_detail$name, "\n")
cat("Gene:", trait_detail$symbol, "\n")
cat("Description:", trait_detail$description, "\n")
cat("Chromosome:", trait_detail$chr, "at", trait_detail$mb, "Mb\n")
cat("Best LRS:", round(trait_detail$lrs, 2),
    "(LOD ~", round(trait_detail$lrs / 4.61, 2), ")\n")
cat("Peak locus:", trait_detail$locus, "\n")
cat("Additive effect:", round(trait_detail$additive, 4), "\n")

## ----run_gemma----------------------------------------------------------------
# Request GEMMA mapping for the Cisd2 probe
# Note: this endpoint may take 30-60 seconds to return results
gemma_url <- paste0(gn2_base, "/mapping?trait_id=", probe_name,
                    "&db=", dataset_name,
                    "&method=gemma")

gemma_raw <- gn_fetch_or_cache(gemma_url, "cisd2_gemma.json", timeout_sec = 120)

# The API returns a nested list; extract the data frame
if (is.data.frame(gemma_raw)) {
  gemma_results <- gemma_raw
} else if (is.list(gemma_raw) && length(gemma_raw) > 0) {
  gemma_results <- gemma_raw[[1]]
} else {
  gemma_results <- data.frame()
}

if (nrow(gemma_results) > 0) {
  eqtl_scan <- data.table(
    chr = gemma_results$chr,
    pos = gemma_results$Mb * 1e6,  # Convert Mb to bp
    pos_mb = gemma_results$Mb,
    lod = gemma_results$lod_score,
    p_value = gemma_results$p_value,
    marker = gemma_results$name,
    gene_symbol = gene_of_interest,
    source = "GeneNetwork2_BxD"
  )

  # Summary of the scan
  cat("Total markers scanned:", nrow(eqtl_scan), "\n")
  cat("Significant positions (LOD > 3):", sum(eqtl_scan$lod > 3), "\n\n")

  # Show the top peaks
  top_peaks <- eqtl_scan[order(-lod)][1:5]
  print(top_peaks[, .(chr, pos_mb, lod, p_value, marker)])
}

## ----list_bxd_phenotypes------------------------------------------------------
# Fetch BXD published phenotype traits
bxd_pheno <- tryCatch({
  gn_fetch_or_cache(
    paste0(gn2_base, "/traits/BXDPublish"),
    cache_file = NULL  # no cache for this endpoint
  )
}, error = function(e) {
  message("BXD published phenotypes unavailable: ", e$message)
  NULL
})

if (!is.null(bxd_pheno) && nrow(bxd_pheno) > 0) {
  cat("Total BXD published phenotype traits:", nrow(bxd_pheno), "\n\n")

  # Show a few example traits
  display_cols <- intersect(names(bxd_pheno), c("Id", "Name", "Description", "Symbol"))
  print(head(bxd_pheno[, display_cols, drop = FALSE], 10))
} else {
  cat("BXD published phenotype data not available via this API endpoint.\n")
  cat("Visit https://genenetwork.org to browse phenotypes interactively.\n")
}

## ----create_project-----------------------------------------------------------
# Create a temporary project with candidate genes
project_dir <- file.path(tempdir(), "gn_demo_project")

gene_input <- data.frame(
  gene_symbol = c("Cisd2", "Pdlim5", "Manba", "Fhod3", "Myh7")
)

initPackRat(
  data = gene_input,
  mode = "gene",
  species = "mouse",
  genome = "mm39",
  project_dir = project_dir,
  force = TRUE
)

## ----integrate_gene_level-----------------------------------------------------
# Build a gene-level summary from the GEMMA scan
# In a real analysis, you would loop over multiple candidate genes
if (exists("eqtl_scan") && nrow(eqtl_scan) > 0) {
  peak_row <- eqtl_scan[which.max(lod)]
  gene_chr <- unique(eqtl_scan[chr == peak_row$chr, chr])

  eqtl_summary <- data.table(
    gene_symbol = gene_of_interest,
    bxd_max_lod = round(peak_row$lod, 2),
    bxd_peak_chr = as.character(peak_row$chr),
    bxd_peak_mb = round(peak_row$pos_mb, 2),
    bxd_cis_eqtl = (as.character(peak_row$chr) == as.character(cisd2_traits$Chr[1])),
    gn2_dataset = dataset_name
  )

  print(eqtl_summary)

  # Add to locusPackRat project
  addRatTable(
    data = eqtl_summary,
    table_name = "bxd_eqtl",
    abbreviation = "bxd",
    link_type = "gene",
    link_by = "gene_symbol",
    project_dir = project_dir
  )
}

## ----integrate_point_level----------------------------------------------------
# Use the full scan data from GEMMA -- keep significant hits only
if (exists("eqtl_scan") && nrow(eqtl_scan) > 0) {
  significant_hits <- eqtl_scan[lod > 3]
  cat("Integrating", nrow(significant_hits), "significant markers\n")

  scan_for_project <- significant_hits[, .(
    chr = as.character(chr),
    start = as.integer(pos),
    end = as.integer(pos),
    lod = round(lod, 3),
    gene_symbol = gene_symbol,
    source = source
  )]

  addRatTable(
    data = scan_for_project,
    table_name = "bxd_eqtl_scan",
    abbreviation = "bxds",
    link_type = "region",
    link_by = "chr,start,end",
    project_dir = project_dir
  )
}

## ----show_project_state-------------------------------------------------------
listPackRatTables(project_dir)

## ----gnapi_example, eval = FALSE----------------------------------------------
# # Install GNapi (https://github.com/kbroman/GNapi)
# # remotes::install_github("kbroman/GNapi")
# library(GNapi)
# 
# # List available species and groups
# species <- list_species()
# groups  <- list_groups("mouse")
# 
# # List BxD datasets
# bxd_datasets <- list_datasets("BXD")
# 
# # Run GEMMA QTL scan for a specific trait
# gemma_scan <- run_gemma("HC_M2_0606_P", "1428441_at")
# 
# # Find correlated traits
# corr <- run_correlation("HC_M2_0606_P", "BXD", "1428441_at")

## ----session_info-------------------------------------------------------------
sessionInfo()

## ----cleanup, include = FALSE-------------------------------------------------
reg.finalizer(environment(), function(e) {
  unlink(project_dir, recursive = TRUE, force = TRUE)
}, onexit = TRUE)


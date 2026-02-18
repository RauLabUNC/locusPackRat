## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)

# Check for API availability
# Skip on CRAN (R CMD check sets _R_CHECK_PACKAGE_NAME_) or when API unreachable
on_cran <- nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))

# Check if we can reach the API endpoint (any response means it's reachable)
api_available <- requireNamespace("httr", quietly = TRUE) && {
  resp <- try(httr::GET("https://api.platform.opentargets.org/api/v4/graphql",
                        httr::timeout(10)), silent = TRUE)
  !inherits(resp, "try-error")
}

skip_api <- on_cran || !api_available

# Check for Bioconductor annotation packages needed for locus zoom plots
has_bioc_annot <- requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
                  requireNamespace("org.Hs.eg.db", quietly = TRUE)
can_plot <- !skip_api && has_bioc_annot

## ----libraries----------------------------------------------------------------
library(locusPackRat)
library(data.table)

## ----define_regions-----------------------------------------------------------
qtl_regions <- data.table(
  chr = c("19", "11", "16"),
  start = c(56769655, 33875225, 18285731),
  end = c(58769655, 35875225, 20285731),
  region_id = c("ZNF586_locus", "APIP_locus", "CLEC19A_locus"),
  peak_lod = c(8.2, 6.5, 7.1),
  trait = c("gene_expression", "gene_expression", "protein_levels")
)

# Display region summary
qtl_regions[, .(
  region_id,
  chr,
  size_kb = round((end - start) / 1000),
  trait
)]

## ----init_project-------------------------------------------------------------
initPackRat(
  data = qtl_regions,
  mode = "region",
  species = "human",
  genome = "hg38",
  project_dir = "qtl_opentargets_demo",
  force = TRUE
)

## ----generate_lod_helper------------------------------------------------------
# Helper function to generate realistic LOD curves
generate_lod_data <- function(regions, points_per_region = 500,
                              noise_sd = 0.3) {
  lod_list <- lapply(seq_len(nrow(regions)), function(i) {
    region <- regions[i, ]
    positions <- seq(region$start, region$end, length.out = points_per_region)

    # Create Gaussian peak centered at midpoint
    peak_pos <- (region$start + region$end) / 2
    peak_width <- (region$end - region$start) / 6

    # Base LOD from Gaussian curve + noise
    base_lod <- region$peak_lod *
      exp(-((positions - peak_pos)^2) / (2 * peak_width^2))
    lod_values <- pmax(0, base_lod + rnorm(length(positions), 0, noise_sd))

    data.table(
      chr = region$chr,
      pos = as.integer(positions),
      lod = round(lod_values, 3),
      region_id = region$region_id,
      trait = region$trait
    )
  })

  rbindlist(lod_list)
}

# Generate LOD data
set.seed(42)
lod_scan <- generate_lod_data(qtl_regions)

# Preview
head(lod_scan)

## ----add_lod_table------------------------------------------------------------
# Add scan data to project
addRatTable(
  data = lod_scan,
  table_name = "qtl_scan",
  link_type = "point",
  project_dir = "qtl_opentargets_demo"
)

## ----query_open_targets, eval=!skip_api---------------------------------------
queryOpenTargets(
  project_dir = "qtl_opentargets_demo",
  data_types = c("diseases", "constraints", "tractability"),
  disease_limit = 100,  # Limit for vignette performance
  verbose = TRUE
)

## ----skip_message_ot, eval=skip_api, echo=FALSE-------------------------------
# message("Skipping Open Targets Platform query (no internet or running on CRAN)")

## ----query_open_targets_qtl, eval=!skip_api, message=TRUE---------------------
qtl_results <- queryOpenTargetsQTL(
  project_dir = "qtl_opentargets_demo",
  study_types = c("eqtl", "pqtl"),
  include_l2g = TRUE,
  min_l2g_score = 0.5,
  verbose = TRUE
)

## ----skip_message_qtl, eval=skip_api, echo=FALSE------------------------------
# message("Skipping Open Targets QTL query (no internet or running on CRAN)")

## ----qtl_studies_preview, eval=!skip_api, echo=TRUE, results='markup'---------
if (!is.null(qtl_results$studies) && nrow(qtl_results$studies) > 0) {
  studies_dt <- qtl_results$studies
  cat("QTL Studies table:", nrow(studies_dt), "rows\n\n")

  # Select key columns for display
  display_cols <- intersect(
    c("study_id", "study_type", "biosample_name",
      "n_samples", "target_gene_symbol"),
    names(studies_dt)
  )

  cat("Head of studies table:\n")
  print(head(studies_dt[, ..display_cols], 5))

  cat("\nTail of studies table:\n")
  print(tail(studies_dt[, ..display_cols], 5))
} else {
  cat("No QTL studies found for the specified regions.\n")
}

## ----qtl_credsets_preview, eval=!skip_api, echo=TRUE, results='markup'--------
credsets_dt <- qtl_results$credible_sets
if (!is.null(credsets_dt) && nrow(credsets_dt) > 0) {
  cat("QTL Credible Sets table:", nrow(credsets_dt), "rows\n\n")

  # Select key columns for display
  display_cols <- intersect(
    c("study_id", "study_type", "chromosome", "position", "rsid",
      "beta", "pvalue_mantissa", "pvalue_exponent",
      "l2g_gene_symbol", "l2g_score"),
    names(credsets_dt)
  )

  cat("Head of credible sets table:\n")
  print(head(credsets_dt[, ..display_cols], 5))

  cat("\nTail of credible sets table:\n")
  print(tail(credsets_dt[, ..display_cols], 5))

  # Summary statistics
  cat("\n\nSummary by study type:\n")
  print(credsets_dt[, .N, by = study_type])

  if ("l2g_score" %in% names(credsets_dt)) {
    cat("\nL2G score distribution:\n")
    print(summary(credsets_dt$l2g_score))
  }
} else {
  cat("No QTL credible sets found for the specified regions.\n")
}

## ----list_tables, message=TRUE------------------------------------------------
listPackRatTables(project_dir = "qtl_opentargets_demo", full_info = TRUE)

## ----locus_zoom_znf586, eval=can_plot-----------------------------------------
# ZNF586 locus - expression QTL region
generateLocusZoomPlot(
  region_id = "ZNF586_locus",
  project_dir = "qtl_opentargets_demo",
  scan_table = "qtl_scan",
  width = 8,
  height = 5,
  highlight_genes = c("ZNF586"),
  threshold = 5,
  output_file = "ZNF586_locus_zoom.png"
)
knitr::include_graphics(
  file.path("qtl_opentargets_demo", ".locusPackRat", "output", "ZNF586_locus_zoom.png")
)

## ----locus_zoom_apip, eval=can_plot-------------------------------------------
# APIP locus - expression QTL region
generateLocusZoomPlot(
  region_id = "APIP_locus",
  project_dir = "qtl_opentargets_demo",
  scan_table = "qtl_scan",
  width = 8,
  height = 5,
  highlight_genes = c("APIP"),
  threshold = 5,
  output_file = "APIP_locus_zoom.png"
)
knitr::include_graphics(
  file.path("qtl_opentargets_demo", ".locusPackRat", "output", "APIP_locus_zoom.png")
)

## ----locus_zoom_clec19a, eval=can_plot----------------------------------------
# CLEC19A locus - protein QTL region
generateLocusZoomPlot(
  region_id = "CLEC19A_locus",
  project_dir = "qtl_opentargets_demo",
  scan_table = "qtl_scan",
  width = 8,
  height = 5,
  highlight_genes = c("CLEC19A"),
  threshold = 5,
  output_file = "CLEC19A_locus_zoom.png"
)
knitr::include_graphics(
  file.path("qtl_opentargets_demo", ".locusPackRat", "output", "CLEC19A_locus_zoom.png")
)

## ----export_excel, eval=FALSE-------------------------------------------------
# makeGeneSheet(
#   format = "excel",
#   output_file = "qtl_opentargets_results.xlsx",
#   split_by = "criteria",
#   prefix_mode = "collision",
#   split_criteria = list(
#     "All_Genes" = "TRUE",
#     "High_LOD" = "peak_lod > 7",
#     "Druggable" = "!is.na(ott_label)",
#     "Constrained" = "otc_constraintType == 'lof' & otc_oe < 0.35"
#   ),
#   project_dir = "qtl_opentargets_demo"
# )

## ----session_info-------------------------------------------------------------
sessionInfo()

## ----cleanup, include=FALSE---------------------------------------------------
# Deferred so PNGs still exist when pandoc embeds them.
reg.finalizer(
  environment(),
  function(e) unlink("qtl_opentargets_demo", recursive = TRUE),
  onexit = TRUE
)


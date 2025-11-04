#' Generate Locus Zoom Plot
#'
#' Creates a comprehensive visualization of QTL mapping results including
#' Manhattan plot, founder allele effects, overlapping QTL, and gene annotations.
#'
#' @param locus_info Single-row data.frame with columns: chr, peak_pos, start_pos,
#'   end_pos, trait, drug, max_lod
#' @param scan_data QTL scan results list containing LOD scores, positions,
#'   marker names, and allele effects matrix
#' @param threshold_data List of significance thresholds indexed by trait_drug_threshold
#' @param genes_in_locus data.table of genes within the locus region
#' @param top_genes_in_locus data.frame with gene names and highlight colors
#' @param overlapping_loci data.frame of other QTL in the region with columns:
#'   chrom, start, end, strand, traitXdrug
#' @param output_file Path for output PDF file
#' @param assembly Optional genome assembly (defaults to mm39)
#' @param plot_params Optional list of plot dimensions and positions
#'
#' @return Character path to saved PDF (invisible)
#'
#' @details
#' Generates a multi-panel plot with QTL mapping results, founder strain
#' allele effects, and gene annotations. Handles up to 8 founder strains
#' with distinct colors and multiple overlapping QTL regions.
#'
#' @importFrom plotgardener assembly pageCreate pgParams plotManhattan annoYaxis
#'   plotText annoHighlight plotSignal plotRanges plotLegend plotGenes annoGenomeLabel
#' @importFrom RColorBrewer brewer.pal
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with QTL data
#' locus <- data.frame(
#'   chr = 1, peak_pos = 105.5, start_pos = 100,
#'   end_pos = 110, trait = "HR", drug = "Ctrl", max_lod = 8.5
#' )
#'
#' generateLocusZoomPlot(
#'   locus_info = locus,
#'   scan_data = scan_data,
#'   threshold_data = threshold_data,
#'   genes_in_locus = genes,
#'   top_genes_in_locus = top_genes,
#'   overlapping_loci = overlapping,
#'   output_file = "qtl_locus.pdf"
#' )
#' }
generateLocusZoomPlot <- function(
  locus_info,
  scan_data,
  threshold_data,
  genes_in_locus,
  top_genes_in_locus,
  overlapping_loci,
  output_file,
  assembly = NULL,
  plot_params = NULL
) {

  # Check required packages
  if (!requireNamespace("plotgardener", quietly = TRUE)) {
    stop("Package 'plotgardener' is required for plotting. Please install it with:\n",
         "  BiocManager::install('plotgardener')")
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required. Please install it.")
  }

  # Define default color palettes
  # CC strain colors (8 founder strains)
  STRAIN_COLORS <- c(
    "#1B9E77",  # A/J - teal
    "#D95F02",  # C57BL/6J - orange
    "#7570B3",  # 129S1/SvImJ - purple
    "#E7298A",  # NOD/ShiLtJ - magenta
    "#66A61E",  # NZO/HILtJ - green
    "#E6AB02",  # CAST/EiJ - gold
    "#A6761D",  # PWK/PhJ - brown
    "#666666"   # WSB/EiJ - gray
  )

  # Overlapping loci colors (up to 8)
  LOCI_COLORS <- c(
    "#8DD3C7",  # light teal
    "#FFFFB3",  # light yellow
    "#BEBADA",  # light purple
    "#FB8072",  # salmon
    "#80B1D3",  # light blue
    "#FDB462",  # peach
    "#B3DE69",  # light green
    "#FCCDE5"   # light pink
  )

  # Gene highlight colors
  GENE_HIGHLIGHT_COLOR <- "#e34a33"  # red for highlighted genes
  GENE_BACKGROUND_COLOR <- "#fdbb84"  # light orange for background genes

  # Set default assembly if not provided
  if (is.null(assembly)) {
    assembly <- plotgardener::assembly(
      Genome = "mm39_GRCm39",
      TxDb   = "TxDb.Mmusculus.UCSC.mm39.knownGene",
      OrgDb  = "org.Mm.eg.db"
    )
  }

  # Set default plot parameters
  if (is.null(plot_params)) {
    plot_params <- list(
      page_width = 10.5,
      page_height = 5.5,
      x = 4.25,
      plot_width = 8,
      plot_height = 1,
      plot_y = 0.5
    )
  }

  # Define founder strains
  founders <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ",
                "NZO/H1LtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")

  # Define plot parameters
  current_chr <- as.character(locus_info$chr)

  # Calculate plot region (add 0.5 Mb padding on each side)
  # Handle both old naming (upper/lower_pos_lod_drop) and new naming (start_pos/end_pos)
  if ("start_pos" %in% names(locus_info) && "end_pos" %in% names(locus_info)) {
    region_start_mb <- locus_info$start_pos
    region_end_mb <- locus_info$end_pos
  } else {
    # Old naming - use min/max to handle either order
    region_start_mb <- min(locus_info$upper_pos_lod_drop, locus_info$lower_pos_lod_drop)
    region_end_mb <- max(locus_info$upper_pos_lod_drop, locus_info$lower_pos_lod_drop)
  }

  plot_start_bp <- max(0, floor(region_start_mb * 1e6) - 5e5)
  plot_end_bp   <- ceiling(region_end_mb * 1e6) + 5e5
  bounds_bp <- c(region_start_mb * 1e6, region_end_mb * 1e6)

  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  message("Generating locus zoom plot: ", output_file)

  # --- Start plotgardener ---
  pdf(output_file, width = plot_params$page_width, height = plot_params$page_height)
  plotgardener::pageCreate(
    width = plot_params$page_width,
    height = plot_params$page_height,
    default.units = "inches",
    showGuides = FALSE
  )

  # Set up genome parameters
  params_genome <- plotgardener::pgParams(
    assembly   = assembly,
    chrom      = paste0("chr", current_chr),
    chromstart = plot_start_bp,
    chromend   = plot_end_bp
  )

  # --- Plot miQTL LOD scores ---
  scan_key <- paste0(locus_info$trait, "_", locus_info$drug)
  current_scan_data <- scan_data[[scan_key]]

  miqtl_df_for_plot <- data.frame(
    marker = names(current_scan_data$LOD),
    chr    = as.character(current_scan_data$chr),
    pos    = current_scan_data$pos$Mb * 1e6,  # Keep as numeric (plotManhattan prefers this)
    lod    = current_scan_data$LOD,
    stringsAsFactors = FALSE
  )

  # Filter for current chromosome
  miqtl_df_for_plot <- miqtl_df_for_plot[
    miqtl_df_for_plot$chr == current_chr &
    !is.na(miqtl_df_for_plot$pos) &
    !is.na(miqtl_df_for_plot$lod),
  ]
  miqtl_df_for_plot$chrom <- paste0("chr", miqtl_df_for_plot$chr)
  miqtl_df_for_plot$p <- 10^(-miqtl_df_for_plot$lod)

  # Keep only columns plotManhattan needs (like transmute in working code)
  miqtl_df_for_plot <- miqtl_df_for_plot[, c("chrom", "pos", "p")]

  # Determine y-axis limits
  threshold_key <- paste0(locus_info$trait, "_", locus_info$drug, "_threshold")
  miqtl_threshold_val <- threshold_data[[threshold_key]]
  miqtl_ylim <- c(0, max(c(-log10(miqtl_df_for_plot$p), miqtl_threshold_val, 5), na.rm = TRUE) + 1)

  # Plot Manhattan
  miqtl_plot <- plotgardener::plotManhattan(
    data = miqtl_df_for_plot,
    params = params_genome,
    range = miqtl_ylim,
    trans = "-log10",
    sigVal = 10^(-miqtl_threshold_val),
    x = plot_params$x,
    y = plot_params$plot_y,
    width = plot_params$plot_width,
    height = plot_params$plot_height,
    just = c("center", "top"),
    fill = "#a6cee3",
    sigCol = "#1f78b4",
    sigLine = TRUE,
    baseline = TRUE,
    default.units = "inches"
  )

  # Add Y axis
  plotgardener::annoYaxis(
    plot = miqtl_plot,
    at = pretty(miqtl_ylim),
    axisLine = TRUE,
    fontsize = 8,
    main = FALSE
  )

  # Add Y axis label
  plotgardener::plotText(
    label = "LOD (miQTL)",
    x = 2 * plot_params$x + 0.1,
    y = plot_params$plot_y + plot_params$plot_height / 2,
    rot = 270,
    fontsize = 8,
    just = "center",
    default.units = "in"
  )

  # Highlight significant region
  plotgardener::annoHighlight(
    plot = miqtl_plot,
    chrom = paste0("chr", current_chr),
    chromstart = floor(min(bounds_bp)),
    chromend = ceiling(max(bounds_bp)),
    fill = "#fb9a99",
    y = plot_params$plot_y,
    height = plot_params$plot_height,
    just = c("left", "top"),
    default.units = "inches",
    alpha = 0.2,
    params = params_genome
  )

  # --- Plot allele effects ---
  pos_in_range_logical <- current_scan_data$pos$Mb * 1e6 >= plot_start_bp &
                          current_scan_data$pos$Mb * 1e6 <= plot_end_bp
  chr_logical <- current_scan_data$chr == locus_info$chr
  marker_positions_bp <- current_scan_data$pos$Mb[pos_in_range_logical & chr_logical] * 1e6

  markers <- data.frame(
    marker = current_scan_data$loci[pos_in_range_logical & chr_logical],
    start  = marker_positions_bp,
    stringsAsFactors = FALSE
  )
  markers <- markers[order(markers$start), ]

  # Check if we have markers in range
  if (nrow(markers) == 0) {
    stop("No markers found in plot region for allele effects. Check plot boundaries.")
  }

  allele_effects_transposed <- t(current_scan_data$allele.effects)
  allele_effects_transposed <- as.data.frame(allele_effects_transposed)
  allele_effects_transposed$marker <- rownames(allele_effects_transposed)

  num_strains <- ncol(allele_effects_transposed) - 1
  chromosome_name <- paste0("chr", locus_info$chr)
  founder_strains <- rownames(current_scan_data$allele.effects)

  # Create plot data for each strain
  plot_data_list <- lapply(1:num_strains, function(strain_idx) {
    curr_strain <- founder_strains[strain_idx]
    temp_df <- allele_effects_transposed[, c("marker", curr_strain), drop = FALSE]
    temp_df <- temp_df[temp_df$marker %in% markers$marker, ]
    temp_df <- merge(temp_df, markers, by = "marker")
    temp_df$chrom <- chromosome_name
    temp_df <- temp_df[order(temp_df$start), ]
    colnames(temp_df)[2] <- "score"

    # Make ranges continuous
    temp_df$end <- c(temp_df$start[2:nrow(temp_df)] - 1L,
                     temp_df$start[nrow(temp_df)] + 1)
    temp_df <- temp_df[, c("chrom", "start", "end", "score")]


    return(temp_df)
  })
  names(plot_data_list) <- founder_strains

  # Get max allele effect for y-axis scaling
  max_allele_effect <- max(abs(unlist(lapply(plot_data_list, function(x) x$score))))
  max_allele_effect <- round(max_allele_effect * 1.1, digits = 2)

  # Use predefined strain colors
  strain_colors <- rep(STRAIN_COLORS, length.out = num_strains)

  # Plot first strain with axes
  signal_y_pos <- plot_params$plot_y + plot_params$plot_height + 0.2
  signalPlots <- c()
signalPlots[[1]] <- plotgardener::plotSignal(
    data = plot_data_list[[1]],
    params = params_genome,
    range = c(-max_allele_effect, max_allele_effect),
    linecolor = strain_colors[1],
    fill = NA,
    x = plot_params$x,
    y = signal_y_pos,
    width = plot_params$plot_width,
    height = plot_params$plot_height,
    just = c("center", "top"),
    default.units = "inches",
    baseline = TRUE,
    baseline.color = "grey"
  )

  plotgardener::annoYaxis(
    plot =   signalPlots[[1]],
    at = c(-max_allele_effect, 0, max_allele_effect),
    axisLine = TRUE,
    fontsize = 8,
    main = FALSE
  )

  # Plot remaining strains
  for (strain_idx in 2:8) {
    signalPlots[[strain_idx]] <- plotgardener::plotSignal(
      data = plot_data_list[[strain_idx]],
      params = params_genome,
      range = c(-max_allele_effect, max_allele_effect),
      linecolor = strain_colors[strain_idx],
      fill = NA,
      x = plot_params$x,
      y = signal_y_pos,
      width = plot_params$plot_width,
      height = plot_params$plot_height,
      just = c("center", "top"),
      default.units = "inches",
      baseline = TRUE,
      baseline.color = "grey"
    )
  }

  # Add founder effects label
  plotgardener::plotText(
    label = "Founder Effects",
    x = 2 * plot_params$x + 0.1,
    y = signal_y_pos + plot_params$plot_height / 2,
    rot = 270,
    fontsize = 8,
    just = c("center", "center"),
    default.units = "in"
  )

  # --- Plot overlapping loci ranges ---
  # Get unique loci and assign colors
  unique_loci <- unique(overlapping_loci$traitXdrug)
  n_loci <- length(unique_loci)
  loci_color_map <- LOCI_COLORS[1:min(n_loci, length(LOCI_COLORS))]

  # Set ranges_y_pos for both single and multiple loci cases
  ranges_y_pos <- signal_y_pos + plot_params$plot_height + 0.2

  if (nrow(overlapping_loci) > 1) {

    # Create color palette function for plotRanges
    loci_palette <- function(n) {
      if (n <= length(LOCI_COLORS)) {
        return(LOCI_COLORS[1:n])
      } else {
        return(grDevices::colorRampPalette(LOCI_COLORS)(n))
      }
    }

    plotgardener::plotRanges(
      data = overlapping_loci,
      params = params_genome,
      order = "random",
      fill = plotgardener::colorby("traitXdrug", palette = loci_palette),
      x = plot_params$x,
      y = ranges_y_pos,
      width = plot_params$plot_width,
      height = plot_params$plot_height,
      just = c("center", "top"),
      default.units = "inches"
    )

    plotgardener::plotText(
      label = "Other Loci",
      x = 2 * plot_params$x + 0.1,
      y = ranges_y_pos + plot_params$plot_height / 2,
      rot = 270,
      fontsize = 8,
      just = "center",
      default.units = "in"
    )
  }

  # --- Plot genes ---
  # Set gene position based on whether we plotted overlapping loci
  if (nrow(overlapping_loci) > 1) {
    genes_y_pos <- ranges_y_pos + plot_params$plot_height + 0.2
  } else {
    genes_y_pos <- signal_y_pos + plot_params$plot_height + 0.2
  }

  # Extract gene names for geneOrder (needs character vector)
  gene_order <- if(nrow(top_genes_in_locus) > 0) {
    top_genes_in_locus$gene
  } else {
    NULL
  }

  # Try to plot genes, but handle errors gracefully
  gene_plot <- tryCatch({
    plotgardener::plotGenes(
      params = params_genome,
      x = plot_params$x,
      y = genes_y_pos,
      width = plot_params$plot_width,
      height = 1,
      just = c("center", "top"),
      default.units = "inches",
      geneOrder = gene_order,
      fontsize = 6,
      geneHighlights = top_genes_in_locus,
      geneBackground = GENE_BACKGROUND_COLOR
    )
  }, error = function(e) {
    # If gene highlighting fails, try without it
    message("Note: Gene highlighting failed, plotting without highlights")
    plotgardener::plotGenes(
      params = params_genome,
      x = plot_params$x,
      y = genes_y_pos,
      width = plot_params$plot_width,
      height = 1,
      just = c("center", "top"),
      default.units = "inches",
      fontsize = 6
    )
  })

  # --- Add genome label ---
  plotgardener::annoGenomeLabel(
    plot = gene_plot,
    params = params_genome,
    x = plot_params$x,
    y = genes_y_pos + 1,
    scale = "Mb",
    fontsize = 10,
    just = c("center", "top"),
    default.units = "inches"
  )

  # --- Add legends ---
  # Founder strain legend
  plotgardener::plotLegend(
    legend = founders,
    fill = strain_colors,
    border = FALSE,
    x = 8.7,
    y = signal_y_pos,
    width = 2,
    height = 1,
    fontsize = 7,
    just = c("left", "center"),
    orientation = "v",
    default.units = "inches"
  )

  # Overlapping loci legend (if present)
  if (nrow(overlapping_loci) > 0) {
    plotgardener::plotLegend(
      legend = unique_loci,
      fill = loci_color_map,
      border = FALSE,
      x = 8.7,
      y = ranges_y_pos,
      width = 2,
      height = 1,
      fontsize = 7,
      just = c("left", "center"),
      orientation = "v",
      default.units = "inches"
    )
  }

  # Top genes legend
  plotgardener::plotLegend(
    legend = c("Protein coding, human orthologue, and >5 CPM in NRVMs"),
    fill = GENE_HIGHLIGHT_COLOR,
    border = FALSE,
    x = 0.25,
    y = genes_y_pos - 0.1,
    width = 5,
    height = 0.15,
    fontsize = 7,
    just = c("left", "top"),
    orientation = "v",
    default.units = "inches"
  )

  # --- Add title ---
  plotgardener::plotText(
    label = paste("QTL:", locus_info$trait, "(", locus_info$drug, ") - Chr",
                  current_chr, "Peak:", round(locus_info$peak_pos, 2), "Mb"),
    x = plot_params$x,
    y = 0.1,
    just = c("center", "top"),
    fontface = "bold",
    fontsize = 12,
    default.units = "inches"
  )

  dev.off()

  return(invisible(output_file))
}


#' Generate Multi-Sheet Gene Information Excel Workbook
#'
#' Creates an Excel file with comprehensive gene annotations across multiple sheets:
#' Sheet 1 (AllGenesInCluster) contains main gene summary with locus presence,
#' Sheet 2 (AllMousePhenotypes) has MGI phenotype data, and Sheet 3 (AllHumanDiseases)
#' contains Open Targets disease associations.
#'
#' @param genes_in_locus data.table with genes in the locus region. Required columns:
#'   \code{mouse_ensembl_id}
#' @param loci_info data.frame with ALL loci in this cluster. Required columns:
#'   \describe{
#'     \item{chr}{Chromosome}
#'     \item{upper_pos_lod_drop}{Upper LOD drop boundary in Mb}
#'     \item{lower_pos_lod_drop}{Lower LOD drop boundary in Mb}
#'     \item{peak_pos}{Peak position in Mb}
#'     \item{trait}{Trait name}
#'     \item{drug}{Treatment condition}
#'   }
#' @param merged_gene_info data.table with gene annotations. Required columns:
#'   \describe{
#'     \item{mouse_ensembl_id}{Ensembl gene ID}
#'     \item{mouse_gene_symbol}{Gene symbol}
#'     \item{n_trait_drug}{Number of associated traits}
#'     \item{trait_drug}{Associated trait:drug combinations}
#'     \item{avgenes_cpm}{Average CPM in control NRVMs}
#'     \item{disease}{Human disease summary string}
#'     \item{ontology}{Mouse phenotype ontology terms}
#'     \item{pubmedID}{PubMed IDs for phenotype evidence}
#'   }
#' @param genes_mouse data.table with gene coordinates. Required columns:
#'   \code{mouse_ensembl_id}, \code{chr}, \code{start_bp}, \code{end_bp}
#' @param ortho_mouse2h data.table with orthology. Required columns:
#'   \code{mouse_ensembl_id}, \code{human_ensembl_id}, \code{human_gene_symbol}
#' @param mouse_pheno data.table with MGI phenotypes. Required columns:
#'   \code{mouse_gene_symbol}, \code{OntologyAnnotation.ontologyTerm.identifier},
#'   \code{OntologyAnnotation.ontologyTerm.name},
#'   \code{OntologyAnnotation.evidence.publications.pubMedId},
#'   \code{OntologyAnnotation.evidence.comments.description}
#' @param associations data.table with human disease associations. Required columns:
#'   \code{human_ensembl_id}, \code{symbol}, \code{disease_id}, \code{disease_name},
#'   \code{association_score}
#' @param rna_info data.table with bulk RNA-seq expression. Must have columns:
#'   \code{Drug}, \code{Sex}, and one column per gene symbol with expression values
#' @param founder_mutations data.table with CC founder mutations (optional).
#'   If provided, must have columns: \code{Gene}, \code{Mutations}
#' @param output_file Character path to save Excel workbook
#'
#' @return data.table of main gene summary table (invisible)
#'
#' @details
#' The AllGenesInCluster sheet includes:
#' \itemize{
#'   \item Gene identifiers and coordinates
#'   \item Expression levels (CPM and VST by drug/sex)
#'   \item Number and names of associated traits
#'   \item Boolean columns indicating presence in each locus
#'   \item Disease and phenotype summaries
#'   \item Founder mutations (if provided)
#' }
#'
#' All sheets have frozen header rows, auto-filter, and auto-sized columns.
#'
#' @importFrom openxlsx createWorkbook addWorksheet writeData writeDataTable
#'   createStyle addStyle setColWidths freezePane addFilter saveWorkbook
#' @importFrom data.table data.table setnames merge.data.table setDT :=
#' @importFrom dplyr left_join group_by summarize across pivot_longer unite
#'   pivot_wider select filter
#' @export
#'
#' @examples
#' \dontrun{
#' # Load required data
#' genes_mouse <- fread("data/processed/joinLoci/relational_tables/genes_mouse.csv")
#' merged_info <- fread("data/processed/joinLoci/geneTables/multTrait_summary.csv")
#' ortho <- fread("data/processed/joinLoci/relational_tables/orthology.csv")
#' pheno <- fread("data/processed/joinLoci/relational_tables/mouseGenePhenotypes.csv")
#' assoc <- fread("data/processed/joinLoci/relational_tables/associations.csv")
#' rna <- fread("data/processed/joinLoci/bulk_exp/VST_Info.csv", drop = 1)
#'
#' # Define loci
#' loci <- data.frame(
#'   chr = 1, upper_pos_lod_drop = 100, lower_pos_lod_drop = 110,
#'   peak_pos = 105, trait = "HR", drug = "Ctrl"
#' )
#'
#' # Get genes in region
#' genes_in_locus <- genes_mouse[chr == 1 & start_bp < 110e6 & end_bp > 100e6]
#'
#' # Generate Excel file
#' generateGeneInfoExcel(
#'   genes_in_locus = genes_in_locus,
#'   loci_info = loci,
#'   merged_gene_info = merged_info,
#'   genes_mouse = genes_mouse,
#'   ortho_mouse2h = ortho,
#'   mouse_pheno = pheno,
#'   associations = assoc,
#'   rna_info = rna,
#'   founder_mutations = NULL,
#'   output_file = "gene_info.xlsx"
#' )
#' }
generateGeneInfoExcel <- function(
  genes_in_locus,
  loci_info,
  merged_gene_info,
  genes_mouse,
  ortho_mouse2h,
  mouse_pheno,
  associations,
  rna_info,
  founder_mutations = NULL,
  output_file
) {

  # Check if genes in locus
  if (nrow(genes_in_locus) == 0) {
    message("No genes found in the specified locus for Excel generation.")
    return(invisible(NULL))
  }

  # Create main gene summary
  ensembl_ids_in_locus <- genes_in_locus$mouse_ensembl_id
  locus_gene_summary <- merged_gene_info[mouse_ensembl_id %in% ensembl_ids_in_locus, ]

  # Add founder mutations if provided
  if (!is.null(founder_mutations) && nrow(founder_mutations) > 0) {
    founder_muts_gene_level <- founder_mutations[
      Gene %in% locus_gene_summary$mouse_gene_symbol,
      .(mouse_gene_symbol = Gene, founder_gene_mutations = Mutations)
    ]

    if (nrow(founder_muts_gene_level) > 0) {
      locus_gene_summary <- merge(locus_gene_summary, founder_muts_gene_level,
                                  by = "mouse_gene_symbol", all.x = TRUE)
    } else {
      locus_gene_summary[, founder_gene_mutations := NA_character_]
    }
  } else {
    locus_gene_summary[, founder_gene_mutations := NA_character_]
  }

  # Create main genes table
  main_genes_table <- locus_gene_summary[, .(
    `Mouse Gene Symbol` = mouse_gene_symbol,
    `Mouse Ensembl ID` = mouse_ensembl_id,
    `# Traits Associated (miQTL)` = n_trait_drug,
    `Associated Traits (Drug)` = trait_drug,
    `Avg NRVM CPM (Ctrl)` = round(avgenes_cpm, 2),
    `Human Disease Summary` = disease,
    `Mouse Phenotype Summary (MGI)` = ontology,
    `Supporting Publications (MGI)` = pubmedID,
    `Founder Gene Mutations (CC)` = founder_gene_mutations
  )]

  # Add gene coordinates
  gene_coords <- genes_mouse[
    mouse_ensembl_id %in% main_genes_table$`Mouse Ensembl ID`,
    .(mouse_ensembl_id, chr, start_bp, end_bp)
  ]
  setnames(gene_coords, "mouse_ensembl_id", "Mouse Ensembl ID")
  main_genes_table <- merge(main_genes_table, gene_coords, by = "Mouse Ensembl ID", all.x = TRUE)

  # Get genes with expression data
  measured_genes <- main_genes_table[`Mouse Gene Symbol` %in% colnames(rna_info)]

  if (nrow(measured_genes) == 0) {
    message("No measured genes found in expression data.")
    return(invisible(NA))
  }

  # Calculate mean expression by Drug and Sex
  gene_cols <- measured_genes$`Mouse Gene Symbol`
  means <- rna_info %>%
    dplyr::select(Drug, Sex, dplyr::all_of(gene_cols)) %>%
    dplyr::group_by(Drug, Sex) %>%
    dplyr::summarize(dplyr::across(dplyr::all_of(gene_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  # Reshape to get Drug_Sex combinations as columns
  means_long <- means %>%
    tidyr::pivot_longer(
      cols = -c(Drug, Sex),
      names_to = "Gene",
      values_to = "Mean"
    ) %>%
    dplyr::mutate(Mean = round(Mean, digits = 1)) %>%
    tidyr::unite("Drug_Sex", Drug, Sex, sep = "_") %>%
    tidyr::pivot_wider(
      id_cols = Gene,
      names_from = Drug_Sex,
      values_from = Mean,
      names_prefix = "Ave_Exp_"
    )

  # Join with main genes table
  main_genes_table <- main_genes_table %>%
    dplyr::left_join(means_long, by = c("Mouse Gene Symbol" = "Gene"))

  # Set initial column order
  initial_cols <- c(
    "Mouse Gene Symbol", "Mouse Ensembl ID", "chr", "start_bp", "end_bp",
    "Avg NRVM CPM (Ctrl)", "Ave_Exp_Ctrl_F", "Ave_Exp_Ctrl_M",
    "Ave_Exp_Iso_F", "Ave_Exp_Iso_M", "# Traits Associated (miQTL)"
  )

  # --- Add locus presence columns (T/F) for each locus in the cluster ---
  loci_info$locus_id <- paste0(
    loci_info$trait, "_", loci_info$drug,
    "_chr", loci_info$chr, "_", round(loci_info$peak_pos, 2), "Mb"
  )

  # For each locus, create a column indicating if each gene is present
  for (i in 1:nrow(loci_info)) {
    current_locus <- loci_info[i, ]
    locus_col_name <- paste0("In_", current_locus$locus_id)

    # Get locus boundaries
    locus_start_bp <- min(current_locus$upper_pos_lod_drop, current_locus$lower_pos_lod_drop) * 1e6
    locus_end_bp <- max(current_locus$upper_pos_lod_drop, current_locus$lower_pos_lod_drop) * 1e6

    # Check which genes are in this locus
    data.table::setDT(main_genes_table)
    main_genes_table[, (locus_col_name) := "No"]
    genes_in_this_locus <- main_genes_table[
      chr == current_locus$chr & start_bp < locus_end_bp & end_bp > locus_start_bp,
      `Mouse Ensembl ID`
    ]

    # Set to TRUE for genes in this locus
    main_genes_table[`Mouse Ensembl ID` %in% genes_in_this_locus, (locus_col_name) := "Yes"]
  }

  # Create final column order
  locus_cols <- grep("^In_", names(main_genes_table), value = TRUE)
  remaining_cols <- setdiff(names(main_genes_table), c(initial_cols, locus_cols))
  col_order <- c(initial_cols, locus_cols, remaining_cols)
  main_genes_table <- main_genes_table[, ..col_order]

  # --- Create Excel Workbook ---
  wb <- openxlsx::createWorkbook()

  # Add main summary sheet
  openxlsx::addWorksheet(wb, "AllGenesInCluster")
  openxlsx::writeDataTable(
    wb,
    sheet = "AllGenesInCluster",
    x = main_genes_table,
    tableName = "AllGenesTbl",
    withFilter = TRUE
  )

  # Wrap text style
  wrapStyle <- openxlsx::createStyle(wrapText = TRUE)
  openxlsx::addStyle(
    wb,
    sheet = "AllGenesInCluster",
    style = wrapStyle,
    rows = 1:(nrow(main_genes_table) + 1),
    cols = 1:ncol(main_genes_table),
    gridExpand = TRUE
  )

  # Auto-fit columns
  openxlsx::setColWidths(
    wb,
    sheet = "AllGenesInCluster",
    cols = 1:ncol(main_genes_table),
    widths = "auto"
  )

  # Freeze header row
  openxlsx::freezePane(
    wb,
    sheet = "AllGenesInCluster",
    firstActiveRow = 2
  )

  # --- Create mouse phenotypes sheet ---
  all_mouse_pheno <- data.table::data.table()
  all_human_disease <- data.table::data.table()

  for (idx in 1:nrow(main_genes_table)) {
    current_gene_symbol <- main_genes_table$`Mouse Gene Symbol`[idx]
    current_mouse_ensembl_id <- main_genes_table$`Mouse Ensembl ID`[idx]

    # Mouse phenotypes
    gene_mouse_pheno_data <- mouse_pheno[
      mouse_gene_symbol == current_gene_symbol,
      .(
        Gene = mouse_gene_symbol,
        OntologyTermID = `OntologyAnnotation.ontologyTerm.identifier`,
        OntologyTermName = `OntologyAnnotation.ontologyTerm.name`,
        PubMedID = `OntologyAnnotation.evidence.publications.pubMedId`,
        Description = `OntologyAnnotation.evidence.comments.description`
      )
    ]

    if (nrow(gene_mouse_pheno_data) > 0) {
      all_mouse_pheno <- data.table::rbindlist(list(all_mouse_pheno, gene_mouse_pheno_data), fill = TRUE)
    } else {
      all_mouse_pheno <- data.table::rbindlist(
        list(
          all_mouse_pheno,
          data.table::data.table(
            Gene = current_gene_symbol,
            OntologyTermID = NA_character_,
            OntologyTermName = "No phenotypes found",
            PubMedID = NA_character_,
            Description = NA_character_
          )
        ),
        fill = TRUE
      )
    }

    # Human disease associations
    human_ortho_info <- ortho_mouse2h[mouse_ensembl_id == current_mouse_ensembl_id, ]

    if (nrow(human_ortho_info) > 0 && human_ortho_info$human_ensembl_id[1] != "") {
      current_human_ensembl_id <- human_ortho_info$human_ensembl_id[1]

      gene_human_disease_data <- associations[
        human_ensembl_id == current_human_ensembl_id,
        .(
          MouseGene = current_gene_symbol,
          human_ensembl_id,
          human_gene_symbol = symbol,
          disease_id,
          disease_name,
          association_score
        )
      ]

      if (nrow(gene_human_disease_data) > 0) {
        all_human_disease <- data.table::rbindlist(list(all_human_disease, gene_human_disease_data), fill = TRUE)
      }
    }
  }

  # Add consolidated sheets
  if (nrow(all_mouse_pheno) > 0) {
    openxlsx::addWorksheet(wb, "AllMousePhenotypes")
    openxlsx::writeData(wb, "AllMousePhenotypes", all_mouse_pheno)
    openxlsx::addFilter(wb, "AllMousePhenotypes", rows = 1, cols = 1:ncol(all_mouse_pheno))

    openxlsx::addStyle(
      wb,
      sheet = "AllMousePhenotypes",
      style = wrapStyle,
      rows = 1:(nrow(all_mouse_pheno) + 1),
      cols = 1:ncol(all_mouse_pheno),
      gridExpand = TRUE
    )

    openxlsx::setColWidths(
      wb,
      sheet = "AllMousePhenotypes",
      cols = 1:ncol(all_mouse_pheno),
      widths = "auto"
    )

    openxlsx::freezePane(
      wb,
      sheet = "AllMousePhenotypes",
      firstActiveRow = 2
    )
  }

  if (nrow(all_human_disease) > 0) {
    openxlsx::addWorksheet(wb, "AllHumanDiseases")
    openxlsx::writeData(wb, "AllHumanDiseases", all_human_disease)
    openxlsx::addFilter(wb, "AllHumanDiseases", rows = 1, cols = 1:ncol(all_human_disease))

    openxlsx::addStyle(
      wb,
      sheet = "AllHumanDiseases",
      style = wrapStyle,
      rows = 1:(nrow(all_human_disease) + 1),
      cols = 1:ncol(all_human_disease),
      gridExpand = TRUE
    )

    openxlsx::setColWidths(
      wb,
      sheet = "AllHumanDiseases",
      cols = 1:ncol(all_human_disease),
      widths = "auto"
    )

    openxlsx::freezePane(
      wb,
      sheet = "AllHumanDiseases",
      firstActiveRow = 2
    )
  }

  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save workbook
  message("Generating gene information Excel file: ", output_file)
  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)

  return(invisible(main_genes_table))
}


#' Generate Complete Locus Packet
#'
#' Master function that creates a complete packet directory for a locus cluster,
#' including Excel workbook, locus zoom plots, README summary, and SNP table.
#'
#' @param locus_cluster data.frame with one or more overlapping loci. Required columns:
#'   \describe{
#'     \item{chr}{Chromosome}
#'     \item{upper_pos_lod_drop}{Upper LOD drop boundary in Mb}
#'     \item{lower_pos_lod_drop}{Lower LOD drop boundary in Mb}
#'     \item{peak_pos}{Peak position in Mb}
#'     \item{max_lod}{Maximum LOD score}
#'     \item{trait}{Trait name}
#'     \item{drug}{Treatment condition}
#'   }
#' @param input_path Character path to data folder (e.g., "data/").
#'   Used as base path for constructing default file paths.
#' @param output_path Character path to save packets (default: "results/qtl_packets/")
#' @param scan_data List with QTL scan results (from miqtl or similar)
#' @param threshold_data List with significance thresholds per trait
#' @param assembly Genome assembly object from \code{plotgardener::assembly()}.
#'   If NULL, defaults to mm39 (GRCm39).
#' @param heart_pattern Character vector of heart-related terms for phenotype filtering.
#'   If NULL, uses default cardiac terms. Set to empty string to skip filtering.
#' @param genes_mouse_file Path to genes CSV. Default: \code{input_path/processed/joinLoci/relational_tables/genes_mouse.csv}
#' @param orthology_file Path to orthology CSV. Default: \code{input_path/processed/joinLoci/relational_tables/orthology.csv}
#' @param trait_loci_file Path to trait loci CSV. Default: \code{input_path/processed/joinLoci/relational_tables/traitLoci.csv}
#' @param associations_file Path to disease associations CSV. Default: \code{input_path/processed/joinLoci/relational_tables/associations.csv}
#' @param mouse_pheno_file Path to MGI phenotypes CSV. Default: \code{input_path/processed/joinLoci/relational_tables/mouseGenePhenotypes.csv}
#' @param merged_gene_info_file Path to merged gene summary CSV. Default: \code{input_path/processed/joinLoci/geneTables/multTrait_cis-eQTL_nrvmExp.csv}
#' @param rna_info_file Path to expression data CSV. Default: \code{input_path/processed/joinLoci/bulk_exp/rna_expression.csv}
#' @param founder_mutations_file Path to CC founder mutations CSV (optional). Default: \code{input_path/processed/joinLoci/relational_tables/ccVariants/gene_mutations.csv}
#' @param founder_snps_file Path to CC founder SNPs CSV (optional). Default: \code{input_path/processed/joinLoci/relational_tables/ccVariants/snp_mutations.csv}
#'
#' @return Character path to output directory (invisible)
#'
#' @details
#' This function orchestrates the complete packet generation workflow:
#' \enumerate{
#'   \item Loads all annotation tables from specified paths
#'   \item Identifies genes in the locus region
#'   \item Creates output directory named by locus coordinates
#'   \item Generates multi-sheet Excel workbook with gene annotations
#'   \item Generates locus zoom plots for each QTL in the cluster
#'   \item Creates README summary highlighting cardiac-related genes
#'   \item Exports founder SNP table (if SNP data available)
#' }
#'
#' \strong{Required Data Files and Schemas:}
#'
#' All paths can be customized via parameters. Defaults assume standard structure:
#'
#' \itemize{
#'   \item \strong{genes_mouse.csv}: Mouse gene coordinates
#'     \itemize{
#'       \item Required columns: \code{mouse_ensembl_id}, \code{mouse_gene_symbol}, \code{chr}, \code{start_bp}, \code{end_bp}
#'     }
#'   \item \strong{orthology.csv}: Mouse-human gene orthology
#'     \itemize{
#'       \item Required columns: \code{mouse_ensembl_id}, \code{human_ensembl_id}, \code{human_gene_symbol}
#'     }
#'   \item \strong{associations.csv}: Human gene-disease associations (from Open Targets)
#'     \itemize{
#'       \item Required columns: \code{human_ensembl_id}, \code{symbol}, \code{disease_id}, \code{disease_name}, \code{association_score}
#'     }
#'   \item \strong{mouseGenePhenotypes.csv}: MGI phenotype annotations
#'     \itemize{
#'       \item Required columns: \code{mouse_gene_symbol}, MGI ontology terms, publications
#'     }
#'   \item \strong{multTrait_cis-eQTL_nrvmExp.csv}: Merged gene summary
#'     \itemize{
#'       \item Required columns: \code{mouse_ensembl_id}, \code{mouse_gene_symbol}, \code{trait_drug}, \code{avgenes_cpm}
#'     }
#'   \item \strong{rna_expression.csv}: Bulk RNA-seq VST values
#'     \itemize{
#'       \item Required columns: \code{Drug}, \code{Sex}, plus columns named by gene symbols
#'     }
#'   \item \strong{gene_mutations.csv} (optional): CC founder gene mutations
#'     \itemize{
#'       \item Required columns: \code{Gene}, \code{Mutations}
#'     }
#'   \item \strong{snp_mutations.csv} (optional): CC founder SNP mutations
#'     \itemize{
#'       \item Required columns: \code{SNP}, \code{rsID}, \code{CHR}, \code{mm39}, \code{Major}, \code{Minor}, \code{MAF}, \code{SNPStrains}
#'     }
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with default paths
#' generateLocusPacket(
#'   locus_cluster = my_locus_df,
#'   input_path = "data/",
#'   scan_data = scan_results,
#'   threshold_data = thresholds
#' )
#'
#' # Custom RNA-seq file (e.g., dated file)
#' generateLocusPacket(
#'   locus_cluster = my_locus_df,
#'   input_path = "data/",
#'   scan_data = scan_results,
#'   threshold_data = thresholds,
#'   rna_info_file = "data/processed/bulk_exp/VST_Info_250429.csv"
#' )
#' }
generateLocusPacket <- function(
  locus_cluster,
  input_path = "data/",
  output_path = "results/qtl_packets/",
  scan_data,
  threshold_data,
  assembly = NULL,
  heart_pattern = NULL,
  genes_mouse_file = NULL,
  orthology_file = NULL,
  trait_loci_file = NULL,
  associations_file = NULL,
  mouse_pheno_file = NULL,
  merged_gene_info_file = NULL,
  rna_info_file = NULL,
  founder_mutations_file = NULL,
  founder_snps_file = NULL
) {

  # Set default assembly
  if (is.null(assembly)) {
    assembly <- plotgardener::assembly(
      Genome = "mm39_GRCm39",
      TxDb   = "TxDb.Mmusculus.UCSC.mm39.knownGene",
      OrgDb  = "org.Mm.eg.db"
    )
  }

  # Set default heart pattern
  if (is.null(heart_pattern)) {
    heart_pattern <- .getDefaultHeartPattern()
  }

  message("\n--- Generating Locus Packet ---")

  # Set default file paths
  base_relational <- file.path(input_path, "processed/joinLoci/relational_tables")

  if (is.null(genes_mouse_file)) {
    genes_mouse_file <- file.path(base_relational, "genes_mouse.csv")
  }
  if (is.null(orthology_file)) {
    orthology_file <- file.path(base_relational, "orthology.csv")
  }
  if (is.null(trait_loci_file)) {
    trait_loci_file <- file.path(base_relational, "traitLoci.csv")
  }
  if (is.null(associations_file)) {
    associations_file <- file.path(base_relational, "associations.csv")
  }
  if (is.null(mouse_pheno_file)) {
    mouse_pheno_file <- file.path(base_relational, "mouseGenePhenotypes.csv")
  }
  if (is.null(merged_gene_info_file)) {
    merged_gene_info_file <- file.path(input_path, "processed/joinLoci/geneTables/multTrait_cis-eQTL_nrvmExp.csv")
  }
  if (is.null(rna_info_file)) {
    rna_info_file <- file.path(input_path, "processed/joinLoci/bulk_exp/rna_expression.csv")
  }
  if (is.null(founder_mutations_file)) {
    founder_mutations_file <- file.path(base_relational, "ccVariants/gene_mutations.csv")
  }
  if (is.null(founder_snps_file)) {
    founder_snps_file <- file.path(base_relational, "ccVariants/snp_mutations.csv")
  }

  # Load all relational tables
  message("Loading annotation tables...")
  genes_mouse <- data.table::fread(genes_mouse_file)
  ortho_mouse2h <- data.table::fread(orthology_file)
  trait_loci <- data.table::fread(trait_loci_file)
  associations <- data.table::fread(associations_file)
  mouse_pheno <- data.table::fread(mouse_pheno_file)

  # Load merged gene info
  merged_gene_info <- data.table::fread(merged_gene_info_file)

  # Load expression data (drop first column if it's row numbers)
  rna_info <- data.table::fread(rna_info_file, drop = 1)

  # Load founder mutations (optional)
  if (file.exists(founder_mutations_file)) {
    founder_mutations <- data.table::fread(founder_mutations_file)
  } else {
    message("Founder mutations file not found: ", founder_mutations_file)
    founder_mutations <- NULL
  }

  # Load SNP data (optional)
  if (file.exists(founder_snps_file)) {
    founder_snps <- data.table::fread(founder_snps_file)
  } else {
    message("Founder SNPs file not found: ", founder_snps_file)
    founder_snps <- NULL
  }

  # Determine locus boundaries
  locus_start_bp <- min(locus_cluster$upper_pos_lod_drop, locus_cluster$lower_pos_lod_drop) * 1e6
  locus_end_bp <- max(locus_cluster$upper_pos_lod_drop, locus_cluster$lower_pos_lod_drop) * 1e6
  locus_chr <- locus_cluster$chr[1]

  # Create output directory
  locus_name <- paste0(
    "locus_chr", locus_chr, "_",
    round(min(locus_cluster$upper_pos_lod_drop), 0), "-",
    round(max(locus_cluster$lower_pos_lod_drop), 0), "Mb"
  )
  output_dir <- file.path(output_path, locus_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  message("Output directory: ", output_dir)

  # Get genes in locus
  genes_in_locus <- .getGenesInRegion(
    target_chr = locus_chr,
    target_start = locus_start_bp,
    target_end = locus_end_bp,
    gene_data = genes_mouse
  )

  message("Found ", nrow(genes_in_locus), " genes in locus region")

  if (nrow(genes_in_locus) == 0) {
    warning("No genes found in locus, skipping packet generation")
    return(invisible(output_dir))
  }

  # Get top genes for highlighting
  top_genes <- merged_gene_info[
    mouse_gene_symbol %in% genes_in_locus$mouse_gene_symbol,
    mouse_gene_symbol
  ]
  if (length(top_genes) < 1 || all(is.na(top_genes))) {
    top_genes_df <- data.frame(gene = NA, color = NA)
  } else {
    top_genes_df <- data.frame(gene = top_genes, color = "#e34a33", stringsAsFactors = FALSE)
  }

  # Create overlapping loci dataframe for plotting
  overlaps <- locus_cluster
  overlaps$chrom <- paste0("chr", overlaps$chr)
  overlaps$start <- overlaps$upper_pos_lod_drop * 1e6
  overlaps$end <- overlaps$lower_pos_lod_drop * 1e6
  overlaps$strand <- "-"
  overlaps$traitXdrug <- paste0(overlaps$trait, ": ", overlaps$drug)
  overlaps <- overlaps[, c("chrom", "start", "end", "strand", "trait", "drug", "traitXdrug")]

  # Generate Excel workbook
  excel_file <- file.path(
    output_dir,
    paste0("gene_info_cluster_chr", locus_chr, "_",
           round(min(locus_cluster$upper_pos_lod_drop), 0), "-",
           round(max(locus_cluster$lower_pos_lod_drop), 0), "Mb.xlsx")
  )

  gene_table <- generateGeneInfoExcel(
    genes_in_locus = genes_in_locus,
    loci_info = locus_cluster,
    merged_gene_info = merged_gene_info,
    genes_mouse = genes_mouse,
    ortho_mouse2h = ortho_mouse2h,
    mouse_pheno = mouse_pheno,
    associations = associations,
    rna_info = rna_info,
    founder_mutations = founder_mutations,
    output_file = excel_file
  )

  # If no gene table generated (e.g., no measured genes), create minimal README and exit
  if (is.null(gene_table) || all(is.na(gene_table))) {
    summary_text <- paste(
      "# LOCUS PACKET SUMMARY\n",
      "Locus:", locus_name,
      "\nAssociated Trait(s):",
      paste(unique(paste(locus_cluster$trait, locus_cluster$drug, sep = " - ")), collapse = ", "),
      "\nChromosome:", locus_chr,
      "\nRegion:", round(locus_start_bp / 1e6, 2), "-", round(locus_end_bp / 1e6, 2), "Mb",
      "\n\nWARNING: No compelling genes found with expression data. Locus may be erroneous."
    )
    writeLines(summary_text, file.path(output_dir, "README_summary.txt"))
    return(invisible(output_dir))
  }

  # Generate locus zoom plots
  zoom_dir <- file.path(output_dir, "zoomPlots")
  if (!dir.exists(zoom_dir)) {
    dir.create(zoom_dir)
  }

  for (i in 1:nrow(locus_cluster)) {
    current_locus <- locus_cluster[i, ]
    plot_file <- file.path(
      zoom_dir,
      paste0("locus_zoom_", current_locus$trait, "_", current_locus$drug,
             "_chr", current_locus$chr, "_", round(current_locus$peak_pos, 2), "Mb.pdf")
    )

    generateLocusZoomPlot(
      locus_info = current_locus,
      scan_data = scan_data,
      threshold_data = threshold_data,
      genes_in_locus = genes_in_locus,
      top_genes_in_locus = top_genes_df,
      overlapping_loci = overlaps,
      output_file = plot_file,
      assembly = assembly
    )
  }

  # Generate README summary (if heart_pattern provided)
  if (nchar(heart_pattern) > 0) {
    .generateCardiacSummary(
      gene_table = gene_table,
      locus_cluster = locus_cluster,
      locus_name = locus_name,
      heart_pattern = heart_pattern,
      output_dir = output_dir
    )
  }

  # Generate founder SNP table (if SNP data available)
  if (!is.null(founder_snps) && nrow(founder_snps) > 0) {
    snps_in_locus <- .getFounderSnpsInRegion(
      target_chr = locus_chr,
      target_start = locus_start_bp,
      target_end = locus_end_bp,
      snp_data = founder_snps
    )

    if (nrow(snps_in_locus) > 0) {
      snp_table <- snps_in_locus[, .(
        SNP_ID = SNP,
        rsID,
        CHR,
        Position_mm39 = mm39,
        MajorAllele = Major,
        MinorAllele = Minor,
        MAF,
        Strains_with_Minor_Allele = SNPStrains
      )]

      snp_file <- file.path(output_dir, paste0("founder_snp_table_", locus_name, ".csv"))
      data.table::fwrite(snp_table, snp_file)
      message("Generated founder SNP table: ", snp_file)
    }
  }

  message("Locus packet complete: ", output_dir)
  return(invisible(output_dir))
}

#' Creates a LocusZoom-style visualization plot from project data
#'
#' Generates a multi-panel locus plot (Manhattan/scan track, optional signal
#' track, and gene track) for a single region in a locusPackRat project.
#'
#' @param region_id Character: ID of the region to plot. If \code{NULL},
#'   the first region in \code{regions.csv} is used.
#' @param project_dir Character: Path to the locusPackRat project directory
#'   containing the \code{.locusPackRat} folder.
#' @param assembly Plotgardener assembly object. If \code{NULL}, an assembly
#'   is constructed automatically based on the project \code{genome} in
#'   \code{config.json}.
#' @param scan_table Character: Name of the supplementary table (without
#'   \code{.csv}) containing scan/QTL data to use for the Manhattan track.
#' @param signal_table Character (optional): Name of the supplementary table
#'   (without \code{.csv}) to plot as a signal track below the Manhattan plot
#'   (e.g. founder effects, other scores). If \code{NULL}, the signal panel
#'   is skipped.
#' @param width Numeric: Width of the output PDF in inches.
#' @param height Numeric: Height of the output PDF in inches.
#' @param layout_ratios Named numeric vector giving the relative heights of
#'   the panels, typically \code{c(manhattan = ..., signal = ..., genes = ...)}.
#'   Values are rescaled internally to fill the available page height.
#' @param highlight_genes Character vector of gene symbols to highlight in
#'   the gene track (optional).
#' @param threshold Numeric: LOD significance threshold used to draw the
#'   horizontal line on the Manhattan plot (e.g. 5 for a 5-LOD threshold).
#' @param output_file Character: Output filename. Extension sets the device
#'   (\code{.png} for PNG, \code{.pdf} for PDF; default \code{"locus_zoom.pdf"}).
#'   Written to \code{.locusPackRat/output/} under \code{project_dir}.
#' @param text_scale Numeric (optional): Scaling factor for all text elements
#'   in the plot. If \code{NULL}, a default scaling based on the
#'   \code{width} and \code{height} is applied (6x4 inches is considered
#'   "normal" size with no scaling).
#' @param font_family Character (optional): Font family to use for all text
#'   elements in the plot. If \code{NULL}, the default system font is used.
#' @return Invisible file path to the saved plot.
#'
#' @details
#' ## Input File Formats
#'
#' ### scan_table (Manhattan plot)
#' The scan table CSV must contain:
#' \itemize{
#'   \item \code{chr}: Chromosome (character, e.g., "3" not 3)
#'   \item \code{pos}: Position in base pairs (or use \code{start} or \code{bp})
#'   \item \code{p}: P-value (or use \code{lod}, which auto-converts via \code{10^(-lod)})
#' }
#'
#' ### signal_table (optional signal track)
#' Two formats are supported:
#'
#' **Founder Effects (8-way crosses like CC/DO mice):**
#' \itemize{
#'   \item \code{chr}: Chromosome (character)
#'   \item \code{start}: Position in base pairs
#'   \item \code{A}, \code{B}, \code{C}, \code{D}, \code{E}, \code{F}, \code{G}, \code{H}:
#'     Effect values for each founder strain
#' }
#'
#' **Generic Signal:**
#' \itemize{
#'   \item \code{chr}: Chromosome (character)
#'   \item \code{start}: Position in base pairs
#'   \item \code{score}: Signal value (or use \code{value}, \code{lod},
#'     \code{intensity}, or \code{effect})
#' }
#'
#' For both formats, the \code{end} column is auto-calculated from \code{start}
#' positions if not provided.
#'
#' @importFrom data.table fread fwrite setDT setnames as.data.table setorderv copy
#' @importFrom jsonlite read_json
#' @importFrom grDevices pdf png dev.off
#' @importFrom tools file_ext
#' @importFrom plotgardener pageCreate plotManhattan plotGenes plotRanges plotText plotSignal annoYaxis annoHighlight
#' @importFrom utils capture.output
#'
#' @examples
#' \dontrun{
#' # Minimal example: Manhattan + genes
#' generateLocusZoomPlot(
#'   region_id   = "EF21_Iso_Peak_Zoom",
#'   project_dir = "EF21_Iso_Analysis",
#'   scan_table  = "full_iso_scan",
#'   output_file = "basictest.pdf"
#' )
#'
#' # Full example: with founder effects, custom layout, and highlights
#' generateLocusZoomPlot(
#'   region_id    = "EF21_Iso_Peak_Zoom",
#'   project_dir  = "EF21_Iso_Analysis",
#'   scan_table   = "full_iso_scan",
#'   signal_table = "full_iso_founders",
#'   width        = 10,
#'   height       = 6,
#'   layout_ratios = c(manhattan = 0.35, signal = 0.40, genes = 0.25),
#'   highlight_genes = c("Tspan5", "Stpg2"),
#'   threshold       = 3.1,
#'   text_scale      = 1.5,
#'   font_family     = "Helvetica",
#'   output_file     = "locus_zoom_region3_full.pdf"
#' )
#' }
#'
#' @export
generateLocusZoomPlot <- function(
    region_id = NULL,
    project_dir = ".",
    assembly = NULL,
    # Data Sources
    scan_table = NULL,
    signal_table = NULL,
    # Aesthetics / layout
    width = 6,
    height = 4,
    layout_ratios = c(manhattan = 0.4, signal = 0.3, genes = 0.3),
    highlight_genes = NULL,
    threshold = 5,
    output_file = "locus_zoom.pdf",
    text_scale = NULL,
    font_family = NULL
) {
    # 1. SETUP & CONFIG LOAD ----------------------------------------------
    packrat_dir <- file.path(project_dir, ".locusPackRat")
    if (!dir.exists(packrat_dir)) stop("Project not initialized.")
    
    config <- jsonlite::read_json(file.path(packrat_dir, "config.json"))
    
    # Load Region info
    regions <- data.table::fread(file.path(packrat_dir, "input", "regions.csv"))
    if (is.null(region_id)) region_id <- regions$region_id[1]
    
    # IMPORTANT: correct data.table filtering
    temp=region_id # otherwise it won't play nice.
    target_region <- regions[which(regions$region_id == temp), ]
    if (nrow(target_region) != 1L) {
        stop("Expected exactly 1 row for region_id = ", region_id,
             ", got ", nrow(target_region))
    }
    
    target_chr <- target_region$chr
    
    # 1a. TEXT SCALING SETUP ----------------------------------------------
    # Take 6x4 as "design size" for defaults
    if (is.null(text_scale)) {
        base_w <- 6
        base_h <- 4
        text_scale <- min(width / base_w, height / base_h)
    }
    
    base_sizes <- list(
        axis_tick    = 8,
        axis_label   = 10,
        panel_label  = 11,
        gene_label   = 8,
        genome_label = 8,
        title        = 12,
        signal_label = 10
    )
    text_sizes <- lapply(base_sizes, `*`, text_scale)
    
    # 1b. DETERMINE ASSEMBLY ----------------------------------------------
    if (is.null(assembly)) {
        if (config$genome == "mm39") {
            req_tx <- "TxDb.Mmusculus.UCSC.mm39.knownGene"
            req_org <- "org.Mm.eg.db"
            genome_build <- "mm39"
        } else if (config$genome == "hg38") {
            req_tx <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
            req_org <- "org.Hs.eg.db"
            genome_build <- "hg38"
        } else if (config$genome == "mm10") {
            req_tx <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
            req_org <- "org.Mm.eg.db"
            genome_build <- "mm10"
        } else if (config$genome == "hg19") {
            req_tx <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
            req_org <- "org.Hs.eg.db"
            genome_build <- "hg19"
        } else {
            stop("Unsupported genome in config: ", config$genome)
        }
        
        missing_pkgs <- character()
        if (!requireNamespace(req_tx, quietly = TRUE))  missing_pkgs <- c(missing_pkgs, req_tx)
        if (!requireNamespace(req_org, quietly = TRUE)) missing_pkgs <- c(missing_pkgs, req_org)
        if (length(missing_pkgs) > 0) {
            install_cmd <- paste0(
                "BiocManager::install(c('",
                paste(missing_pkgs, collapse = "', '"),
                "'))"
            )
            stop(
                "To plot genes for ", config$genome,
                ", you must install the following annotation package(s): ",
                paste(missing_pkgs, collapse = ", "),
                "\nRun this command: ", install_cmd
            )
        }
        
        assembly <- plotgardener::assembly(
            Genome = genome_build,
            TxDb   = req_tx,
            OrgDb  = req_org
        )
    }
    
    # Set file to save to the 'output' directory of the project_dir
    if (basename(output_file) != output_file) {
        message("Note: output_file should be a filename, not a path. Using basename: ",
                basename(output_file))
        output_file <- basename(output_file)
    }
    output_path <- file.path(project_dir, ".locusPackRat", "output", output_file)
    ext <- tolower(tools::file_ext(output_file))
    if (ext == "png") {
        grDevices::png(output_path, width = width, height = height,
                       units = "in", res = 150)
    } else {
        grDevices::pdf(output_path, width = width, height = height)
    }
    .quiet_pg(plotgardener::pageCreate(
        width = width,
        height = height,
        default.units = "inches",
        showGuides = FALSE
    ))
    
    # Margins and gaps as FRACTIONS of the page
    margin_left_frac   <- 0.12
    margin_right_frac  <- 0.06
    margin_top_frac    <- 0.08
    margin_bottom_frac <- 0.10
    gap_frac           <- 0.04   # between panels
    genome_label_frac  <- 0.06   # space reserved under gene panel
    
    margin_left   <- width  * margin_left_frac
    margin_right  <- width  * margin_right_frac
    margin_top    <- height * margin_top_frac
    margin_bottom <- height * margin_bottom_frac
    panel_gap     <- height * gap_frac
    genome_extra  <- height * genome_label_frac
    
    plot_width <- width - (margin_left + margin_right)
    plot_height <- height - (margin_top + margin_bottom + genome_extra + 2 * panel_gap)
    
    if (plot_width <= 0 || plot_height <= 0) {
        grDevices::dev.off()
        stop("Non-positive plotting area; reduce margins or increase width/height.")
    }
    
    # 1d. pgParams for the region -----------------------------------------
    # Ensure target_chr is character (plotgardener handles prefix internally)
    target_chr <- as.character(target_chr)

    pg_params <- plotgardener::pgParams(
        assembly   = assembly,
        chrom      = target_chr,
        chromstart = as.numeric(target_region$start),
        chromend   = as.numeric(target_region$end)
    )
    
    # 2. DATA LOADING ------------------------------------------------------
    if (is.null(scan_table))
        stop("Must provide a scan_table name")
    
    scan_path <- file.path(packrat_dir, "supplementary", paste0(scan_table, ".csv"))
    if (!file.exists(scan_path))
        stop("Scan table not found: ", scan_path)
    
    scan_dt <- data.table::fread(scan_path)
    scan_dt <- .standardize_chrom(scan_dt)

    # Signal table is truly optional
    signal_dt <- NULL
    if (!is.null(signal_table)) {
        signal_path <- file.path(packrat_dir, "supplementary", paste0(signal_table, ".csv"))
        if (!file.exists(signal_path)) {
            grDevices::dev.off()
            stop("Signal table not found: ", signal_path)
        }
        signal_dt <- data.table::fread(signal_path)
        signal_dt <- .standardize_chrom(signal_dt)
    }
    
    # Genes / highlights
    gene_list_vector <- NULL
    if (config$mode == "region") {
        if ("genes" %in% names(target_region) && !is.na(target_region$genes)) {
            raw_genes <- strsplit(as.character(target_region$genes), ",")[[1]]
            gene_list_vector <- trimws(raw_genes)
        }
    } else {
        gene_file <- file.path(packrat_dir, "input", "genes.csv")
        if (file.exists(gene_file)) {
            gene_dt <- data.table::fread(gene_file)
            gene_list_vector <- gene_dt$gene_symbol
        }
    }
    
    hl_df <- NULL
    if (!is.null(highlight_genes)) {
        hl_df <- data.frame(gene = highlight_genes, color = "#a03b60")
    }
    
    # 3. PANEL HEIGHTS (normalize ratios over drawn panels) ----------------
    # We always draw manhattan + genes; signals only if we have data rows.
    has_signal_data <- !is.null(signal_dt) && nrow(signal_dt) > 0
    
    active_panels <- c("manhattan", if (has_signal_data) "signal", "genes")
    ratios <- layout_ratios[active_panels]
    if (any(is.na(ratios))) {
        grDevices::dev.off()
        stop("layout_ratios must have named entries for: ",
             paste(active_panels, collapse = ", "))
    }
    ratios <- ratios / sum(ratios)
    
    h_manhattan <- plot_height * ratios["manhattan"]
    h_genes     <- plot_height * ratios["genes"]
    h_signals   <- if (has_signal_data) plot_height * ratios["signal"] else 0
    
    # 4. RENDER MODULES ----------------------------------------------------
    current_y <- margin_top
    
    # A. Manhattan / scan
    current_y <- .render_manhattan(
        scan_data   = scan_dt,
        params      = pg_params,
        x           = margin_left,
        y           = current_y,
        w           = plot_width,
        h           = h_manhattan,
        threshold   = threshold,
        text_sizes  = text_sizes,
        font_family = font_family
    )
    current_y <- current_y + panel_gap
    
    # B. Signals (if available)
    if (has_signal_data && h_signals > 0) {
        new_y <- .render_signals(
            signal_dt  = signal_dt,
            params     = pg_params,
            x          = margin_left,
            y          = current_y,
            w          = plot_width,
            h          = h_signals,
            text_sizes = text_sizes,
            font_family = font_family
        )
        # Only add gap if something was actually drawn
        if (new_y > current_y) {
            current_y <- new_y + panel_gap
        }
    }
    
    # C. Genes + genome label
    current_y <- .render_genes(
        gene_list    = gene_list_vector,
        highlights   = hl_df,
        params       = pg_params,
        x            = margin_left,
        y            = current_y,
        w            = plot_width,
        h            = h_genes,
        genome_extra = genome_extra,
        text_sizes   = text_sizes,
       font_family  = font_family
    )
    
    grDevices::dev.off()
    message("Plot saved to ", output_file)
    message("Citation: Kramer NE, Davis ES, Wenger CD, Deoudes EM, Parker SM, Love MI, Phanstiel DH. Plotgardener: cultivating precise multi-panel figures in R. Bioinformatics, 2022.")
    invisible(output_path)
}

.quiet_pg <- function(expr) {
  capture.output(
    res <- suppressMessages(expr),
    file = NULL
  )
  res
}

#' Standardize chromosome column to character
#' @noRd
.standardize_chrom <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return(dt)

  # Ensure chr column is character (plotgardener requires this)
  if ("chr" %in% names(dt)) {
    dt[, chr := as.character(chr)]
  }

  # Create chrom column for internal filtering (same value, no prefix needed)
  if (!"chrom" %in% names(dt) && "chr" %in% names(dt)) {
    dt[, chrom := chr]
  }

  return(dt)
}

#' Internal module to render Genes
#' @noRd
.render_genes <- function(gene_list, highlights, params, x, y, w, h) {
  
  # Determine labeling priority:
  # 1. Highlighted genes (if any)
  # 2. Genes listed in the region (if available)
  # 3. Plotgardener defaults
  
  gene_order <- NULL
  
  if (!is.null(highlights)) {
    # Start with highlights
    gene_order <- highlights$gene
    
    # Append the rest of the known genes in the region (if we have them)
    if (!is.null(gene_list)) {
      remaining <- setdiff(gene_list, highlights$gene)
      gene_order <- c(gene_order, remaining)
    }
  } else if (!is.null(gene_list)) {
    # No highlights, just use the region list order
    gene_order <- gene_list
  }
  
  # Add "chr" prefix if needed
  params$chrom <- if (grepl("^chr", params$chrom)) {
    params$chrom
  } else {
    paste0("chr", params$chrom)
  }
  # Plot Genes
  # plotgardener pulls coordinates from the 'assembly' in params,
  # so we don't need to pass a dataframe of coordinates here.
  plt <- plotgardener::plotGenes(
    params = params,
    x = x, y = y, width = w, height = h,
    just = c("left", "top"),
    geneOrder = gene_order,
    geneHighlights = highlights,
    geneBackground = "grey",
    fontsize = 8,
    default.units = "inches"
  )
  
  # Add Genome Label below genes
  plotgardener::annoGenomeLabel(
    plot = plt, x = x, y = y + h + 0.2,
    scale = "Mb", just = c("left", "top"),
    default.units = "inches"
  )
  
  return(y + h + 0.5)
}

.render_manhattan <- function(
    scan_data,
    params,
    x, y, w, h,
    threshold,
    text_sizes,
    font_family = NULL
) {
    # Empty?
    if (is.null(scan_data) || nrow(scan_data) == 0) return(y)

    # Standardize pos column if needed
    if (!"pos" %in% names(scan_data)) {
        if ("start" %in% names(scan_data))      scan_data[, pos := start]
        else if ("bp" %in% names(scan_data))    scan_data[, pos := bp]
    }
    if ("lod" %in% names(scan_data) && !"p" %in% names(scan_data)) {
        scan_data[, p := 10^(-lod)]
    }
    
    # Validate required
    required <- c("chrom", "pos", "p")
    missing  <- setdiff(required, names(scan_data))
    if (length(missing)) {
        warning("scan_data is missing required columns: ",
                paste(missing, collapse = ", "),
                "; skipping Manhattan panel.")
        return(y)
    }
    
    # Filter for chromosome
    if (isTRUE(!is.na(params$chrom))) {
        scan_data <- scan_data[chrom == params$chrom]
    }
    if (nrow(scan_data) == 0) return(y)
    
    max_y <- max(c(-log10(scan_data$p), threshold, 5), na.rm = TRUE) + 1
    ylim <- c(0, max_y)

    # Convert to data.frame for plotgardener compatibility
    scan_data <- as.data.frame(scan_data)

    plt <- .quiet_pg(plotgardener::plotManhattan(
        data    = scan_data,
        params  = params,
        range   = ylim,
        trans   = "-log10",
        sigVal  = 10^(-threshold),
        x       = x,
        y       = y,
        width   = w,
        height  = h,
        just    = c("left", "top"),
        fill    = "#a1c6d1",
        sigCol  = "#53add2",
        sigLine = TRUE,
        baseline = TRUE,
        default.units = "inches"
    ))
    
    # Y axis with scaled tick text
    if (requireNamespace("grid", quietly = TRUE)) {
        gp_args <- list(fontsize = text_sizes$axis_tick)
        if (!is.null(font_family)) gp_args$fontfamily <- font_family

        .quiet_pg(do.call(
            plotgardener::annoYaxis,
            c(
                list(
                    plot     = plt,
                    at       = pretty(ylim),
                    axisLine = TRUE
                ),
                gp_args
            )
        ))
    } else {
        .quiet_pg(plotgardener::annoYaxis(
            plot     = plt,
            at       = pretty(ylim),
            axisLine = TRUE
        ))
    }
    
    # Axis label ("LOD"), size-aware and offset as fraction of width
    axis_offset <- w * 0.06
    .quiet_pg(plotgardener::plotText(
        label  = "LOD",
        x      = x - axis_offset,
        y      = y + (h / 2),
        rot    = 90,
        fontsize = text_sizes$axis_label,
        just   = "center",
        default.units = "inches",
        fontfamily = font_family
    ))
    
    y + h
}

.render_signals <- function(
    signal_dt,
    params,
    x, y, w, h,
    text_sizes,
    font_family = NULL
) {
    if (is.null(signal_dt) || nrow(signal_dt) == 0) return(y)

    target_chrom <- params$chrom

    if (!"start" %in% colnames(signal_dt)) {
        warning("signal_dt must include 'start' column; skipping signal panel.")
        return(y)
    }
    
    # Filter by chromosome
    if (isTRUE(!is.na(target_chrom))) {
        signal_dt <- signal_dt[chrom == target_chrom]
    }
    if (nrow(signal_dt) == 0) return(y)
    
    data.table::setorderv(signal_dt, "start")
    signal_dt[, end := data.table::shift(start, type = "lead", fill = start[.N] + 2) - 1]
    
    founder_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
    is_founder_data <- all(founder_cols %in% colnames(signal_dt))
    
    if (is_founder_data) {
        strain_colors <- c(
            "A" = "#1B9E77", "B" = "#D95F02", "C" = "#7570B3", "D" = "#E7298A",
            "E" = "#66A61E", "F" = "#E6AB02", "G" = "#A6761D", "H" = "#666666"
        )
        
        all_vals <- as.vector(as.matrix(signal_dt[, ..founder_cols]))
        max_eff <- round(max(abs(all_vals), na.rm = TRUE) * 1.1, digits = 0)
        if (!is.finite(max_eff) || max_eff == 0) {
            return(y)  # nothing meaningful to plot
        }
        sig_range <- c(-max_eff, max_eff)
        
        first_plot <- NULL
        for (strain_let in founder_cols) {
            sub_dt <- signal_dt[, .(chrom, start, end, score = get(strain_let))]
            # Convert to data.frame for plotgardener compatibility
            sub_dt <- as.data.frame(sub_dt)

            plt <- .quiet_pg(
                plotgardener::plotSignal(
                    data    = sub_dt,
                    params  = params,
                    range   = sig_range,
                    negData = TRUE,
                    linecolor = strain_colors[strain_let],
                    fill    = NA,
                    x       = x,
                    y       = y,
                    width   = w,
                    height  = h,
                    just    = c("left", "top"),
                    baseline = TRUE,
                    baseline.color = "grey",
                    default.units  = "inches"
                )
            )

            if (is.null(first_plot)) first_plot <- plt
        }
        
        if (!is.null(first_plot)) {
            # Y axis
            if (requireNamespace("grid", quietly = TRUE)) {
                gp_args <- list(fontsize = text_sizes$axis_tick)
                if (!is.null(font_family)) gp_args$fontfamily <- font_family

                .quiet_pg(do.call(
                    plotgardener::annoYaxis,
                    c(
                        list(
                            plot     = first_plot,
                            at       = c(-max_eff, 0, max_eff),
                            axisLine = TRUE
                        ),
                        gp_args
                    )
                ))
            } else {
                .quiet_pg(plotgardener::annoYaxis(
                    plot     = first_plot,
                    at       = c(-max_eff, 0, max_eff),
                    axisLine = TRUE
                ))
            }
        }
        
        # Side label
        axis_offset <- w * 0.08
        .quiet_pg(plotgardener::plotText(
            label   = "Founder Effects",
            x       = x - axis_offset,
            y       = y + (h / 2),
            rot     = 90,
            fontsize = text_sizes$signal_label,
            just    = "center",
            default.units = "inches",
            fontfamily = font_family
        ))
        
    } else {
        # Standard single signal
        if (!"score" %in% names(signal_dt)) {
            val_col <- intersect(names(signal_dt), c("value", "lod", "intensity", "effect"))[1]
            if (!is.na(val_col)) {
                signal_dt[, score := get(val_col)]
            } else {
                warning("Could not find a value column for signal_dt; skipping signal panel.")
                return(y)
            }
        }
        
        max_val <- max(abs(signal_dt$score), na.rm = TRUE)
        if (!is.finite(max_val) || max_val == 0) return(y)
        sig_range <- c(-max_val * 1.1, max_val * 1.1)

        # Convert to data.frame for plotgardener compatibility
        signal_dt <- as.data.frame(signal_dt)

        plt <- .quiet_pg(
        plotgardener::plotSignal(
            data    = signal_dt,
            params  = params,
            range   = sig_range,
            negData = TRUE,
            linecolor = "#377eb8",
            fill    = NA,
            x       = x,
            y       = y,
            width   = w,
            height  = h,
            just    = c("left", "top"),
            baseline = TRUE,
            baseline.color = "grey",
            default.units  = "inches"
        )
        )
        
        # Y axis
        if (requireNamespace("grid", quietly = TRUE)) {
            gp_args <- list(fontsize = text_sizes$axis_tick)
            if (!is.null(font_family)) gp_args$fontfamily <- font_family

            .quiet_pg(do.call(
                plotgardener::annoYaxis,
                c(
                    list(
                        plot     = plt,
                        at       = pretty(sig_range, n = 3),
                        axisLine = TRUE
                    ),
                    gp_args
                )
            ))
        } else {
            .quiet_pg(plotgardener::annoYaxis(
                plot     = plt,
                at       = pretty(sig_range, n = 3),
                axisLine = TRUE
            ))
        }

        # Panel label on the side (generic)
        axis_offset <- w * 0.08
        .quiet_pg(plotgardener::plotText(
            label   = "Signal",
            x       = x - axis_offset,
            y       = y + (h / 2),
            rot     = 90,
            fontsize = text_sizes$signal_label,
            just    = "center",
            default.units = "inches",
            fontfamily = font_family
        ))
    }
    
    y + h
}

.render_genes <- function(
    gene_list,
    highlights,
    params,
    x, y, w, h,
    genome_extra,
    text_sizes,
    font_family = NULL
) {
    # Gene labeling priority:
    # 1. highlights
    # 2. region gene_list
    gene_order <- NULL
    if (!is.null(highlights)) {
        gene_order <- highlights$gene
        if (!is.null(gene_list)) {
            remaining <- setdiff(gene_list, highlights$gene)
            gene_order <- c(gene_order, remaining)
        }
    } else if (!is.null(gene_list)) {
        gene_order <- gene_list
    }
    
    # plotGenes wants chrom with "chr" prefix if annotations are like that
    params$chrom <- if (grepl("^chr", params$chrom)) params$chrom else paste0("chr", params$chrom)
    
    # Gene track
    plt <- .quiet_pg(plotgardener::plotGenes(
        params        = params,
        x             = x,
        y             = y,
        width         = w,
        height        = h,
        just          = c("left", "top"),
        geneOrder     = gene_order,
        geneHighlights = highlights,
        geneBackground = "grey",
        fontsize      = text_sizes$gene_label,
        default.units = "inches"
        # font_family is not explicitly supported here; fontcolor
        # can be customized separately if you want.
    ))
    
    # Genome label underneath, scaled
    .quiet_pg(
    suppressWarnings(
        plotgardener::annoGenomeLabel(
        plot        = plt,
        x           = x,
        y           = y + h + (genome_extra * 0.3),
        scale       = "Mb",
        fontsize    = text_sizes$genome_label,
        fontcolor   = "black",
        just        = c("left", "top"),
        default.units = "inches"
        )
    )
    )

    
    y + h + genome_extra
}

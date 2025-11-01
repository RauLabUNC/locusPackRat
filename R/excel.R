#' Create multi-sheet Excel workbook for gene curation
#' 
#' @description
#' Generate comprehensive Excel workbooks for manual gene curation.
#' Each workbook contains multiple sheets with different data views.
#' 
#' @param gene_data Primary gene data table
#' @param output_file Path for output Excel file
#' @param sheet_config List defining which sheets to create
#' @param styles List of style specifications for formatting
#' @param freeze_panes Freeze first row and/or column (default: TRUE)
#' @param conditional_format Apply conditional formatting rules (default: TRUE)
#' @param hyperlinks Add hyperlinks to external resources (default: TRUE)
#' 
#' @return Invisible NULL (creates file as side effect)
#' 
#' @examples
#' # Basic workbook creation
#' genes <- data.frame(
#'   gene = c("Abcb10", "Acsl5"),
#'   biotype = c("protein_coding", "protein_coding"),
#'   expression = c(150, 90)
#' )
#' create_gene_workbook(genes, "output.xlsx")
#' 
#' # Custom sheet configuration
#' config <- list(
#'   overview = TRUE,
#'   detailed = TRUE,
#'   variants = TRUE,
#'   phenotypes = TRUE,
#'   custom_sheets = list(
#'     high_priority = function(df) df[df$expression > 100, ]
#'   )
#' )
#' create_gene_workbook(genes, "output.xlsx", sheet_config = config)
#' 
#' @export
#' @importFrom openxlsx createWorkbook addWorksheet writeData createStyle addStyle freezePane
create_gene_workbook <- function(gene_data,
                                output_file,
                                sheet_config = NULL,
                                styles = NULL,
                                freeze_panes = TRUE,
                                conditional_format = TRUE,
                                hyperlinks = TRUE) {

  # Helper function to safely subset columns
  subset_columns <- function(data, cols) {
    cols <- intersect(cols, names(data))
    if (length(cols) == 0) return(NULL)

    if (inherits(data, "data.table")) {
      return(data[, ..cols])
    } else {
      return(data[, cols, drop = FALSE])
    }
  }

  # Default sheet configuration
  if (is.null(sheet_config)) {
    sheet_config <- list(
      overview = TRUE,      # High-level summary
      detailed = TRUE,      # All columns
      expression = TRUE,    # Expression focused
      variants = TRUE,      # Variant focused
      phenotypes = TRUE,    # Phenotype focused
      prioritized = TRUE    # Top candidates
    )
  }
  
  # Create workbook
  wb <- openxlsx::createWorkbook()
  
  # Default styles
  if (is.null(styles)) {
    styles <- list(
      header = openxlsx::createStyle(
        fontSize = 11,
        fontColour = "#FFFFFF",
        fgFill = "#4472C4",
        halign = "center",
        valign = "center",
        textDecoration = "bold",
        border = "TopBottomLeftRight",
        borderColour = "#000000"
      ),
      
      highlight = openxlsx::createStyle(
        fgFill = "#FFE699",
        border = "TopBottomLeftRight"
      ),
      
      number = openxlsx::createStyle(
        numFmt = "0.00",
        halign = "right"
      ),
      
      hyperlink = openxlsx::createStyle(
        fontColour = "#0563C1",
        textDecoration = "underline"
      )
    )
  }
  
  # Sheet 1: Overview
  if (sheet_config$overview) {
    sheet_name <- "Overview"
    openxlsx::addWorksheet(wb, sheet_name)
    
    # Select key columns - try different gene column names
    gene_col <- intersect(c("gene_id", "gene", "symbol", "gene_symbol"), names(gene_data))[1]
    if (is.na(gene_col)) gene_col <- names(gene_data)[1]  # fallback to first column

    overview_cols <- c(gene_col, "biotype", "chromosome", "chr", "start", "end",
                      "human_ortholog", "human_symbol", "expression", "has_eqtl",
                      "n_variants", "n_phenotypes")

    overview_data <- subset_columns(gene_data, overview_cols)
    if (is.null(overview_data)) overview_data <- gene_data
    
    openxlsx::writeData(wb, sheet_name, overview_data, 
                       startRow = 1, startCol = 1, 
                       headerStyle = styles$header)
    
    if (freeze_panes) {
      openxlsx::freezePane(wb, sheet_name, firstRow = TRUE, firstCol = TRUE)
    }
    
    # Conditional formatting for expression
    if (conditional_format && "expression" %in% names(overview_data)) {
      expr_col <- which(names(overview_data) == "expression")
      openxlsx::conditionalFormatting(
        wb, sheet_name,
        cols = expr_col,
        rows = 2:(nrow(overview_data) + 1),
        type = "databar",
        style = c("#63BE7B", "#FFEB84", "#F8696B"),
        rule = NULL
      )
    }
  }
  
  # Sheet 2: Detailed view
  if (sheet_config$detailed) {
    sheet_name <- "Detailed"
    openxlsx::addWorksheet(wb, sheet_name)
    
    openxlsx::writeData(wb, sheet_name, gene_data,
                       startRow = 1, startCol = 1,
                       headerStyle = styles$header)
    
    if (freeze_panes) {
      openxlsx::freezePane(wb, sheet_name, firstRow = TRUE, firstCol = TRUE)
    }
    
    # Auto-fit column widths
    openxlsx::setColWidths(wb, sheet_name, 
                          cols = 1:ncol(gene_data), 
                          widths = "auto")
  }
  
  # Sheet 3: Expression-focused
  if (sheet_config$expression && any(grepl("express|CPM|TPM|FPKM", names(gene_data)))) {
    sheet_name <- "Expression"
    openxlsx::addWorksheet(wb, sheet_name)
    
    gene_col <- intersect(c("gene_id", "gene", "symbol", "gene_symbol"), names(gene_data))[1]
    if (is.na(gene_col)) gene_col <- names(gene_data)[1]

    expr_cols <- c(gene_col, names(gene_data)[grepl("express|CPM|TPM|FPKM|eqtl",
                                                    names(gene_data),
                                                    ignore.case = TRUE)])
    expr_data <- subset_columns(gene_data, expr_cols)
    if (is.null(expr_data)) expr_data <- gene_data
    
    # Sort by expression if available
    if (any(grepl("CPM|TPM|FPKM", names(expr_data)))) {
      expr_col <- names(expr_data)[grepl("CPM|TPM|FPKM", names(expr_data))][1]
      expr_data <- expr_data[order(expr_data[[expr_col]], decreasing = TRUE), ]
    }
    
    openxlsx::writeData(wb, sheet_name, expr_data,
                       startRow = 1, startCol = 1,
                       headerStyle = styles$header)
    
    if (freeze_panes) {
      openxlsx::freezePane(wb, sheet_name, firstRow = TRUE, firstCol = TRUE)
    }
    
    # Apply number format to expression columns
    for (col in which(sapply(expr_data, is.numeric))) {
      openxlsx::addStyle(wb, sheet_name, 
                        style = styles$number,
                        rows = 2:(nrow(expr_data) + 1),
                        cols = col,
                        gridExpand = TRUE)
    }
  }
  
  # Sheet 4: Variants
  if (sheet_config$variants && any(grepl("variant|mutation|SNP", names(gene_data)))) {
    sheet_name <- "Variants"
    openxlsx::addWorksheet(wb, sheet_name)
    
    gene_col <- intersect(c("gene_id", "gene", "symbol", "gene_symbol"), names(gene_data))[1]
    if (is.na(gene_col)) gene_col <- names(gene_data)[1]

    var_cols <- c(gene_col, names(gene_data)[grepl("variant|mutation|SNP|missense|nonsense",
                                                   names(gene_data),
                                                   ignore.case = TRUE)])
    var_data <- subset_columns(gene_data, var_cols)
    if (is.null(var_data)) var_data <- gene_data
    
    # Filter to genes with variants if column exists
    if ("n_variants" %in% names(var_data)) {
      var_data <- var_data[var_data$n_variants > 0, ]
    }
    
    openxlsx::writeData(wb, sheet_name, var_data,
                       startRow = 1, startCol = 1,
                       headerStyle = styles$header)
    
    if (freeze_panes) {
      openxlsx::freezePane(wb, sheet_name, firstRow = TRUE, firstCol = TRUE)
    }
  }
  
  # Sheet 5: Phenotypes
  if (sheet_config$phenotypes && any(grepl("phenotype|disease|trait", names(gene_data)))) {
    sheet_name <- "Phenotypes"
    openxlsx::addWorksheet(wb, sheet_name)
    
    gene_col <- intersect(c("gene_id", "gene", "symbol", "gene_symbol"), names(gene_data))[1]
    if (is.na(gene_col)) gene_col <- names(gene_data)[1]

    pheno_cols <- c(gene_col, names(gene_data)[grepl("phenotype|disease|trait|syndrome",
                                                     names(gene_data),
                                                     ignore.case = TRUE)])
    pheno_data <- subset_columns(gene_data, pheno_cols)
    if (is.null(pheno_data)) pheno_data <- gene_data
    
    # Filter to genes with phenotypes if column exists
    if ("n_phenotypes" %in% names(pheno_data)) {
      pheno_data <- pheno_data[pheno_data$n_phenotypes > 0, ]
    }
    
    openxlsx::writeData(wb, sheet_name, pheno_data,
                       startRow = 1, startCol = 1,
                       headerStyle = styles$header)
    
    if (freeze_panes) {
      openxlsx::freezePane(wb, sheet_name, firstRow = TRUE, firstCol = TRUE)
    }
  }
  
  # Sheet 6: Prioritized candidates
  if (sheet_config$prioritized) {
    sheet_name <- "Top_Candidates"
    openxlsx::addWorksheet(wb, sheet_name)
    
    # Apply prioritization logic
    priority_data <- prioritize_genes(gene_data)
    
    openxlsx::writeData(wb, sheet_name, priority_data,
                       startRow = 1, startCol = 1,
                       headerStyle = styles$header)
    
    if (freeze_panes) {
      openxlsx::freezePane(wb, sheet_name, firstRow = TRUE, firstCol = TRUE)
    }
    
    # Highlight top genes
    if (nrow(priority_data) > 0) {
      openxlsx::addStyle(wb, sheet_name,
                        style = styles$highlight,
                        rows = 2:min(6, nrow(priority_data) + 1),
                        cols = 1:ncol(priority_data),
                        gridExpand = TRUE)
    }
  }
  
  # Add hyperlinks if requested
  if (hyperlinks) {
    add_hyperlinks_to_workbook(wb, gene_data)
  }
  
  # Add custom sheets if provided
  if ("custom_sheets" %in% names(sheet_config)) {
    for (sheet_name in names(sheet_config$custom_sheets)) {
      sheet_function <- sheet_config$custom_sheets[[sheet_name]]
      custom_data <- sheet_function(gene_data)
      
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, custom_data,
                         startRow = 1, startCol = 1,
                         headerStyle = styles$header)
      
      if (freeze_panes) {
        openxlsx::freezePane(wb, sheet_name, firstRow = TRUE, firstCol = TRUE)
      }
    }
  }
  
  # Save workbook
  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
  
  invisible(NULL)
}

#' Prioritize genes based on multiple criteria
#' 
#' @description
#' Internal function to prioritize genes for the Top Candidates sheet.
#' 
#' @param gene_data Gene data table
#' @return Prioritized subset of gene data
#' 
#' @keywords internal
prioritize_genes <- function(gene_data) {
  
  # Add priority score
  gene_data$priority_score <- 0
  
  # Scoring criteria
  if ("biotype" %in% names(gene_data)) {
    gene_data$priority_score <- gene_data$priority_score + 
      ifelse(gene_data$biotype == "protein_coding", 2, 0)
  }
  
  if ("human_ortholog" %in% names(gene_data)) {
    gene_data$priority_score <- gene_data$priority_score + 
      ifelse(!is.na(gene_data$human_ortholog), 2, 0)
  }
  
  if ("expression" %in% names(gene_data) || "CPM" %in% names(gene_data)) {
    expr_col <- if ("expression" %in% names(gene_data)) "expression" else "CPM"
    # Score based on expression quartiles
    expr_quartiles <- quantile(gene_data[[expr_col]], 
                              probs = c(0.25, 0.5, 0.75), 
                              na.rm = TRUE)
    gene_data$priority_score <- gene_data$priority_score + 
      ifelse(gene_data[[expr_col]] > expr_quartiles[3], 3,
             ifelse(gene_data[[expr_col]] > expr_quartiles[2], 2,
                    ifelse(gene_data[[expr_col]] > expr_quartiles[1], 1, 0)))
  }
  
  if ("has_eqtl" %in% names(gene_data)) {
    gene_data$priority_score <- gene_data$priority_score + 
      ifelse(gene_data$has_eqtl == TRUE, 2, 0)
  }
  
  if ("n_variants" %in% names(gene_data)) {
    gene_data$priority_score <- gene_data$priority_score + 
      ifelse(gene_data$n_variants > 0, 1, 0) +
      ifelse(gene_data$n_variants > 2, 1, 0)
  }
  
  if ("n_phenotypes" %in% names(gene_data)) {
    gene_data$priority_score <- gene_data$priority_score + 
      ifelse(gene_data$n_phenotypes > 0, 1, 0) +
      ifelse(gene_data$n_phenotypes > 2, 1, 0)
  }
  
  # Sort by priority score
  gene_data <- gene_data[order(gene_data$priority_score, decreasing = TRUE), ]
  
  # Return top candidates (max 50)
  return(head(gene_data, 50))
}

#' Add hyperlinks to gene databases
#' 
#' @description
#' Internal function to add hyperlinks to external resources.
#' 
#' @param wb Workbook object
#' @param gene_data Gene data table
#' 
#' @keywords internal
add_hyperlinks_to_workbook <- function(wb, gene_data) {
  
  sheets <- openxlsx::sheets(wb)
  
  for (sheet in sheets) {
    sheet_data <- openxlsx::readWorkbook(wb, sheet = sheet)
    
    # Add hyperlinks for gene symbols
    if ("gene" %in% names(sheet_data)) {
      gene_col <- which(names(sheet_data) == "gene")
      
      # MGI links for mouse genes
      if (any(grepl("^[A-Z][a-z]", sheet_data$gene))) {
        for (i in 1:nrow(sheet_data)) {
          gene <- sheet_data$gene[i]
          url <- paste0("http://www.informatics.jax.org/searchtool/Search.do?query=", gene)
          openxlsx::writeFormula(wb, sheet,
                                x = paste0('HYPERLINK("', url, '","', gene, '")'),
                                startRow = i + 1, startCol = gene_col)
        }
      }
    }
    
    # Add hyperlinks for human orthologs
    if ("human_ortholog" %in% names(sheet_data)) {
      human_col <- which(names(sheet_data) == "human_ortholog")
      
      for (i in 1:nrow(sheet_data)) {
        if (!is.na(sheet_data$human_ortholog[i])) {
          gene <- sheet_data$human_ortholog[i]
          url <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene)
          openxlsx::writeFormula(wb, sheet,
                                x = paste0('HYPERLINK("', url, '","', gene, '")'),
                                startRow = i + 1, startCol = human_col)
        }
      }
    }
  }
}
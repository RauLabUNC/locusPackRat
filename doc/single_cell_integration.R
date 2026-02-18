## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = requireNamespace("Seurat", quietly = TRUE) &&
         requireNamespace("ggplot2", quietly = TRUE) &&
         dir.exists("my_qtl_project/.locusPackRat")
)

## ----libraries----------------------------------------------------------------
# # Core packages
# library(locusPackRat)
# library(data.table)
# 
# # scRNA-seq analysis
# library(Seurat)     # install from CRAN: install.packages("Seurat")
# library(ggplot2)

## ----extract_candidates-------------------------------------------------------
# # Load the master gene sheet from a completed locusPackRat project
# # This assumes you have already run initPackRat(), addRatTable(), etc.
# master_genes <- fread("my_qtl_project/.locusPackRat/input/genes.csv")
# 
# # Define candidate gene sets based on your prioritization criteria
# # These are example genes from the CC cardiac hypertrophy study
# cardiac_candidates <- c("Cisd2", "Pdlim5", "Manba", "Fhod3", "Myh7")
# 
# # You can also define broader gene sets for comparison
# all_locus_genes <- master_genes$gene_symbol

## ----load_tabula_muris--------------------------------------------------------
# # Download Tabula Muris heart data (10X Genomics)
# # Data available from: https://tabula-muris.ds.czbiohub.org/
# # Or via the TabulaMurisData Bioconductor package
# 
# # Using pre-downloaded .h5 or .rds file:
# heart_data <- Read10X(data.dir = "path/to/tabula_muris/Heart-10X_P7_4/")
# 
# # Create Seurat object
# heart_seurat <- CreateSeuratObject(
#   counts = heart_data,
#   project = "TabulaMuris_Heart",
#   min.cells = 3,
#   min.features = 200
# )
# 
# # Standard preprocessing
# heart_seurat <- NormalizeData(heart_seurat)
# heart_seurat <- FindVariableFeatures(heart_seurat, nfeatures = 2000)
# heart_seurat <- ScaleData(heart_seurat)
# heart_seurat <- RunPCA(heart_seurat)
# heart_seurat <- FindNeighbors(heart_seurat, dims = 1:20)
# heart_seurat <- FindClusters(heart_seurat, resolution = 0.5)
# heart_seurat <- RunUMAP(heart_seurat, dims = 1:20)

## ----load_preprocessed--------------------------------------------------------
# # Load a pre-processed Seurat object with cell type annotations
# heart_seurat <- readRDS("path/to/heart_annotated_seurat.rds")
# 
# # Check available cell type annotations
# table(heart_seurat$cell_type)

## ----module_scores------------------------------------------------------------
# # Define gene lists for module scoring
# # Convert mouse gene symbols to match the Seurat object's feature names
# candidate_list <- list(
#   cardiac_candidates = cardiac_candidates
# )
# 
# # Calculate module scores
# heart_seurat <- AddModuleScore(
#   object = heart_seurat,
#   features = candidate_list,
#   name = "candidate_score"
# )
# 
# # The score is added as a metadata column: "candidate_score1"
# # Higher scores indicate stronger collective expression of the candidate genes
# summary(heart_seurat$candidate_score1)

## ----plot_umap----------------------------------------------------------------
# # Plot module score on UMAP
# FeaturePlot(
#   heart_seurat,
#   features = "candidate_score1",
#   cols = c("lightgrey", "blue", "red"),
#   pt.size = 0.5
# ) +
#   ggtitle("Candidate Gene Module Score") +
#   theme(plot.title = element_text(hjust = 0.5))

## ----plot_violin--------------------------------------------------------------
# # Violin plot of module score by cell type
# VlnPlot(
#   heart_seurat,
#   features = "candidate_score1",
#   group.by = "cell_type",
#   pt.size = 0
# ) +
#   ggtitle("Candidate Gene Enrichment by Cell Type") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ----plot_individual_genes----------------------------------------------------
# # Dot plot of individual candidate gene expression by cell type
# DotPlot(
#   heart_seurat,
#   features = cardiac_candidates,
#   group.by = "cell_type"
# ) +
#   RotatedAxis() +
#   ggtitle("Individual Candidate Gene Expression by Cell Type")

## ----quantify_specificity-----------------------------------------------------
# # Extract expression matrix for candidate genes
# expr_matrix <- GetAssayData(heart_seurat, slot = "data")
# candidates_in_data <- intersect(cardiac_candidates, rownames(expr_matrix))
# 
# if (length(candidates_in_data) > 0) {
#   # Calculate mean expression per cell type per gene
#   cell_type_expr <- data.table(
#     cell_type = heart_seurat$cell_type,
#     as.data.table(t(as.matrix(expr_matrix[candidates_in_data, , drop = FALSE])))
#   )
# 
#   # Aggregate by cell type
#   mean_expr <- cell_type_expr[, lapply(.SD, mean), by = cell_type,
#                                .SDcols = candidates_in_data]
# 
#   # Calculate percentage of cells expressing each gene (> 0) per cell type
#   pct_expr <- cell_type_expr[, lapply(.SD, function(x) mean(x > 0) * 100),
#                               by = cell_type, .SDcols = candidates_in_data]
# 
#   cat("Mean expression by cell type:\n")
#   print(mean_expr)
# 
#   cat("\nPercent cells expressing (> 0) by cell type:\n")
#   print(pct_expr)
# }

## ----add_to_project-----------------------------------------------------------
# # Create a summary table of cell-type specificity for each candidate
# specificity_summary <- data.table(
#   gene_symbol = candidates_in_data,
#   cardiomyocyte_expr = as.numeric(mean_expr[cell_type == "cardiomyocyte",
#                                              ..candidates_in_data]),
#   fibroblast_expr = as.numeric(mean_expr[cell_type == "fibroblast",
#                                           ..candidates_in_data]),
#   endothelial_expr = as.numeric(mean_expr[cell_type == "endothelial cell",
#                                            ..candidates_in_data]),
#   immune_expr = as.numeric(mean_expr[cell_type == "immune cell",
#                                       ..candidates_in_data])
# )
# 
# # Calculate a simple specificity index
# # Ratio of target cell type expression to mean of all other cell types
# specificity_summary[, cardiomyocyte_specificity :=
#   cardiomyocyte_expr / rowMeans(.SD),
#   .SDcols = c("fibroblast_expr", "endothelial_expr", "immune_expr")]
# 
# # Add to locusPackRat project
# addRatTable(
#   data = specificity_summary,
#   table_name = "sc_cell_type_expression",
#   abbreviation = "sc",
#   link_type = "gene",
#   link_by = "gene_symbol",
#   project_dir = "my_qtl_project"
# )

## ----session_info-------------------------------------------------------------
# sessionInfo()


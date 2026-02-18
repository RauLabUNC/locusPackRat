## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  eval = requireNamespace("Seurat", quietly = TRUE) &&
         requireNamespace("ggplot2", quietly = TRUE)
)

## ----libraries----------------------------------------------------------------
# Core packages
library(locusPackRat)
library(data.table)

# scRNA-seq analysis
library(Seurat)     # install from CRAN: install.packages("Seurat")
library(ggplot2)

## ----create_project-----------------------------------------------------------
# Define candidate genes from a cardiac QTL study
cardiac_candidates <- c("Cisd2", "Pdlim5", "Manba", "Fhod3", "Myh7")

# Create a temporary project directory
project_dir <- file.path(tempdir(), "sc_demo_project")

# Initialize a locusPackRat project with these genes
gene_input <- data.frame(gene_symbol = cardiac_candidates)
initPackRat(
  data = gene_input,
  mode = "gene",
  species = "mouse",
  genome = "mm39",
  project_dir = project_dir,
  force = TRUE
)

# Show the initialized gene table
master_genes <- fread(file.path(project_dir, ".locusPackRat", "input", "genes.csv"))
print(master_genes[, .(gene_symbol, chr, start, end)])

## ----load_data----------------------------------------------------------------
# Load the bundled mini Seurat object
seurat_path <- system.file("extdata", "mini_heart_seurat.rds",
                           package = "locusPackRat")
# Fallback for devtools/development builds where inst/ isn't installed yet
if (!nzchar(seurat_path)) {
  seurat_path <- file.path("..", "inst", "extdata", "mini_heart_seurat.rds")
}
heart_seurat <- readRDS(seurat_path)

# Inspect the object
heart_seurat

# Cell type composition
table(heart_seurat$cell_type)

## ----load_tabula_muris, eval = FALSE------------------------------------------
# # Download Tabula Muris heart data (10X Genomics)
# # Data available from: https://tabula-muris.ds.czbiohub.org/
# heart_data <- Read10X(data.dir = "path/to/tabula_muris/Heart-10X_P7_4/")
# heart_seurat <- CreateSeuratObject(counts = heart_data, project = "TabulaMuris_Heart",
#                                    min.cells = 3, min.features = 200)
# heart_seurat <- NormalizeData(heart_seurat)
# heart_seurat <- FindVariableFeatures(heart_seurat, nfeatures = 2000)
# heart_seurat <- ScaleData(heart_seurat)
# heart_seurat <- RunPCA(heart_seurat)
# heart_seurat <- FindNeighbors(heart_seurat, dims = 1:20)
# heart_seurat <- FindClusters(heart_seurat, resolution = 0.5)
# heart_seurat <- RunUMAP(heart_seurat, dims = 1:20)

## ----module_scores------------------------------------------------------------
# Only score genes that are present in the dataset
candidates_in_data <- intersect(cardiac_candidates, rownames(heart_seurat))
cat("Candidates found in dataset:", paste(candidates_in_data, collapse = ", "), "\n")

# Define gene lists for module scoring
candidate_list <- list(cardiac_candidates = candidates_in_data)

# Calculate module scores
# For small gene panels (e.g., bundled demo data), reduce ctrl and nbin
# so that control gene sampling does not exceed available genes per bin
n_features <- nrow(heart_seurat)
use_ctrl <- if (n_features < 200) 5 else 100
use_nbin <- if (n_features < 200) 5 else 24
heart_seurat <- AddModuleScore(
  object = heart_seurat,
  features = candidate_list,
  ctrl = use_ctrl,
  nbin = use_nbin,
  name = "candidate_score"
)

# The score is added as a metadata column: "candidate_score1"
# Higher scores indicate stronger collective expression of the candidate genes
summary(heart_seurat$candidate_score1)

## ----plot_umap, fig.cap = "UMAP colored by candidate gene module score. Warmer colors indicate higher collective expression."----
FeaturePlot(
  heart_seurat,
  features = "candidate_score1",
  cols = c("lightgrey", "blue", "red"),
  pt.size = 1.5
) +
  ggtitle("Candidate Gene Module Score") +
  theme(plot.title = element_text(hjust = 0.5))

## ----plot_violin, fig.cap = "Module score distribution by cell type. Cardiomyocytes show the highest enrichment."----
VlnPlot(
  heart_seurat,
  features = "candidate_score1",
  group.by = "cell_type",
  pt.size = 0
) +
  ggtitle("Candidate Gene Enrichment by Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ----plot_individual_genes, fig.cap = "Dot plot of individual candidate genes by cell type. Size = percent expressing; color = mean expression."----
# Use only candidates present in the data
genes_to_plot <- intersect(cardiac_candidates, rownames(heart_seurat))

DotPlot(
  heart_seurat,
  features = genes_to_plot,
  group.by = "cell_type"
) +
  RotatedAxis() +
  ggtitle("Individual Candidate Gene Expression by Cell Type")

## ----quantify_specificity-----------------------------------------------------
# Extract expression matrix for candidate genes
# Use layer= for Seurat v5 compatibility, fall back to slot= for v4
expr_matrix <- tryCatch(
  GetAssayData(heart_seurat, layer = "data"),
  error = function(e) GetAssayData(heart_seurat, slot = "data")
)
candidates_in_data <- intersect(cardiac_candidates, rownames(expr_matrix))

if (length(candidates_in_data) > 0) {
  # Calculate mean expression per cell type per gene
  cell_type_expr <- data.table(
    cell_type = heart_seurat$cell_type,
    as.data.table(t(as.matrix(expr_matrix[candidates_in_data, , drop = FALSE])))
  )

  # Aggregate by cell type
  mean_expr <- cell_type_expr[, lapply(.SD, mean), by = cell_type,
                               .SDcols = candidates_in_data]

  # Calculate percentage of cells expressing each gene (> 0) per cell type
  pct_expr <- cell_type_expr[, lapply(.SD, function(x) mean(x > 0) * 100),
                              by = cell_type, .SDcols = candidates_in_data]

  cat("Mean expression by cell type:\n")
  print(mean_expr)

  cat("\nPercent cells expressing (> 0) by cell type:\n")
  print(pct_expr)
}

## ----add_to_project-----------------------------------------------------------
# Build a summary table of cell-type expression for each candidate
cell_type_names <- sort(unique(heart_seurat$cell_type))

specificity_summary <- data.table(
  gene_symbol = candidates_in_data,
  cardiomyocyte_expr = as.numeric(mean_expr[cell_type == "cardiomyocyte",
                                             ..candidates_in_data]),
  fibroblast_expr = as.numeric(mean_expr[cell_type == "fibroblast",
                                          ..candidates_in_data]),
  endothelial_expr = as.numeric(mean_expr[cell_type == "endothelial",
                                           ..candidates_in_data]),
  immune_expr = as.numeric(mean_expr[cell_type == "immune",
                                      ..candidates_in_data])
)

# Calculate a simple specificity index:
# ratio of cardiomyocyte expression to mean of all other cell types
specificity_summary[, cardiomyocyte_specificity :=
  cardiomyocyte_expr / rowMeans(.SD),
  .SDcols = c("fibroblast_expr", "endothelial_expr", "immune_expr")]

print(specificity_summary)

# Add to locusPackRat project
addRatTable(
  data = specificity_summary,
  table_name = "sc_cell_type_expression",
  abbreviation = "sc",
  link_type = "gene",
  link_by = "gene_symbol",
  project_dir = project_dir
)

## ----show_project_state-------------------------------------------------------
# Show all tables in the project
listPackRatTables(project_dir)

## ----session_info-------------------------------------------------------------
sessionInfo()

## ----cleanup, include = FALSE-------------------------------------------------
# Clean up the temporary project directory
reg.finalizer(environment(), function(e) {
  unlink(project_dir, recursive = TRUE, force = TRUE)
}, onexit = TRUE)


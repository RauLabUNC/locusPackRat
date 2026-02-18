## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)

## ----install, eval=FALSE------------------------------------------------------
# # Install from GitHub
# #devtools::install_github("RauLabUNC/locusPackRat")
# 
# #Or, install locally
# #devtools::install_local("Path/To/File/locusPackRat.zip")
# 
# # Load required libraries
# library(locusPackRat)
# library(data.table)
# library(jsonlite)

## ----libraries, echo=FALSE----------------------------------------------------
library(locusPackRat)
library(data.table)
library(jsonlite)

## ----gene_mode_init, message=TRUE---------------------------------------------
# Prepare example RNA-seq differential expression results
set.seed(27)
rna_genes <- data.table(
  gene_symbol = c(
    c("Myc", "Tp53", "Egfr", "Vegfa", "Il6", "Tnf", "Apoe", "Lep", "Ins2", "Bdnf"),
    sample(fread(system.file("extdata", "mouse_coords_mm39.csv",
                             package = "locusPackRat"))$gene_symbol, 140)
  ),
  log2FC = round(rnorm(150, mean = 0, sd = 2), 2),
  padj = 10^(-runif(150, 0.5, 10)),
  baseMean = round(10^runif(150, 1, 5), 1),
  direction = ifelse(rnorm(150) > 0, "up", "down")
)

# Initialize project - automatically adds genomic coordinates
initPackRat(
  data = rna_genes,
  mode = "gene",
  species = "mouse",
  genome = "mm39",
  project_dir = "rnaseq_analysis",
  force = TRUE
)

## ----add_tissue_expression----------------------------------------------------
# Load processed genes with coordinates
processed_genes <- fread("rnaseq_analysis/.locusPackRat/input/genes.csv")

# Add tissue expression data (TPM values)
tissue_expression <- data.table(
  gene_symbol = processed_genes$gene_symbol,
  brain = round(abs(rnorm(nrow(processed_genes), 50, 30)), 1),
  liver = round(abs(rnorm(nrow(processed_genes), 40, 25)), 1),
  heart = round(abs(rnorm(nrow(processed_genes), 45, 28)), 1),
  kidney = round(abs(rnorm(nrow(processed_genes), 35, 20)), 1),
  lung = round(abs(rnorm(nrow(processed_genes), 42, 24)), 1),
  muscle = round(abs(rnorm(nrow(processed_genes), 38, 22)), 1)
)

# Calculate tissue specificity
tissue_expression[, max_tissue := names(.SD)[max.col(.SD)],
                 .SDcols = c("brain", "liver", "heart", "kidney", "lung", "muscle")]
tissue_expression[, tissue_specific := apply(.SD, 1, function(x) max(x)/mean(x) > 2),
                 .SDcols = c("brain", "liver", "heart", "kidney", "lung", "muscle")]

addRatTable(
  data = tissue_expression,
  table_name = "tissue_expression",
  abbreviation = "te",
  link_type = "gene",
  link_by = "gene_symbol",
  project_dir = "rnaseq_analysis"
)

## ----add_pathways-------------------------------------------------------------
# Add pathway membership data
pathways <- c("Cell cycle", "Apoptosis", "Immune response", "Metabolism",
              "Signal transduction", "DNA repair", "Protein synthesis")

pathway_data <- data.table(
  gene_symbol = sample(processed_genes$gene_symbol, 80, replace = TRUE),
  pathway = sample(pathways, 80, replace = TRUE),
  evidence = sample(c("experimental", "computational", "literature"), 80, replace = TRUE)
)

addRatTable(
  data = unique(pathway_data),
  table_name = "pathway_annotations",
 # abbreviation = "pa",
  link_type = "gene",
  link_by = "gene_symbol",
  project_dir = "rnaseq_analysis"
)

## ----add_mouseMine_and_Open_Targets_data_gene---------------------------------
queryMouseMine(project_dir = "rnaseq_analysis")
queryOpenTargets(project_dir = "rnaseq_analysis")

## ----gene_examine_package, message=T------------------------------------------
listPackRatTables(project_dir = "rnaseq_analysis", full_info = TRUE)


## ----generate_gene_excel------------------------------------------------------
# Create Excel workbook with multiple filtered views
makeGeneSheet(
  format = "excel",
  output_file = "rnaseq_complete_analysis.xlsx",
  split_by = "criteria",
  prefix_mode = "collision",
    
  split_criteria = list(
    "All_DEGs" = "TRUE",
    "Upregulated" = "direction == 'up'",
    "Downregulated" = "direction == 'down'",
    "High_Confidence" = "padj < 0.001 & abs(log2FC) > 2",
    "Brain_Specific" = "max_tissue == 'brain' & tissue_specific == TRUE",
    "Immune_Genes" = "pathway == 'Immune response'",
    "mouseMine_Immune" = "grepl('immune',phenotype)"
  ),
  highlight_genes = c("Myc", "Tp53", "Egfr", "Il6", "Tnf"),
  include_supplementary = TRUE,
  project_dir = "rnaseq_analysis"
)

## ----region_mode_init_QTL, message=TRUE---------------------------------------
#LocusPackRat includes an example QTL output with a single significant locus on chromsome 3.
#Here we examine two regions, one significant and the other one insignificant

QTL_peaks <- data.table(
  chr = as.character(c(3,6)),
  start=c(126967000,50000000),
  end=c(149860000,60000000),
  peak_id=c("Significant_1","Insignificant_1"),
  max_lod=c(5.64, 3.21)
)

initPackRat(
  data = QTL_peaks,
  mode = "region",
  species = "mouse",
  genome = "mm39",
  project_dir = "qtl_analysis",
  keep_pseudo=F,
  force = TRUE
)



## ----add_QTL_and_Founder_Information------------------------------------------
QTL_file <- system.file("extdata", "sample_scan.csv", package = "locusPackRat")
QTL_Data <- fread(QTL_file)

QTL_Data <- data.table(
  chr = QTL_Data$chr,
  start = QTL_Data$pos,
  end = QTL_Data$pos,
  lod = QTL_Data$lod,
  marker = QTL_Data$marker
)

addRatTable(
  data = QTL_Data,
  table_name = "QTL_results",
  abbreviation = "qr",
  link_type = "region",
  link_by = "chr,start,end",
  project_dir = "qtl_analysis"
)

founder_file <- system.file("extdata", "sample_founders.csv", package = "locusPackRat")
founder_Data <- fread(founder_file)

addRatTable(
  data = founder_Data,
  table_name = "founder_info",
  abbreviation = "fi",
  link_type = "region",
  link_by = "chr,start,end",
  project_dir = "qtl_analysis"
)



## ----add NS_SNP Information---------------------------------------------------

SNP_file <- system.file("extdata", "sample_NS_SNP.csv", package = "locusPackRat")
SNP_Data <- fread(SNP_file)
SNP_Data$start <- SNP_Data$pos
SNP_Data$end <- SNP_Data$pos

addRatTable(
  data = SNP_Data,
  table_name = "NS_SNP_information",
  abbreviation = "NSi",
  link_type = "region",
  link_by = "chr,start,end",
  project_dir = "qtl_analysis"
)



## ----add_mouseMine_and_Open_Targets_data_QTL----------------------------------
queryMouseMine(project_dir = "qtl_analysis")
queryOpenTargets(project_dir = "qtl_analysis")

## ----gene_examine_package_QTL, message=T--------------------------------------
listPackRatTables(project_dir = "qtl_analysis", full_info = TRUE)


## ----generate_QTL_plot--------------------------------------------------------
generateLocusZoomPlot(
     region_id="region_1",
     project_dir="qtl_analysis",
     scan_table="QTL_results",
     signal_table="founder_info",
     width=10,
     height=6,
     threshold=4,
     layout_ratios = c(manhattan = 0.35, signal = 0.40, genes = 0.25))


## ----generate_region_excels_QTL-----------------------------------------------
# Peak-centric analysis workbook
makeGeneSheet(
  format = "excel",
  output_file = "qtl_analysis.xlsx",
  prefix_mode = "abbreviated",
  split_by = "criteria",
  split_criteria = list(
    "All_Results" = "TRUE",
    "Chr 3 Peak" = "peak_id == 'Significant_1'",
    "Mouse_Phenotype" = "grepl('[c|C]ardiac',mm_phenotype) | grepl('[h|H]eart',mm_phenotype) ",
    "Sig OTC" = "(otc_score>1 | otc_score< -1) & peak_id == 'Significant_1'"
  ),
  include_supplementary =  TRUE,
  exclude_tables=c("founder_info","QTL_results"),
  project_dir = "qtl_analysis"
)


## ----region_mode_init_ATAC, message=TRUE--------------------------------------
# Prepare example ATAC-seq peaks
set.seed(123)
chromosomes <- c(as.character(1:19), "X", "Y")

atac_peaks <- data.table(
  chr = sample(chromosomes, 200, replace = TRUE,
               prob = c(rep(0.06, 19), 0.07, 0.07)),
  start = sample(1:150000000, 200),
  peak_id = paste0("peak_", 1:200),
  fold_enrichment = round(runif(200, 2, 50), 2),
  qvalue = 10^(-runif(200, 1, 20))
)

# Calculate end positions (typical ATAC peak width: 200-1000bp)
atac_peaks[, end := start + sample(200:1000, .N, replace = TRUE)]
atac_peaks[, summit := start + round((end - start) / 2)]
atac_peaks[, peak_type := sample(c("promoter", "enhancer", "intergenic"),
                                 .N, replace = TRUE,
                                 prob = c(0.3, 0.4, 0.3))]

# Initialize region-based project
initPackRat(
  data = atac_peaks,
  mode = "region",
  species = "mouse",
  genome = "mm39",
  project_dir = "atacseq_analysis",
  force = TRUE
)

## ----add_chromatin------------------------------------------------------------
# Load processed regions
regions <- fread("atacseq_analysis/.locusPackRat/input/regions.csv")

# Add chromatin state annotations
chromatin_states <- data.table(
  chr = regions$chr,
  start = regions$start,
  end = regions$end,
  state = sample(c("Active_Promoter", "Strong_Enhancer", "Weak_Enhancer",
                  "Poised_Promoter", "Repressed", "Heterochromatin", "Transcribed"),
                size = nrow(regions), replace = TRUE,
                prob = c(0.15, 0.2, 0.15, 0.1, 0.1, 0.1, 0.2)),
  cell_type = sample(c("ES_cells", "Neural", "Hepatocyte", "T_cells"),
                    size = nrow(regions), replace = TRUE)
)

addRatTable(
  data = chromatin_states,
  table_name = "chromatin_states",
  abbreviation = "cs",
  link_type = "region",
  link_by = "chr,start,end",
  project_dir = "atacseq_analysis"
)

## ----add_tf_motifs------------------------------------------------------------
# Add TF binding motif data
tf_list <- c("CTCF", "NFKB", "AP1", "ETS1", "GATA1", "SOX2", "OCT4", "NANOG",
            "MYC", "MAX", "STAT3", "CREB", "SP1", "YY1", "REST")

# Create random starts
random_starts <- sample(regions$start, 120, replace = TRUE)

motif_data <- data.table(
  chr = sample(regions$chr, 120, replace = TRUE),
  start = random_starts,
  end = random_starts + sample(50:200, 120, replace = TRUE),
  tf_name = sample(tf_list, 120, replace = TRUE),
  motif_score = round(runif(120, 6, 15), 2),
  p_value = 10^(-runif(120, 2, 8))
)

# Add additional motifs for some peaks
additional_motifs <- motif_data[sample(1:nrow(motif_data), 30)]
additional_motifs[, tf_name := sample(tf_list, .N, replace = TRUE)]
motif_data <- rbind(motif_data, additional_motifs)

addRatTable(
  data = motif_data,
  table_name = "tf_motifs",
  abbreviation = "tf_motifs",
  link_type = "region",
  link_by = "chr,start,end",
  project_dir = "atacseq_analysis"
)

## ----add_differential---------------------------------------------------------
# Add differential accessibility between conditions
diff_access <- data.table(
  chr = regions$chr,
  start = regions$start,
  end = regions$end,
  condition_A_signal = round(abs(rnorm(nrow(regions), 100, 50)), 1),
  condition_B_signal = round(abs(rnorm(nrow(regions), 100, 50)), 1),
  log2FC = round(rnorm(nrow(regions), 0, 1.5), 2),
  padj = 10^(-runif(nrow(regions), 0, 8)),
  higher_in = ifelse(rnorm(nrow(regions)) > 0, "condition_B", "condition_A"),
  significant = runif(nrow(regions)) < 0.3
)

addRatTable(
  data = diff_access,
  table_name = "differential_accessibility",
  abbreviation = "da",
  link_type = "region",
  link_by = "chr,start,end",
  project_dir = "atacseq_analysis"
)

## ----add_mouseMine_and_Open_Targets_data_ATAC---------------------------------
queryMouseMine(project_dir = "atacseq_analysis")
queryOpenTargets(project_dir = "atacseq_analysis")

## ----gene_examine_package_ATAC, message=T-------------------------------------
listPackRatTables(project_dir = "atacseq_analysis", full_info = TRUE)


## ----generate_region_excels---------------------------------------------------
# Peak-centric analysis workbook
makeGeneSheet(
  format = "excel",
  output_file = "atacseq_peak_analysis.xlsx",
  split_by = "criteria",
  split_criteria = list(
    "All_Peaks" = "TRUE",
    "High_Enrichment" = "fold_enrichment > 20",
    "Promoter_Peaks" = "peak_type == 'promoter'",
    "Enhancer_Peaks" = "peak_type == 'enhancer'",
    "Active_Chromatin" = "state %in% c('Active_Promoter', 'Strong_Enhancer')",
    "Differential" = "significant == TRUE & abs(log2FC) > 1"
  ),
  include_supplementary = TRUE,
  project_dir = "atacseq_analysis"
)

# TF-centric analysis workbook
makeGeneSheet(
  format = "excel",
  output_file = "atacseq_tf_binding.xlsx",
  split_by = "criteria",
  split_criteria = list(
    "CTCF_Binding" = "tf_name == 'CTCF'",
    "Pioneer_TFs" = "tf_name %in% c('SOX2', 'OCT4', 'NANOG')",
    "Immune_TFs" = "tf_name %in% c('NFKB', 'AP1', 'STAT3')",
    "High_Motif_Score" = "motif_score > 12"
  ),
  include_supplementary = c("tf_motifs", "chromatin_states"),
  project_dir = "atacseq_analysis"
)

## ----build_packet_qtl, message=TRUE-------------------------------------------
# Package the QTL analysis for sharing
buildPacket(project_dir = "qtl_analysis", overwrite = TRUE)

## ----build_packet_custom, eval=FALSE------------------------------------------
# # Write to a specific directory with a custom name
# buildPacket(
#   project_dir = "rnaseq_analysis",
#   output_path = "~/shared_results",
#   filename = "rnaseq_v2_packet.zip",
#   overwrite = TRUE
# )
# 
# # Skip the README if you prefer
# buildPacket(project_dir = "atacseq_analysis", include_readme = FALSE)

## ----build_packet_supplementary, eval=FALSE-----------------------------------
# # Include all supplementary CSVs
# buildPacket(
#   project_dir = "qtl_analysis",
#   include_supplementary = TRUE,
#   overwrite = TRUE
# )
# 
# # Include only specific supplementary tables
# buildPacket(
#   project_dir = "qtl_analysis",
#   include_supplementary = c("mouseMine", "openTargets_diseases"),
#   overwrite = TRUE
# )

## ----crossspecies_analysis, message=TRUE--------------------------------------
# Define human cancer genes
cancer_genes <- c("BRCA1", "BRCA2", "TP53", "EGFR", "KRAS", "MYC",
                 "PTEN", "RB1", "VHL", "APC", "ATM", "CDKN2A",
                 "MLH1", "MSH2", "BRAF", "PIK3CA", "ERBB2", "ALK")

human_genes <- data.table(
  gene_symbol = cancer_genes,
  cancer_type = c(rep("breast", 2), "multiple", "lung", "colorectal", "multiple",
                  "brain", "retinoblastoma", "kidney", "colorectal", "breast", "melanoma",
                  rep("colorectal", 2), "melanoma", "multiple", "breast", "lung"),
  mutation_frequency = round(runif(length(cancer_genes), 5, 80), 1)
)

# Initialize with automatic orthology mapping
initPackRat(
  data = human_genes,
  mode = "gene",
  species = "human",
  genome = "hg38",
  project_dir = "cancer_genes",
  force = TRUE
)

# Add clinical significance data
clinical_data <- data.table(
  gene_symbol = cancer_genes,
  pathogenicity = sample(c("pathogenic", "likely_pathogenic", "VUS"),
                        length(cancer_genes), replace = TRUE,
                        prob = c(0.5, 0.3, 0.2)),
  FDA_approved_drug = sample(c(TRUE, FALSE), length(cancer_genes),
                            replace = TRUE, prob = c(0.3, 0.7)),
  clinical_trials = sample(0:10, length(cancer_genes), replace = TRUE)
)

addRatTable(
  data = clinical_data,
  table_name = "clinical_significance",
  link_type = "gene",
  link_by = "gene_symbol",
  project_dir = "cancer_genes"
)

# Generate cross-species analysis workbook
makeGeneSheet(
  format = "excel",
  output_file = "cancer_genes_analysis.xlsx",
  split_by = "criteria",
  split_criteria = list(
    "High_Frequency" = "mutation_frequency > 30",
    "FDA_Targets" = "FDA_approved_drug == TRUE",
    "Breast_Cancer" = "cancer_type == 'breast'",
    "With_Mouse_Ortholog" = "!is.na(mouse_gene_symbol)"
  ),
  include_supplementary = TRUE,
  project_dir = "cancer_genes"
)

## ----liftOver, message=F, eval=requireNamespace("liftOver", quietly=TRUE)-----
# library(liftOver)
# mm9_ncRNA = read.csv(system.file("extdata", "sample_mm9_ncRNA.csv", package = "locusPackRat"))
# chain = import.chain( system.file("extdata", "mm9ToMm39.over.chain", package = "locusPackRat"))
# 
# #convert table to GRanges object.  Note that the chain files normally want
# #chromosomes in 'chr1' format vs '1' format.  So you will need to convert.
# mm9_ncRNA$chr=paste0("chr",mm9_ncRNA$chr)
# mm9_ncRNA=makeGRangesFromDataFrame(mm9_ncRNA,keep.extra.columns = T)
# 
# mm39_ncRNA=unlist(liftOver(mm9_ncRNA,chain = chain))
# mm39_ncRNA=as.data.frame(mm39_ncRNA)
# mm39_ncRNA$seqnames=gsub("chr","",mm39_ncRNA$seqnames)
# 
# #Your data are now ready to be added to a locusPackRat project
# 
# 

## ----list_tables, message=T---------------------------------------------------
# View supplementary tables in each project
listPackRatTables("rnaseq_analysis", full_info = TRUE)
listPackRatTables("atacseq_analysis", full_info = FALSE)

## ----remove_table-------------------------------------------------------------

removeRatTable("chromatin_states", project_dir = "atacseq_analysis")


## ----show_config, eval=FALSE--------------------------------------------------
# # Load project configuration
# config <- jsonlite::read_json("rnaseq_analysis/.locusPackRat/config.json")
# str(config)

## ----csv_output, eval=FALSE---------------------------------------------------
# makeGeneSheet(
#   format = "csv",
#   output_file = "results.csv",
#   include_supplementary = TRUE,
#   project_dir = "rnaseq_analysis"
# )

## ----filtered_output, eval=FALSE----------------------------------------------
# makeGeneSheet(
#   format = "csv",
#   output_file = "significant_genes.csv",
#   filter_expr = "padj < 0.05 & abs(log2FC) > 1",
#   include_supplementary = TRUE,
#   project_dir = "rnaseq_analysis"
# )

## ----selective_tables, eval=FALSE---------------------------------------------
# makeGeneSheet(
#   format = "excel",
#   output_file = "selected_data.xlsx",
#   include_supplementary = c("tissue_expression", "pathway_annotations"),
#   project_dir = "rnaseq_analysis"
# )

## ----liftover_example, eval=FALSE---------------------------------------------
# # Example: Convert mm10 coordinates to mm39
# library(rtracklayer)
# library(GenomicRanges)
# 
# # Load a chain file for the conversion
# # Chain files are available from UCSC:
# # https://hgdownload.cse.ucsc.edu/goldenpath/mm10/liftOver/mm10ToMm39.over.chain.gz
# chain <- import.chain("mm10ToMm39.over.chain")
# 
# # Create a GRanges object from your data
# # Suppose you have eQTL data in mm10 coordinates
# eqtl_mm10 <- data.table(
#   chr = c("chr3", "chr3", "chr3"),
#   start = c(131000000, 135000000, 142000000),
#   end = c(131001000, 135001000, 142001000),
#   gene = c("Cisd2", "Manba", "Pdlim5"),
#   lod = c(5.2, 4.8, 3.1)
# )
# 
# gr_mm10 <- GRanges(
#   seqnames = eqtl_mm10$chr,
#   ranges = IRanges(start = eqtl_mm10$start, end = eqtl_mm10$end),
#   gene = eqtl_mm10$gene,
#   lod = eqtl_mm10$lod
# )
# 
# # Perform liftOver
# gr_mm39 <- liftOver(gr_mm10, chain)
# 
# # Convert back to data.table — liftOver returns a GRangesList
# # (one element per input range, since some ranges may map to multiple locations)
# gr_mm39_unlisted <- unlist(gr_mm39)
# 
# eqtl_mm39 <- data.table(
#   chr = as.character(seqnames(gr_mm39_unlisted)),
#   start = start(gr_mm39_unlisted),
#   end = end(gr_mm39_unlisted),
#   gene = gr_mm39_unlisted$gene,
#   lod = gr_mm39_unlisted$lod
# )
# 
# # Remove "chr" prefix if your project uses UCSC-style names without it
# eqtl_mm39[, chr := gsub("^chr", "", chr)]
# 
# # Now add to your mm39 locusPackRat project
# addRatTable(
#   data = eqtl_mm39,
#   table_name = "eqtl_lifted",
#   abbreviation = "eql",
#   link_type = "region",
#   link_by = "chr,start,end",
#   project_dir = "qtl_analysis"
# )

## ----session------------------------------------------------------------------
sessionInfo()


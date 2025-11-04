# genePackRat

R package for generating standardized gene packets from QTL/GWAS loci.

## What It Does

When you map QTLs in mice, you get broad loci (often 10-50 Mb) containing dozens of genes. genePackRat creates standardized "packets" for each locus - multi-sheet Excel files, locus zoom plots, and README summaries - to help your team systematically evaluate candidates.

Each packet consolidates:
- Gene annotations with expression levels
- Human disease associations (Open Targets)
- Mouse phenotypes (MGI)
- QTL mapping visualizations
- Collaborative Cross founder variants (optional)

## Installation

Requires R >= 4.3.0

```r
# Install BiocManager if needed (for plotgardener dependency)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install from GitHub with all dependencies
BiocManager::install("RauLabUNC/genePackRat")

# Or using devtools
devtools::install_github("RauLabUNC/genePackRat")
```

Core dependencies install automatically: data.table, dplyr, openxlsx, tidyr, plotgardener (from Bioconductor), and RColorBrewer.

## Basic Usage

```r
library(genePackRat)

# Load your QTL results
sig_regions <- read.csv("data/qtl_significant_regions.csv")
scan_data <- readRDS("data/qtl_scans.rds")
thresholds <- readRDS("data/qtl_thresholds.rds")

# Generate packet for first locus
generateLocusPacket(
  locus_cluster = sig_regions[1, ],
  input_path = "data/",
  output_path = "results/locus_packets/",
  scan_data = scan_data,
  threshold_data = thresholds
)
```

This creates:
```
results/locus_packets/locus_chr1_100-110Mb/
├── gene_info_cluster_chr1_100-110Mb.xlsx    # 3-sheet workbook
├── zoomPlots/
│   └── locus_zoom_HR_Ctrl_chr1_105.2Mb.pdf
├── README_summary.txt                        # Cardiac genes highlighted
└── founder_snp_table_chr1_100-110Mb.csv     # Optional
```

## Core Functions

**`generateLocusPacket()`** - Creates complete packet directory
- Excel workbook with gene annotations
- Locus zoom plots for each QTL
- README summary with cardiac associations
- Founder SNP table (optional)

**`generateLocusZoomPlot()`** - Creates multi-panel QTL visualization
- Manhattan plot with LOD scores and significance thresholds
- Founder allele effects for 8 CC strains with preset colors
- Overlapping QTL display for up to 8 loci
- Gene annotations (works without TxDb/org.Mm packages)

**`generateGeneInfoExcel()`** - 3-sheet workbook
- Sheet 1: AllGenesInCluster (main summary)
- Sheet 2: AllMousePhenotypes (MGI data)
- Sheet 3: AllHumanDiseases (Open Targets)

## Required Data Structure

Test data is included in `tests/testthat/fixtures/` for development.

Production data expects CSV files in this structure:

```
data/
└── processed/joinLoci/
    ├── relational_tables/
    │   ├── genes_mouse.csv              # Gene coordinates
    │   ├── orthology.csv                # Mouse-human orthologs
    │   ├── associations.csv             # Human disease associations
    │   ├── mouseGenePhenotypes.csv      # MGI phenotypes
    │   └── traitLoci.csv                # QTL definitions
    ├── geneTables/
    │   └── multTrait_cis-eQTL_nrvmExp.csv  # Merged gene summary
    └── bulk_exp/
        └── rna_expression.csv           # Bulk RNA-seq (VST)
```

See documentation for required column schemas.

## Customizing File Paths

All paths can be overridden:

```r
generateLocusPacket(
  locus_cluster = my_locus,
  input_path = "data/",
  rna_info_file = "data/rnaseq/VST_Info_2024-12-15.csv",  # Custom path
  founder_mutations_file = NULL,  # Skip CC variants
  heart_pattern = NULL  # Use default cardiac keywords
)
```

## Additional Utilities

**`filterGenes()`** - Simple relational filtering
```r
# Keep only genes present in QTL results
filtered <- filterGenes(
  inputTable = all_genes,
  referenceTable = qtl_genes,
  by = "gene_id",
  joinType = "semi"
)

# Add filters
filtered <- filterGenes(
  inputTable = filtered,
  filters = list(
    makeFilter("biotype", "==", "protein_coding"),
    makeFilter("expression", ">", 50)
  )
)
```

## Citation

Gural B, Kimball T, Luu A, Rau CD. genePackRat: Standardized gene packets for QTL follow-up. *In preparation* (2025).

## License

MIT

## Contact

bgural@unc.edu

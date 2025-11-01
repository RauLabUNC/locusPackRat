# genePackRat

R package for gene prioritization from GWAS/QTL studies using relational filtering.

## Problem

GWAS in model organisms often produces large loci (10-100+ Mb) containing dozens to hundreds of genes. Researchers need to integrate multiple data sources - expression, variants, phenotypes, orthologs - to identify likely causal genes. genePackRat provides tools for flexible data integration and filtering based on relational database principles.

## Installation

Requires R >= 4.3.0

```r
# Install from GitHub
install.packages("devtools")
devtools::install_github("RauLabUNC/genePackRat")

# Or using remotes
install.packages("remotes")
remotes::install_github("RauLabUNC/genePackRat")
```

Dependencies (data.table, dplyr, openxlsx) will be installed automatically.

## Usage

```r
library(genePackRat)

# Load your data
genes <- read.csv("my_genes.csv")
expression <- read.csv("expression_data.csv")
gwas_hits <- read.csv("gwas_significant_genes.csv")

# Filter for GWAS hits only
candidate_genes <- filterGenes(
  inputTable = genes,
  referenceTable = gwas_hits,
  by = "gene_id",
  joinType = "semi"  # keeps only matching genes
)

# Add expression data and filter
expressed_candidates <- filterGenes(
  inputTable = candidate_genes,
  referenceTable = expression,
  by = "gene_id",
  joinType = "left",
  filters = list(
    makeFilter("TPM_heart", condition = ">", value = 100),
    makeFilter("biotype", condition = "==", value = "protein_coding")
  )
)

# Build gene table with multiple annotations
gene_table <- build_gene_table(
  genes_df = expressed_candidates,
  expression_df = expression,
  orthology_df = orthologs
)

# Export to Excel
create_gene_workbook(
  gene_data = gene_table,
  output_file = "candidate_genes.xlsx"
)
```

## Functions

### Filtering (filterGenes.R)
- `filterGenes()` - Filter data using joins and custom criteria
- `makeFilter()` - Create filter specifications
- `applyFilter()` - Apply filter to data

### Data Integration (join_tables.R)
- `build_gene_table()` - Combine annotations from multiple sources
- `aggregate_genes_by_locus()` - Group genes by genomic regions
- `filter_genes()` - Legacy filtering function

### Excel Export (excel.R)
- `create_gene_workbook()` - Create multi-sheet Excel files

### Phenotypes (phenotype.R)
- `extract_phenotypes()` - Extract phenotypes by keyword
- `score_phenotype_relevance()` - Rank phenotypes
- `summarize_phenotypes_by_gene()` - Aggregate by gene

### Reporting (summary.R)
- `generate_locus_report()` - Create summary reports
- `identify_top_candidates()` - Prioritize genes
- `generate_recommendations()` - Suggest follow-ups

### Plotting (plotting.R)
- `plot_locus()` - Create locus plots
- `plot_locus_dashboard()` - Multi-panel plots


## Join Types

- `semi` - Keep rows that match reference (like `%in%`)
- `anti` - Keep rows that DON'T match reference (like `!%in%`)
- `inner` - Merge matching rows from both tables
- `left` - Keep all input rows, add reference columns

Column mapping for different names:
```r
filtered <- filterGenes(
  inputTable = mouse_genes,
  referenceTable = human_orthologs,
  by = c("mouse_symbol" = "gene"),  # Map columns
  joinType = "inner"
)
```



## Example Workflow

```r
library(genePackRat)

# Load data
all_genes <- read.csv("all_genes.csv")
gwas_results <- read.csv("gwas_significant.csv")
rnaseq <- read.csv("heart_rnaseq.csv")
eqtl_data <- read.csv("heart_eqtl.csv")

# Step 1: Get GWAS hits
gwas_genes <- filterGenes(
  inputTable = all_genes,
  referenceTable = gwas_results,
  by = "gene_id",
  joinType = "semi"
)

# Step 2: Add expression and filter
expressed_hits <- filterGenes(
  inputTable = gwas_genes,
  referenceTable = rnaseq,
  by = "gene_id",
  joinType = "left",
  filters = list(
    makeFilter("TPM", condition = ">", value = 10),
    makeFilter("padj", condition = "<", value = 0.05)
  )
)

# Step 3: Add eQTL data
with_eqtl <- filterGenes(
  inputTable = expressed_hits,
  referenceTable = eqtl_data,
  by = c("gene_id" = "gene"),
  joinType = "left",
  filters = list(
    makeFilter("eqtl_pvalue", condition = "<", value = 1e-5)
  )
)

# Export results
create_gene_workbook(with_eqtl, "candidates.xlsx")
```


## Testing

```r
# Run tests
cd tests
Rscript run_tests.R

# Or generate HTML report
R -e "rmarkdown::render('test_package_functionality.Rmd')"
```

## Citation

Gural B, Kimball T, Luu A, Rau CD. genePackRat: A framework for prioritizing candidate genes from GWAS intervals. *In preparation* (2025).

## License

MIT License

## Contact

bgural@unc.edu
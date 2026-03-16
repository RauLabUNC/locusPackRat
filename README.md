# locusPackRat

An R package for organizing genomic analysis projects. You give it a set of genes or regions, attach any supplementary data you want, query external databases, and get back a merged Excel sheet (or CSV) ready for downstream work.

Built for QTL/GWAS follow-up, but works for any gene- or region-level study.

## Installation

```r
# From GitHub (requires devtools)
devtools::install_github("RauLabUNC/locusPackRat")
```

Optional Bioconductor packages for locus zoom plots (install for your genome):

```r
# Mouse mm39
BiocManager::install(c("plotgardener", "TxDb.Mmusculus.UCSC.mm39.knownGene", "org.Mm.eg.db"))

# Mouse mm10
BiocManager::install(c("plotgardener", "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db"))

# Human hg38
BiocManager::install(c("plotgardener", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"))

# Human hg19
BiocManager::install(c("plotgardener", "TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"))
```

## Quick start

All functions have short `lpr_*` aliases (e.g., `lpr_init` for `initPackRat`) for easier discovery and tab-completion. Both names work identically.

```r
library(locusPackRat)

# Check your setup
lpr_check()

# Initialize a project
lpr_init(
  data = data.frame(gene_symbol = c("Myc", "Tp53", "Egfr"),
                    log2FC = c(2.5, -1.8, 3.2)),
  mode = "gene", species = "mouse", genome = "mm39",
  project_dir = "my_analysis"
)

# Attach your own data
lpr_add_table(my_annotations, table_name = "pathway_data",
              link_type = "gene", link_by = "gene_symbol",
              project_dir = "my_analysis")

# Pull from external databases
lpr_query_mousemine(project_dir = "my_analysis")
lpr_query_ot(project_dir = "my_analysis")

# Export everything
lpr_export(format = "excel", output_file = "results.xlsx",
           project_dir = "my_analysis")

# Check project status
lpr_check("my_analysis")
```

### Vignettes

See the package vignettes for full walkthroughs:
- `locusPackRat_workflow` — main tutorial (gene mode with RNA-seq data)
- `region_qtl_opentargets` — QTL analysis with Open Targets
- `genenetwork_qtl` — QTL data from GeneNetwork2
- `single_cell_integration` — cell-type-specific evaluation with scRNA-seq

## Supported genomes

Human (hg38, hg19) and mouse (mm39, mm10). We plan to add other species — let us know what would be useful.

## Citation

Gural B, Kimball T, Luu A, Rau CD. locusPackRat: A Flexible Framework for Prioritizing Candidate Genes from GWAS and other Gene-Level Studies. *Under Review* (2025).

## License

MIT

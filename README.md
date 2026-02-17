# locusPackRat

An R package for organizing genomic analysis projects. You give it a set of genes or regions, attach any supplementary data you want, query external databases, and get back a merged Excel sheet (or CSV) ready for downstream work.

Built for QTL/GWAS follow-up, but works for any gene- or region-level study.

## Installation

```r
# From GitHub (requires devtools)
devtools::install_github("RauLabUNC/locusPackRat")
```

Bioconductor packages for locus zoom plots (optional):

```r
BiocManager::install(c("plotgardener", "TxDb.Mmusculus.UCSC.mm39.knownGene", "org.Mm.eg.db"))
```

## Quick start

```r
library(locusPackRat)

# Initialize a project
initPackRat(
  data = data.frame(gene_symbol = c("Myc", "Tp53", "Egfr"),
                    log2FC = c(2.5, -1.8, 3.2)),
  mode = "gene", species = "mouse", genome = "mm39",
  project_dir = "my_analysis"
)

# Attach your own data
addRatTable(my_annotations, table_name = "pathway_data",
            link_type = "gene", link_by = "gene_symbol",
            project_dir = "my_analysis")

# Pull from external databases
queryMouseMine(project_dir = "my_analysis")
queryOpenTargets(project_dir = "my_analysis")

# Export everything
makeGeneSheet(format = "excel", output_file = "results.xlsx",
              project_dir = "my_analysis")
```

See the package vignettes for full walkthroughs (`region_qtl_opentargets`, `genenetwork_qtl`, `single_cell_integration`).

## Supported genomes

Human (hg38, hg19) and mouse (mm39, mm10). We plan to add other species — let us know what would be useful.

## Citation

Gural B, Kimball T, Luu A, Rau CD. locusPackRat: A Flexible Framework for Prioritizing Candidate Genes from GWAS and other Gene-Level Studies. *Under Review* (2025).

## License

MIT

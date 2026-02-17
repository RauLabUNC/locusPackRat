# locusPackRat 0.6.1

## Bug Fixes
- **Known Drugs query**: Fixed invalid `mechanismsOfAction { actionType }` field in GraphQL query; now uses `mechanismOfAction` at the KnownDrug row level (matches current Open Targets API schema)
- **Pharmacogenomics query**: Fixed `variantFunctionalConsequence` (now `SequenceOntologyTerm` object) and drug fields (now nested under `drugs: [DrugWithIdentifiers!]!`); rewrote parser for correct extraction
- **Expression parser**: Fixed `lapply` iterating over data.frame columns instead of rows when `fromJSON` auto-simplifies; now produces ~180 tissue rows per gene instead of 3 empty rows
- **Mouse Phenotypes parser**: Applied same data.frame-to-list-of-rows guard to prevent column iteration bug
- **DepMap parser**: Fixed column explosion (177 wide-format columns instead of 6) by extracting only expected columns from auto-simplified data.frame
- **QTL credible set gene linking**: Expanded `id_map` to include all ensembl IDs from QTL response (not just project genes), fixing 0/179 match issue in region mode

## Improvements
- Reduced batch size to 5 genes when querying heavy data types (mouse_phenotypes, known_drugs) to avoid API timeouts
- Increased API timeout from 60s to 120s for Open Targets queries
- Audit script now uses name-based file lookup instead of newest-file heuristic, includes row-count validation and locus zoom plot test
- Vignette `genenetwork_qtl`: defensive `intersect(names(...), expected_cols)` column subsetting
- Vignette `region_qtl_opentargets`: locus zoom chunks conditionally evaluated when Bioconductor annotation packages are available
- Vignette `locusPackRat_workflow`: synced VignetteIndexEntry with YAML title

## Package Metadata
- Added `Seurat`, `TxDb.Hsapiens.UCSC.hg38.knownGene`, `org.Hs.eg.db` to Suggests
- Added `chr_clean` and `study_locus_id` to `globalVariables()` (fixes R CMD check NOTEs)
- Added `.specify` and `.locusPackRat` patterns to `.Rbuildignore`

---

# locusPackRat 0.6.0

## New Features
- `queryOpenTargetsQTL()`: Query eQTL, pQTL, and sQTL credible sets from the Open Targets Genetics portal, with locus-to-gene (L2G) prediction scores
- Expanded `queryOpenTargets()` data types: now supports `diseases`, `constraints`, `tractability`, `expression`, `mouse_phenotypes`, `homologues`, `interactions`, `known_drugs`, `depmap`, and `pharmacogenomics`
- Genome build messaging: `initPackRat()` now displays a prominent note reminding users that all data must match the project's genome build, with guidance on using `rtracklayer::liftOver()` for coordinate conversion

## New Vignettes
- `region_qtl_opentargets`: Region-level QTL analysis using Open Targets eQTL/pQTL data
- `genenetwork_qtl`: GeneNetwork2 integration for BXD/DO/HMDP QTL studies
- `single_cell_integration`: Single-cell RNA-seq integration via Tabula Muris and Seurat

## Removed
- Removed `ComplexUpset` dependency (replaced with base R alternatives)

## Improvements
- Added plotgardener citation to locus zoom plot output
- `listPackRatTables(full_info = TRUE)` now reports column types, value ranges, and sample counts

---

# locusPackRat 0.5.0

## Breaking Changes
- `listPackRatTables()`: parameter `fullInfo` renamed to `full_info` for API consistency

## Improvements
- Testing infrastructure with testthat
- Documentation examples for all exported functions
- Input validation improvements in `addRatTable()` (validates table_name characters)
- Match statistics reporting in `initPackRat()` (shows matched/unmatched genes)
- GitHub Actions CI for automated R CMD check

---

# locusPackRat 0.4.0

## New Features
- Overlap detection in `addRatTable()` with `overlap_mode` parameter for region-based linking
- `removeRatTable()` function for removing supplementary tables from projects
- Column abbreviation feature for Excel output prefixes (avoids column name collisions)
- Legacy genome support: hg19 and mm10 now available alongside hg38 and mm39
- `queryMouseMine()` now uses httr for more reliable API communication
- Automatic formatting for MouseMine and OpenTargets data in Excel sheet outputs

## Bug Fixes
- Row stripping in `makeGeneSheet()` only occurs with >= 2 rows present
- Fix for handling multiple regions in `generateLocusZoomPlot()`
- Rare error fix in OpenTargets query response handling

## Documentation
- Updated QTL mapping vignette
- Improved function documentation and column name descriptions

---

# locusPackRat 0.3.0

- Initial support for locus zoom plot generation via plotgardener
- Query functions for external databases (MouseMine, OpenTargets)
- Project-based workflow with persistent storage

---

# locusPackRat 0.2.0

- Project directory structure generation
- Vignette for basic workflow

---

# locusPackRat 0.1.0

- Initial package release
- Core functions: `initPackRat()`, `addRatTable()`, `makeGeneSheet()`
- Support for human and mouse genomes (hg38, mm39)

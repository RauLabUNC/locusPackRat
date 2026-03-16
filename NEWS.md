# locusPackRat 0.7.0

## New Features
- `lpr_check()` / `checkPackRat()`: diagnostic function that validates project state, checks file integrity, and reports optional dependency status
- `lpr_*` aliases for all 12 exported functions (e.g., `lpr_init`, `lpr_export`, `lpr_query_ot`) — easier tab-completion and discoverability; both CamelCase and `lpr_*` names work identically
- New vignette: `new_ux_features` demonstrating check and alias usage

## Improvements
- Moved `plotgardener` and `httr` from Imports to Suggests with `requireNamespace()` guards, reducing the install footprint for users who don't need locus zoom plots or MouseMine queries
- `buildPacket()` zip creation now tries the `zip` R package first, then falls back to system `zip`, with an informative error if neither is available
- Expanded README with install instructions for all 4 genome builds (hg38, hg19, mm39, mm10) and `lpr_*` quick-start examples

## Internal
- Consolidated duplicate helpers (`.render_genes`, `.standardize_chrom`) and `%||%` operator into single definitions
- Added `zip` to Suggests

---

# locusPackRat 0.6.3

## Bug Fixes
- **Locus zoom chromosome matching**: Added `.standardize_chrom()` helper that normalizes chromosome columns in scan and signal tables before rendering. Strips "chr" prefix from `chrom` values so they match `pgParams` convention (derived from `regions.csv`). Previously, data with "chr3" in `chrom` would silently produce an empty Manhattan/signal panel because it didn't match "3" from `pgParams`.
- **Missing `.standardize_chrom()` function**: The function was referenced in `generateLocusZoomPlot()` but absent from source, causing errors on fresh installs. Now defined in `plotZoom.R`.

## Improvements
- `generateLocusZoomPlot()` now checks for `GenomicFeatures` package availability before attempting gene track rendering, with a clear install message on failure
- Added `GenomicFeatures` to Suggests in DESCRIPTION

## Notes
- Users with `Seqinfo` >= 1.0.0 (Bioconductor 3.21+) should ensure `GenomeInfoDb` >= 1.46.0 is installed. Older `GenomeInfoDb` versions conflict with the `Seqinfo` package split, causing `plotGenes()` to fail with an opaque `seqinfo` method dispatch error. Upgrade via: `mamba install -c bioconda bioconductor-genomeinfodb=1.46.2` or `BiocManager::install("GenomeInfoDb")`

---

# locusPackRat 0.6.2

## New Features
- `buildPacket()`: new exported function that zips the contents of `.locusPackRat/output/` into a single shareable archive with an auto-generated `README.md` summarizing project metadata, supplementary data sources, and included files
- `buildPacket(include_supplementary)`: optionally bundles raw supplementary CSVs from `.locusPackRat/supplementary/` into the zip. Accepts `FALSE` (default), `TRUE` (all tables), or a character vector of specific table names. Files are placed in a `supplementary/` subdirectory. Warns on missing table names.
- Auto-generated README in packets includes an "Included Supplementary Files" table when supplementary CSVs are bundled

## Improvements
- Vignette `genenetwork_qtl`: chunks now evaluate by default with cached API fallback when GeneNetwork2 is unavailable
- Vignette `single_cell_integration`: conditional evaluation no longer requires a pre-existing project directory; only Seurat and ggplot2 are needed

## Documentation
- Vignette `locusPackRat_workflow`: added "Building Shareable Packets" section with `buildPacket()` usage and "Including Supplementary CSVs" examples

## Package Metadata
- Exported `buildPacket()` in NAMESPACE
- Added `data-raw` to `.Rbuildignore`

---

# locusPackRat 0.6.1

## Bug Fixes
- **Known Drugs query**: Fixed invalid `mechanismsOfAction { actionType }` field in GraphQL query; now uses `mechanismOfAction` at the KnownDrug row level (matches current Open Targets API schema)
- **Known Drugs query**: Added `size: 500` parameter for complete results (default page size truncated output)
- **Pharmacogenomics query**: Fixed `variantFunctionalConsequence` (now `SequenceOntologyTerm` object) and drug fields (now nested under `drugs: [DrugWithIdentifiers!]!`); rewrote parser for correct extraction
- **Expression parser**: Fixed `lapply` iterating over data.frame columns instead of rows when `fromJSON` auto-simplifies; now produces ~180 tissue rows per gene instead of 3 empty rows
- **Mouse Phenotypes parser**: Applied same data.frame-to-list-of-rows guard to prevent column iteration bug; renamed `targetFromSourceId` to `targetInModelEnsemblId` to match current schema
- **DepMap parser**: Fixed column explosion (177 wide-format columns instead of 6) by extracting only expected columns from auto-simplified data.frame
- **QTL credible set gene linking**: Expanded `id_map` to include all ensembl IDs from QTL response (not just project genes), fixing 0/179 match issue in region mode
- **Locus zoom plot path**: `output_file` now accepts only a filename; full paths are reduced to `basename()` with a warning, preventing nested `.locusPackRat/output/<path>/file.pdf` output

## Improvements
- `generateLocusZoomPlot()` supports PNG output via file extension; locus zoom plots display inline in `region_qtl_opentargets` vignette
- Reduced batch size to 5 genes when querying heavy data types (mouse_phenotypes, known_drugs) to avoid API timeouts
- Increased API timeout from 60s to 120s for Open Targets queries
- Improved parser dual-path handling: data.frame vs list-of-lists detection is now consistent across all data type parsers
- Audit script now uses name-based file lookup instead of newest-file heuristic, includes row-count validation and locus zoom plot test
- Vignette `genenetwork_qtl`: defensive `intersect(names(...), expected_cols)` column subsetting
- Vignette `region_qtl_opentargets`: locus zoom chunks conditionally evaluated when Bioconductor annotation packages are available
- Vignette `single_cell_integration`: all chunks conditionally evaluated when Seurat and ggplot2 are installed (fixes R CMD check failure on systems without Seurat)
- Vignette `locusPackRat_workflow`: synced VignetteIndexEntry with YAML title

## Package Metadata
- Added `Seurat`, `TxDb.Hsapiens.UCSC.hg38.knownGene`, `org.Hs.eg.db` to Suggests
- Added `chr_clean`, `study_locus_id`, `study_target_gene_id`, and `qtl_gene_id` to `globalVariables()` (fixes R CMD check NOTEs)
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

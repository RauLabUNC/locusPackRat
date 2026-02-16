<!--
Sync Impact Report
===================
Version change: 0.0.0 (template) → 1.0.0 (initial ratification)

Modified principles: N/A (first ratification — no prior principles existed)

Added sections:
  - Principle I: Project-as-Truth
  - Principle II: Genome Build Fidelity
  - Principle III: External API Resilience
  - Principle IV: Minimal Required Dependencies
  - Principle V: Test-Isolated Reproducibility
  - Principle VI: Stable Public API
  - Dependency Policy
  - R Version Policy
  - Data Format Policy
  - Development Workflow
  - Governance

Removed sections: None

Templates requiring updates:
  - .specify/templates/plan-template.md — ✅ no updates needed
    (Constitution Check section is generic; principles are project-specific gates)
  - .specify/templates/spec-template.md — ✅ no updates needed
    (User stories + acceptance criteria align with Principle V test requirements)
  - .specify/templates/tasks-template.md — ✅ no updates needed
    (Phase structure accommodates R package workflow; test-first pattern matches Principle V)

Follow-up TODOs: None
-->

# locusPackRat Constitution

## Core Principles

### I. Project-as-Truth

All analysis state lives in `.locusPackRat/` with `config.json` as the single
source of truth. Every exported function re-reads from disk via the `project_dir`
path parameter — no in-memory state persists between calls. A user MUST be able
to close R and resume exactly where they left off.

- The `.locusPackRat/` directory is the canonical project root.
- `config.json` stores species, genome build, region coordinates, and all
  user-defined parameters.
- Functions MUST NOT cache state across calls; each invocation reads fresh
  configuration from disk.
- Project portability: moving the `.locusPackRat/` directory to another machine
  MUST preserve full reproducibility.

### II. Genome Build Fidelity

A project is locked to one genome build at `initPackRat()` time.
Species-genome compatibility is validated (e.g., mouse + hg38 = error). All
reference files, coordinates, and orthology mappings MUST match the declared
build. No mixing across builds is permitted.

- Supported builds are declared per species and validated at init.
- Coordinate systems, annotation files, and API queries MUST use the locked
  build exclusively.
- Orthology lookups MUST respect the declared species pair and builds.
- Changing a genome build requires creating a new project.

### III. External API Resilience

All external calls (MouseMine REST, Open Targets GraphQL) MUST implement batch
processing, retry with exponential backoff, rate limiting, and graceful
degradation. Failures emit `warning()`, not `stop()`. Partial results are
preserved.

- Network-dependent functions MUST NOT halt execution on transient failures.
- Batch endpoints MUST be preferred over per-item calls.
- Retry logic MUST use exponential backoff with a configurable maximum.
- Rate limiting MUST respect upstream API constraints.
- Partial results MUST be returned with warnings describing what failed.

### IV. Minimal Required Dependencies

Imports MUST stay small (currently 6 CRAN packages, hard ceiling of 8).
Bioconductor packages stay in `Suggests`, gated by `requireNamespace()`.
Dependencies used in fewer than 2 exported functions with base R alternatives
are candidates for removal.

- The `Imports` field in DESCRIPTION has a ceiling of 8 packages.
- Bioconductor packages MUST NOT appear in `Imports`; they belong in `Suggests`
  and MUST be guarded by `requireNamespace()` with an informative message.
- Before adding a new dependency, verify no base R alternative exists and that
  at least 2 exported functions will use it.
- Active pruning: dependencies that fall below the 2-export threshold after
  refactoring MUST be evaluated for removal.

### V. Test-Isolated Reproducibility

Every test uses `tempdir()`-based projects with `on.exit()` cleanup via
`helper-fixtures.R` helpers. No hardcoded paths, no leftover files, no network
dependency for core tests. Every export needs happy-path and error-path
coverage.

- Test projects MUST be created in `tempdir()` and cleaned up with `on.exit()`.
- Shared test setup lives in `tests/testthat/helper-fixtures.R`.
- Tests MUST NOT depend on network availability for core functionality.
- Network-dependent tests MUST be skipped when offline (`skip_if_offline()`).
- Every exported function MUST have at least one happy-path and one error-path
  test.

### VI. Stable Public API

The exported functions are a stable contract. Parameter names use `snake_case`.
Internal helpers use `.` prefix and `@noRd`. New capabilities are added as
parameters on existing functions before considering new exports. Breaking
changes MUST be documented in NEWS.md.

- Exported functions form the public API; additions require justification.
- Parameter names MUST use `snake_case` consistently.
- Internal helper functions MUST use a `.` prefix (e.g., `.validate_build()`)
  and MUST include `@noRd` to suppress documentation generation.
- Prefer extending existing functions with new parameters over creating new
  exports.
- Any breaking change (renamed parameter, removed function, changed return
  structure) MUST be documented in NEWS.md with migration guidance.

## Dependency Policy

- **Imports ceiling**: 8 CRAN packages maximum in the `Imports` field.
- **Bioconductor isolation**: All Bioconductor packages in `Suggests`, accessed
  via `requireNamespace()` with user-facing install instructions on failure.
- **Active pruning**: After any refactor, audit `Imports` for packages used in
  fewer than 2 exported functions. Replace with base R or move to `Suggests`.
- **New dependency checklist**: (1) no base R alternative, (2) used by >= 2
  exports, (3) actively maintained, (4) does not exceed the ceiling of 8.

## R Version Policy

- Minimum supported version: R >= 4.3.0.
- The minimum version MUST only be raised when the package actively uses a
  language feature or base function introduced in a newer R release.
- Version bumps MUST be documented in NEWS.md with the specific feature
  requiring the newer R version.

## Data Format Policy

- Tabular data stored as CSV files via `data.table::fread()` / `fwrite()`.
- Project configuration stored as JSON (`config.json`).
- Column names in all stored data MUST use `snake_case`.
- No binary or proprietary formats for user-facing data.

## Development Workflow

- `R CMD check` MUST pass with zero errors, zero warnings, and zero notes
  before any merge to `main`.
- Test coverage: every exported function requires happy-path + error-path tests.
- Branching: feature development on `dev/vX.Y.Z` branches, merged to `main`
  via pull request.
- Commit messages: conventional style preferred (`feat:`, `fix:`, `docs:`,
  `test:`, `refactor:`).

## Governance

This constitution supersedes informal conventions and ad-hoc decisions.
All code changes, reviews, and architectural decisions MUST comply with the
principles above. Amendments to this constitution require lead author agreement
and MUST follow semantic versioning:

- **MAJOR**: Removal or incompatible redefinition of a principle.
- **MINOR**: Addition of a new principle or material expansion of an existing
  one.
- **PATCH**: Clarifications, wording improvements, non-semantic refinements.

Amendments MUST update the version, last-amended date, and include a Sync
Impact Report as an HTML comment at the top of this file.

**Version**: 1.0.0 | **Ratified**: 2026-02-16 | **Last Amended**: 2026-02-16

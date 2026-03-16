#' Filter data using relational joins and flexible criteria
#'
#' @description
#' Filter any data frame/table using relational joins with reference tables
#' and/or custom filter criteria. This function provides a flexible way to
#' filter genes (or any data) based on relationships with other tables.
#'
#' @param inputTable Data frame or data.table to filter
#' @param referenceTable Optional reference table to join for filtering
#' @param by Column(s) to join on. Can be:
#'   - A character vector of column names (used for both tables)
#'   - A named vector (names for inputTable, values for referenceTable)
#'   - NULL to auto-detect common column names
#' @param filters Optional list of filter criteria to apply after joining. Each can be:
#'   - A list with 'column', 'condition', and 'value'
#'   - A function that takes the data and returns a logical vector
#'   - A character expression string
#' @param joinType Type of join if referenceTable provided: "inner", "left", "semi", "anti"
#'   - "inner": Keep only rows with matches in reference
#'   - "left": Keep all input rows, add reference columns
#'   - "semi": Keep input rows with matches (no reference columns added)
#'   - "anti": Keep input rows WITHOUT matches
#' @param keepColumns Columns to retain in output (NULL = all columns)
#' @param verbose Print progress and statistics
#'
#' @return Filtered data frame/table
#'
#' @details
#' This function operates in two modes:
#'
#' 1. **Relational filtering**: When a referenceTable is provided, filter based on
#'    relationships between tables (e.g., "keep genes that appear in this list")
#'
#' 2. **Direct filtering**: When only filters are provided, apply them directly
#'    to the input table
#'
#' These modes can be combined - first join with reference, then apply filters.
#'
#' @examples
#' # Example 1: Filter genes by a reference list (e.g., expressed genes)
#' allGenes <- data.frame(
#'   gene_symbol = c("Abc1", "Def2", "Ghi3", "Jkl4"),
#'   expression = c(10, 50, 200, 0)
#' )
#' 
#' expressedGenes <- data.frame(gene_symbol = c("Abc1", "Def2", "Ghi3"))
#' 
#' filtered <- filterGenes(
#'   inputTable = allGenes,
#'   referenceTable = expressedGenes,
#'   by = "gene_symbol",
#'   joinType = "semi"  # Keep only genes in the expressed list
#' )
#'
#' # Example 2: Anti-join to exclude genes
#' excludeList <- data.frame(gene = c("Def2", "Jkl4"))
#' filtered <- filterGenes(
#'   inputTable = allGenes,
#'   referenceTable = excludeList,
#'   by = c("gene_symbol" = "gene"),  # Different column names
#'   joinType = "anti"  # Exclude genes in the list
#' )
#'
#' @export
filterGenes <- function(inputTable,
                       referenceTable = NULL,
                       by = NULL,
                       filters = NULL,
                       joinType = "inner",
                       keepColumns = NULL,
                       verbose = FALSE) {

  # Input validation
  if (!is.data.frame(inputTable)) {
    stop("inputTable must be a data frame or data.table")
  }

  # Convert to data.table for efficient operations
  if (!inherits(inputTable, "data.table")) {
    result <- data.table::as.data.table(inputTable)
  } else {
    result <- data.table::copy(inputTable)
  }

  originalRows <- nrow(result)

  # Step 1: Apply relational filtering if reference table provided
  if (!is.null(referenceTable)) {

    if (!is.data.frame(referenceTable)) {
      stop("referenceTable must be a data frame or data.table")
    }

    if (!inherits(referenceTable, "data.table")) {
      referenceTable <- data.table::as.data.table(referenceTable)
    }

    # Determine join columns
    if (is.null(by)) {
      # Auto-detect common columns
      by <- intersect(names(result), names(referenceTable))
      if (length(by) == 0) {
        stop("No common columns found between tables. Please specify 'by' parameter.")
      }
      if (verbose) {
        cat(sprintf("Auto-detected join columns: %s\n", paste(by, collapse = ", ")))
      }
    }

    # Perform the join based on type
    if (verbose) {
      cat(sprintf("Performing %s join on: %s\n", joinType,
                  if (is.character(by) && is.null(names(by))) paste(by, collapse = ", ")
                  else paste(names(by), "=", by, collapse = ", ")))
    }

    if (joinType == "inner") {
      # Inner join - keep only matching rows with all columns
      result <- merge(result, referenceTable, by = by, all = FALSE)

    } else if (joinType == "left") {
      # Left join - keep all input rows, add reference columns
      result <- merge(result, referenceTable, by = by, all.x = TRUE)

    } else if (joinType == "semi") {
      # Semi join - keep only matching rows, no new columns
      if (is.null(names(by))) {
        # Same column names
        matchKeys <- unique(referenceTable[, by, with = FALSE])
        result <- merge(result, matchKeys, by = by)
      } else {
        # Different column names
        matchKeys <- unique(referenceTable[, by, with = FALSE])
        data.table::setnames(matchKeys, by, names(by))
        result <- merge(result, matchKeys, by = names(by))
      }

    } else if (joinType == "anti") {
      # Anti join - keep only non-matching rows
      if (is.null(names(by))) {
        # Same column names
        matchKeys <- unique(referenceTable[, by, with = FALSE])
        matchKeys$..match.. <- TRUE
        temp <- merge(result, matchKeys, by = by, all.x = TRUE)
        result <- temp[is.na(..match..)][, ..match.. := NULL]
      } else {
        # Different column names
        matchKeys <- unique(referenceTable[, by, with = FALSE])
        data.table::setnames(matchKeys, by, names(by))
        matchKeys$..match.. <- TRUE
        temp <- merge(result, matchKeys, by = names(by), all.x = TRUE)
        result <- temp[is.na(..match..)][, ..match.. := NULL]
      }

    } else {
      stop(sprintf("Invalid joinType '%s'. Must be one of: inner, left, semi, anti", joinType))
    }

    if (verbose) {
      cat(sprintf("  After join: %d rows\n", nrow(result)))
    }
  }

  # Step 2: Apply additional filters if provided
  if (!is.null(filters) && length(filters) > 0) {

    if (verbose) {
      cat("Applying additional filters:\n")
    }

    filterResults <- list()

    for (i in seq_along(filters)) {
      filter <- filters[[i]]

      if (is.function(filter)) {
        # Function filter
        if (verbose) cat(sprintf("  Filter %d: custom function\n", i))
        filterResults[[i]] <- filter(result)

      } else if (is.character(filter) && length(filter) == 1) {
        # Expression string
        if (verbose) cat(sprintf("  Filter %d: %s\n", i, filter))
        filterResults[[i]] <- eval(parse(text = filter), envir = result)

      } else if (is.list(filter)) {
        # Structured filter
        filterResults[[i]] <- applyFilter(result, filter, verbose)

      } else {
        stop(sprintf("Invalid filter type at position %d", i))
      }
    }

    # Combine all filter results with AND logic
    finalFilter <- Reduce("&", filterResults)
    result <- result[finalFilter, ]

    if (verbose) {
      cat(sprintf("  After filters: %d rows\n", nrow(result)))
    }
  }

  # Step 3: Select columns if specified
  if (!is.null(keepColumns)) {
    availableCols <- intersect(keepColumns, names(result))
    if (length(availableCols) < length(keepColumns)) {
      missing <- setdiff(keepColumns, names(result))
      warning(sprintf("Columns not found: %s", paste(missing, collapse = ", ")))
    }
    result <- result[, availableCols, with = FALSE]
  }

  # Report summary
  if (verbose) {
    cat("\n--- Filtering Summary ---\n")
    cat(sprintf("Original rows: %d\n", originalRows))
    cat(sprintf("Final rows: %d\n", nrow(result)))
    cat(sprintf("Removed: %d (%.1f%%)\n",
                originalRows - nrow(result),
                100 * (originalRows - nrow(result)) / originalRows))
  }

  return(result)
}

#' Create a filter specification
#'
#' @description
#' Helper function to create filter specifications for use with filterGenes
#'
#' @param column Column name to filter on
#' @param condition Condition to apply (==, !=, >, <, >=, <=, \%in\%, etc.)
#' @param value Value(s) to compare against
#' @param na.rm Remove NAs when filtering (default: FALSE)
#'
#' @return Filter specification list
#'
#' @examples
#' # Numeric comparison
#' f1 <- makeFilter("expression", ">", 100)
#'
#' # Exact match
#' f2 <- makeFilter("biotype", "==", "protein_coding")
#'
#' # Multiple values
#' f3 <- makeFilter("chromosome", "%in%", c("chr1", "chr2", "chr3"))
#'
#' # Pattern matching
#' f4 <- makeFilter("gene_name", "matches", "^ABC")
#'
#' @export
makeFilter <- function(column, condition, value = NULL, na.rm = FALSE) {

  validConditions <- c("==", "!=", ">", ">=", "<", "<=",
                      "%in%", "%nin%", "matches", "not_matches",
                      "is_na", "not_na", "between")

  if (!condition %in% validConditions) {
    stop(sprintf("Invalid condition '%s'. Valid conditions: %s",
                condition, paste(validConditions, collapse = ", ")))
  }

  # Validate value requirements
  if (condition %in% c("is_na", "not_na")) {
    if (!is.null(value)) {
      warning(sprintf("Condition '%s' doesn't use a value parameter", condition))
    }
    value <- NULL
  } else if (is.null(value)) {
    stop(sprintf("Condition '%s' requires a value parameter", condition))
  }

  if (condition == "between" && length(value) != 2) {
    stop("Condition 'between' requires exactly 2 values")
  }

  list(
    column = column,
    condition = condition,
    value = value,
    na.rm = na.rm
  )
}

#' Apply a single filter to data
#'
#' @description
#' Internal function to apply a structured filter
#'
#' @param data Data to filter
#' @param filter Filter specification
#' @param verbose Whether to print details
#'
#' @return Logical vector
#'
#' @keywords internal
applyFilter <- function(data, filter, verbose = FALSE) {

  column <- filter$column
  condition <- filter$condition
  value <- filter$value
  na.rm <- filter$na.rm %||% FALSE

  if (!column %in% names(data)) {
    warning(sprintf("Column '%s' not found, skipping filter", column))
    return(rep(TRUE, nrow(data)))
  }

  if (verbose) {
    valueStr <- if (is.null(value)) "NULL"
                else if (length(value) > 3) paste0(paste(head(value, 3), collapse = ", "), "...")
                else paste(value, collapse = ", ")
    cat(sprintf("  Filter: %s %s %s\n", column, condition, valueStr))
  }

  # Apply the condition
  result <- switch(condition,
    "==" = data[[column]] == value,
    "!=" = data[[column]] != value,
    ">" = data[[column]] > value,
    ">=" = data[[column]] >= value,
    "<" = data[[column]] < value,
    "<=" = data[[column]] <= value,
    "%in%" = data[[column]] %in% value,
    "%nin%" = !(data[[column]] %in% value),
    "matches" = grepl(value, data[[column]], perl = TRUE),
    "not_matches" = !grepl(value, data[[column]], perl = TRUE),
    "is_na" = is.na(data[[column]]),
    "not_na" = !is.na(data[[column]]),
    "between" = data[[column]] >= value[1] & data[[column]] <= value[2],
    stop(sprintf("Unknown condition: %s", condition))
  )

  # Handle NAs if requested
  if (na.rm && any(is.na(result))) {
    result[is.na(result)] <- FALSE
  }

  result
}

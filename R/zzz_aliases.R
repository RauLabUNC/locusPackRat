# Ergonomic aliases for locusPackRat functions
#
# All lpr_* aliases are simple function assignments — identical behavior
# to the original CamelCase functions. Both names remain exported.

#' @rdname initPackRat
#' @export
lpr_init <- initPackRat

#' @rdname addRatTable
#' @export
lpr_add_table <- addRatTable

#' @rdname removeRatTable
#' @export
lpr_remove_table <- removeRatTable

#' @rdname listPackRatTables
#' @export
lpr_list_tables <- listPackRatTables

#' @rdname makeGeneSheet
#' @export
lpr_export <- makeGeneSheet

#' @rdname buildPacket
#' @export
lpr_packet <- buildPacket

#' @rdname generateLocusZoomPlot
#' @export
lpr_zoom <- generateLocusZoomPlot

#' @rdname queryMouseMine
#' @export
lpr_query_mousemine <- queryMouseMine

#' @rdname queryOpenTargets
#' @export
lpr_query_ot <- queryOpenTargets

#' @rdname queryOpenTargetsQTL
#' @export
lpr_query_ot_qtl <- queryOpenTargetsQTL

#' @rdname filterGenes
#' @export
lpr_filter <- filterGenes

#' @rdname makeFilter
#' @export
lpr_make_filter <- makeFilter

# Reference-code indexing for the classroom examples.

#' Search the supplied upstream source-code catalog
#'
#' Search named R functions, S4 methods and registered cpp11 entry points.
#' The accompanying file catalog also lists native libraries that are not
#' registered directly with R. An index entry is a locator, not validation.
#' @param query Literal search text, matched without case sensitivity.
#' @param package Optional upstream package name.
#' @param table One of symbols, files or packages.
#' @return Data frame containing the matching source records.
#' @export
sourceCatalog <- function(query = "", package = NULL,
                           table = c("symbols", "files", "packages")) {
    table <- match.arg(table)
    path <- system.file("catalog", paste0(table, ".tsv"), package = "hpgR")
    if (!nzchar(path)) stop("The upstream catalog has not been generated", call. = FALSE)
    result <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!is.null(package)) result <- result[result$package %in% package, , drop = FALSE]
    if (!is.character(query) || length(query) != 1L || is.na(query)) {
        stop("query must be one string", call. = FALSE)
    }
    if (nzchar(query)) {
        text <- apply(result, 1L, paste, collapse = " ")
        result <- result[grepl(tolower(query), tolower(text), fixed = TRUE), , drop = FALSE]
    }
    rownames(result) <- NULL
    result
}

# Source inspection for the classroom examples.

#' Inspect the sources used by hpgR
#' @return Data frame of source files, SHA-256 hashes and adaptation scope.
#' @export
sourceProvenance <- function() {
    path <- system.file("provenance", "source-map.tsv", package = "hpgR")
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

#!/usr/bin/env Rscript
# Index supplied source code. This script does not execute upstream algorithms.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript tools/index_upstream.R POPGENTOOLS")
root <- normalizePath(args[1L], mustWork = TRUE)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
destination <- file.path(dirname(normalizePath(script)), "..", "inst", "catalog")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
packages <- list.dirs(root, recursive = FALSE, full.names = TRUE)
packages <- packages[file.exists(file.path(packages, "DESCRIPTION"))]
package_rows <- file_rows <- symbol_rows <- issue_rows <- list()

bodyKey <- function(expr) paste(deparse(expr, width.cutoff = 500L), collapse = "\n")
symbolRow <- function(package, file, kind, symbol, line, end_line, args = "", calls = "") {
    data.frame(package = package, file = file, kind = kind, symbol = symbol,
               line = as.integer(line), end_line = as.integer(end_line),
               arguments = args, calls = calls, stringsAsFactors = FALSE)
}

for (package_dir in packages) {
    desc <- read.dcf(file.path(package_dir, "DESCRIPTION"))
    field <- function(name) if (name %in% colnames(desc)) desc[1L, name] else ""
    pkg <- field("Package")
    package_rows[[length(package_rows) + 1L]] <- data.frame(
        package = pkg, version = field("Version"), title = field("Title"),
        authors = field("Authors@R"), license = field("License"),
        url = field("URL"), depends = field("Depends"), imports = field("Imports"),
        linking_to = field("LinkingTo"), stringsAsFactors = FALSE)
    files <- list.files(package_dir, full.names = TRUE, recursive = TRUE)
    keep <- grepl("[.](R|r|c|cc|cpp|h|hpp|py)$", files) |
        basename(files) %in% c("DESCRIPTION", "NAMESPACE", "COPYING", "LICENSE", "LICENSE.md")
    files <- files[keep]
    hashes <- character(length(files))
    for (start in seq.int(1L, length(files), by = 100L)) {
        idx <- seq.int(start, min(length(files), start + 99L))
        result <- system2("shasum", c("-a", "256", shQuote(files[idx])), stdout = TRUE)
        if (length(result) != length(idx)) stop("Source hash enumeration failed")
        hashes[idx] <- substr(result, 1L, 64L)
    }
    for (i in seq_along(files)) {
        file <- files[i]
        relative <- substring(file, nchar(package_dir) + 2L)
        lines <- readLines(file, warn = FALSE)
        file_rows[[length(file_rows) + 1L]] <- data.frame(
            package = pkg, file = relative, sha256 = hashes[i],
            lines = length(lines), bytes = file.info(file)$size,
            area = if (startsWith(relative, "tests/")) "test" else
                if (startsWith(relative, "R/")) "R" else
                if (startsWith(relative, "src/")) "native" else "other",
            stringsAsFactors = FALSE)
        if (startsWith(relative, "R/") && grepl("[.][Rr]$", relative)) {
            parsed <- tryCatch(parse(file, keep.source = TRUE), error = function(e) e)
            if (inherits(parsed, "error")) {
                issue_rows[[length(issue_rows) + 1L]] <- data.frame(
                    package = pkg, file = relative, issue = conditionMessage(parsed))
                next
            }
            locations <- attr(parsed, "srcref")
            for (j in seq_along(parsed)) {
                expr <- parsed[[j]]
                if (!is.call(expr)) next
                loc <- as.integer(locations[[j]])
                if (identical(expr[[1L]], as.name("<-")) && is.symbol(expr[[2L]]) &&
                    is.call(expr[[3L]]) && identical(expr[[3L]][[1L]], as.name("function"))) {
                    fun <- eval(expr[[3L]], envir = baseenv())
                    calls <- tryCatch(codetools::findGlobals(fun, merge = FALSE)$functions,
                                      error = function(e) character())
                    symbol_rows[[length(symbol_rows) + 1L]] <- symbolRow(
                        pkg, relative, "R_function", as.character(expr[[2L]]), loc[1L], loc[3L],
                        paste(names(formals(fun)), collapse = ";"), paste(calls, collapse = ";"))
                } else if (identical(expr[[1L]], as.name("setMethod"))) {
                    symbol <- paste0(bodyKey(expr[[2L]]), " [", bodyKey(expr[[3L]]), "]")
                    symbol_rows[[length(symbol_rows) + 1L]] <- symbolRow(
                        pkg, relative, "R_S4_method", symbol, loc[1L], loc[3L])
                }
            }
        } else if (grepl("[.](c|cc|cpp|h|hpp)$", relative)) {
            for (j in grep("[[cpp11::register]]", lines, fixed = TRUE)) {
                fragment <- paste(sub("//.*$", "", lines[seq.int(j + 1L,
                    min(length(lines), j + 15L))]), collapse = " ")
                hit <- regmatches(fragment, regexpr("[A-Za-z_][A-Za-z_0-9]*[[:space:]]*[(]",
                                                    fragment))
                if (length(hit) && nzchar(hit)) {
                    name <- sub("[[:space:]]*[(]$", "", hit)
                    symbol_rows[[length(symbol_rows) + 1L]] <- symbolRow(
                        pkg, relative, "cpp11_entry", name, j, j)
                }
            }
        }
    }
}
clean <- function(x) {
    x[] <- lapply(x, function(v) if (is.character(v)) gsub("[\t\r\n]+", " ", v) else v)
    x
}
writeCatalog <- function(rows, filename) {
    table <- clean(do.call(rbind, rows))
    write.table(table, file.path(destination, filename), sep = "\t", quote = TRUE,
                row.names = FALSE, na = "")
    table
}
pkg <- writeCatalog(package_rows, "packages.tsv")
files <- writeCatalog(file_rows, "files.tsv")
symbols <- writeCatalog(symbol_rows, "symbols.tsv")
if (length(issue_rows)) writeCatalog(issue_rows, "issues.tsv")
cat("Indexed", nrow(pkg), "packages,", nrow(files), "source files and",
    nrow(symbols), "R functions/S4 methods/cpp11 entries.\n")
cat("R parse failures:", length(issue_rows), "\n")
print(with(symbols, table(package, kind)))

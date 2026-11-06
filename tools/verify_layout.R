#!/usr/bin/env Rscript
# Check source organization without executing scientific algorithms.
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
root <- if (length(args)) normalizePath(args[1L], mustWork = TRUE) else
    normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
manifest <- utils::read.delim(file.path(root, "inst/provenance/modules.tsv"),
                            stringsAsFactors = FALSE, check.names = FALSE)
sources <- utils::read.delim(file.path(root, "inst/provenance/source-map.tsv"),
                           stringsAsFactors = FALSE, check.names = FALSE)
files <- list.files(file.path(root, "R"), pattern = "[.]R$")
stopifnot(!anyDuplicated(manifest$file), setequal(files, manifest$file))
stopifnot(all(manifest$kind %in% c("upstream", "pedagogical", "utility", "documentation")))
approved_sources <- c("haplotypeTools", "ancestryTools", "demographyTools",
                     "ARGtools", "tskitR", "simuPOP")
stopifnot(!anyDuplicated(sources$id), all(sources$upstream %in% approved_sources),
          all(grepl("^[[:xdigit:]]{64}$", sources$sha256)))
functions <- character()
pedagogical_functions <- character()
for (i in seq_len(nrow(manifest))) {
    file <- manifest$file[i]
    kind <- manifest$kind[i]
    if (kind == "documentation") {
        stopifnot(identical(file, "hpgR-package.R"))
    } else {
        stopifnot(grepl("^[a-z][a-z0-9_]*[.]R$", file))
        stopifnot(identical(startsWith(file, "pedagogical_"), kind == "pedagogical"))
        if (kind == "utility") stopifnot(identical(file, "input_validation.R"))
        stopifnot(!grepl("^(ported|adopted|adapted|adapter|simple|utils)_", file))
    }
    if (kind == "upstream") {
        ids <- strsplit(manifest$source_ids[i], ";", fixed = TRUE)[[1L]]
        stopifnot(length(ids) > 0L, all(nzchar(ids)), all(ids %in% sources$id))
    } else {
        stopifnot(is.na(manifest$source_ids[i]) || !nzchar(manifest$source_ids[i]))
    }
    for (expr in parse(file.path(root, "R", file))) {
        if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
            is.symbol(expr[[2L]]) && is.call(expr[[3L]]) &&
            identical(expr[[3L]][[1L]], as.name("function"))) {
            functions <- c(functions, as.character(expr[[2L]]))
            if (kind == "pedagogical") pedagogical_functions <-
                c(pedagogical_functions, as.character(expr[[2L]]))
        }
    }
}
simupop_ids <- sources$id[sources$upstream == "simuPOP"]
simupop_files <- vapply(strsplit(manifest$source_ids, ";", fixed = TRUE),
    function(ids) any(ids %in% simupop_ids), logical(1))
stopifnot(all(grepl("^simupop(_[a-z0-9]+)*[.]R$", manifest$file[simupop_files])),
          !dir.exists(file.path(root, "docs")),
          file.exists(file.path(root, "inst", "PORTING_SOP.md")))
stopifnot(!anyDuplicated(functions),
          all(grepl("^[a-z][A-Za-z0-9]*$", functions)),
          !any(grepl("lecture|example[_.]?[0-9]", c(files, functions), ignore.case = TRUE)))
review <- utils::read.delim(file.path(root,"inst/catalog/teaching-functions.tsv"))
stopifnot(!anyDuplicated(review$function_name),
          setequal(review$function_name,pedagogical_functions),
          all(review$file %in% manifest$file[manifest$kind=="pedagogical"]),
          all(nzchar(review$decision)),all(nzchar(review$candidate)))
cat("PASS:", length(files), "registered R files;", length(functions),
    "unique function definitions;", nrow(review), "individual pedagogical reviews.\n")

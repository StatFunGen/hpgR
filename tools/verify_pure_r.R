#!/usr/bin/env Rscript
# Source/payload gate and independent R probability-model checks.
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
root <- if (length(args)) normalizePath(args[1L], mustWork = TRUE) else
    if (length(script)) normalizePath(file.path(dirname(script), ".."), mustWork = TRUE) else
    normalizePath(getwd(), mustWork = TRUE)
library(hpgR)
checks <- 0L
assert <- function(ok, message) {
    if (!isTRUE(ok)) stop(message, call. = FALSE)
    checks <<- checks + 1L
}

native <- c(".Call", ".External", ".External2", ".C", ".Fortran",
    "dyn.load", "library.dynam", "system", "system2", "shell", "pipe")
external <- c("reticulate", "Rcpp", "rPython", "basilisk")
violations <- character()
visit <- function(x, file) {
    if (!is.call(x)) {
        if (is.character(x) && length(x) &&
            any(grepl("(^|\n)[[:space:]]*(import[[:space:]]+[[:alnum:]_]+|from[[:space:]]+[[:alnum:]_.]+[[:space:]]+import|def[[:space:]]+[[:alnum:]_]+[[:space:]]*\\()", x)))
            violations <<- c(violations, paste(basename(file), "embedded Python source", sep = ": "))
        if (is.pairlist(x) || is.expression(x) || is.list(x))
            invisible(lapply(as.list(x), visit, file = file))
        return(invisible(NULL))
    }
    h <- x[[1L]]
    bad <- is.symbol(h) && as.character(h) %in% native
    if (is.symbol(h) && as.character(h) %in% c("::", ":::") &&
        length(x) >= 3L && (as.character(x[[2L]]) %in% external ||
            as.character(x[[3L]]) %in% native)) bad <- TRUE
    if (is.symbol(h) && as.character(h) %in% c("library", "require",
        "requireNamespace", "getNamespace", "asNamespace", "getExportedValue") &&
        length(x) >= 2L && as.character(x[[2L]])[1L] %in% external) bad <- TRUE
    if (is.symbol(h) && as.character(h) == "do.call" && length(x) >= 2L &&
        is.character(x[[2L]]) && x[[2L]] %in% native) bad <- TRUE
    if (bad) violations <<- c(violations,
        paste(basename(file), paste(deparse(x, width.cutoff = 120L), collapse = " "), sep = ": "))
    invisible(lapply(as.list(x), visit, file = file))
}
for (file in list.files(file.path(root, "R"), "[.]R$", full.names = TRUE))
    invisible(lapply(as.list(parse(file)), visit, file = file))
assert(!length(violations), paste(c("R source calls an external runtime:", violations), collapse = "\n"))

description <- read.dcf(file.path(root, "DESCRIPTION"))
deps <- paste(description[, intersect(c("Imports", "Depends", "LinkingTo", "Suggests"),
    colnames(description))], collapse = ",")
assert(!grepl("(^|[,[:space:]])reticulate([[:space:](,]|$)", deps),
    "DESCRIPTION still declares reticulate")
files <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE)
ignore_file <- file.path(root, ".Rbuildignore")
if (file.exists(ignore_file)) {
    patterns <- readLines(ignore_file, warn = FALSE)
    patterns <- patterns[nzchar(patterns)]
    files <- files[!vapply(files, function(file) {
        parts <- strsplit(file, "/", fixed = TRUE)[[1L]]
        ancestors <- vapply(seq_along(parts), function(i)
            paste(parts[seq_len(i)], collapse = "/"), character(1))
        any(vapply(patterns, function(pattern)
            any(grepl(pattern, ancestors, ignore.case = TRUE)), logical(1)))
    }, logical(1))]
}
assert(!any(grepl("[.]py[cow]?$", files, ignore.case = TRUE)),
    paste("Python files are still shipped:", paste(files[grepl("[.]py[cow]?$", files,
        ignore.case = TRUE)], collapse = ", ")))
assert(!any(grepl("[.](so|dll|dylib|pyd)$", files, ignore.case = TRUE)),
    "A native library is shipped with the supposed pure-R package")

implemented <- c("inferAncestry", "dateAncestry")
assert(all(implemented %in% getNamespaceExports("hpgR")),
    "Sequence inference and node dating must be exported")
manifest <- sourceProvenance()
for (name in implemented) {
    rows <- which(vapply(strsplit(manifest$functions, ";", fixed = TRUE),
        function(functions) name %in% functions, logical(1)))
    text <- if (length(rows)) paste(unlist(manifest[rows,
        intersect(c("status", "adaptation", "validation_scope"), names(manifest)),
        drop = FALSE]), collapse = " ") else ""
    assert(length(rows) && grepl("r (translation|implementation|adaptation)",
        tolower(text)), paste(name, "lacks an implemented source record"))
}

Ne <- 37
waits <- vapply(seq_len(400L), function(i) {
    x <- simulateAncestry(1L, Ne, 10, 0, seed = 19500L + i)
    max(x$nodes$time)
}, numeric(1))
assert(all(waits > 0), "A neutral pair did not have a positive ancestor time")
assert(abs(mean(waits) / (2 * Ne) - 1) < .20,
    "The neutral pair waiting-time mean differs from 2Ne")
assert(abs(mean(waits <= 2 * Ne) - (1 - exp(-1))) < .12,
    "The neutral pair waiting-time distribution differs from its exponential law")

tree <- simulateAncestry(1L, 50, 100, 0, seed = 19601L)
age <- setNames(tree$nodes$time, tree$nodes$id)
exposure <- sum((tree$edges$right - tree$edges$left) *
    (age[as.character(tree$edges$parent)] - age[as.character(tree$edges$child)]))
mu <- 20 / exposure
counts <- vapply(seq_len(300L), function(i) {
    x <- ancestryData(mutateAncestry(tree, mu, seed = 19700L + i))
    assert(nrow(x$observed$haplotypes) == 2L &&
        all(colSums(x$observed$haplotypes) == 1L),
        "Fixed two-copy tree mutations do not mark a single sampled branch")
    assert(!anyDuplicated(x$observed$positions) &&
        all(x$observed$positions >= 0 & x$observed$positions < 100),
        "Continuous mutation positions violate infinite-sites support")
    as.numeric(x$nSites)
}, numeric(1))
assert(abs(mean(counts) - 20) < 1.25, "The conditional Poisson mean was not recovered")
assert(abs(stats::var(counts) - 20) < 6, "The conditional Poisson variance was not recovered")

fixture <- structure(list(nodes = data.frame(id = 0:6,
    time = c(0, 0, 0, 2, 5, 3, 7), flags = c(1L, 1L, 1L, rep(0L, 4L))),
    edges = data.frame(left = rep(c(0, 5), each = 4L),
        right = rep(c(5, 10), each = 4L),
        parent = c(3L, 3L, 4L, 4L, 5L, 5L, 6L, 6L),
        child = c(0L, 1L, 3L, 2L, 0L, 2L, 5L, 1L)),
    sampleIds = 0:2, sequenceLength = 10, timeUnits = "generations"),
    class = c("hpgR_ancestry", "list"))
placed <- hpgR:::.msprimeMutations(fixture, 0, 19901L,
    events = data.frame(edge = c(3L, 7L, 8L), position = c(1, 5, 9),
        time = c(3, 5, 2)))
extracted <- ancestryData(placed)
assert(identical(extracted$observed$haplotypes,
    rbind(c(1L, 1L, 0L), c(1L, 0L, 1L), c(0L, 1L, 0L))),
    "Hand-checked interval mutation carriers differ")
assert(identical(extracted$observed$positions, c(1, 5, 9)),
    "A mutation at the exact interval boundary was not retained correctly")
assert(isTRUE(all.equal(extracted$intervals$pair_time, c(2, 7))) &&
    isTRUE(all.equal(extracted$intervals$tmrca, c(5, 7))),
    "Hand-checked pair/root ages differ between the two local trees")

cat("PASS:", checks, "pure-R boundary and probability-law checks.\n")
cat(sprintf("Pair mean %.5f, expected %.5f; mutation mean %.5f, variance %.5f, expected20.\n",
    mean(waits), 2 * Ne, mean(counts), stats::var(counts)))
cat("Sequence inference and inside-outside node dating use source-mapped R functions.\n")

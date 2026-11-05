#!/usr/bin/env Rscript
# Optional: --truth-file PRIVATE_P1_TRUTH_RDS validates all nine existing
# no-recombination P1 regions without regenerating their full ancestry records.
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
lecture <- Sys.getenv("HPGR_COURSE_DIR")
library(hpgR)
library(hpgR)

failure <- function(expr) stopifnot(inherits(try(expr, silent = TRUE), "try-error"))
compareRecord <- function(record) {
    stopifnot(record$map_length == 0, length(record$local_trees) == 1L)
    n <- length(record$sample_nodes)
    tree <- record$local_trees[[1L]]
    quick <- pgForwardExposure(record$size_history$N, n,
                                 record$ancestry_seed,
                                 sample_indices = record$sample_indices)
    reference <- tabulate(lengths(tree$descendants), nbins = n - 1L)
    stopifnot(identical(as.numeric(reference), quick$branch_generations),
              isTRUE(all.equal(as.numeric(tree$root_age), as.numeric(quick$root_age))),
              identical(record$extension_generations, quick$extension_generations))
    invisible(quick)
}

# Tiny populations exercise simultaneous mergers, changing pool sizes, early
# coalescence, an empty forward horizon, and nondefault preselected samples.
designs <- list(list(sizes = rep(1L, 3L), n = 2L, sample = c(1L, 2L)),
                list(sizes = c(3L, 2L, 4L, 2L), n = 4L, sample = 1:4),
                list(sizes = rep(4L, 8L), n = 4L, sample = c(2L, 4L, 6L, 8L)),
                list(sizes = 3L, n = 3L, sample = c(2L, 4L, 6L)))
for (i in seq_along(designs)) {
    d <- designs[[i]]
    record <- pgRecordAncestry(d$sizes, 1000, d$n, ancestry_seed = 81800L + i,
                                 sample_indices = d$sample)
    compareRecord(record)
}
cat("PASS: exact full-map exposure, root age and extension on four tiny designs\n")

truth_index <- match("--truth-file", args)
if (!is.na(truth_index)) {
    stopifnot(truth_index < length(args), file.exists(args[truth_index + 1L]))
    truth <- readRDS(args[truth_index + 1L])
    compared <- 0L
    for (dataset in truth) for (region in dataset$regions) {
        record <- region$ancestry
        if (record$map_length != 0) next
        compareRecord(record)
        compared <- compared + 1L
    }
    stopifnot(compared == 9L)
    cat("PASS: exact exposure/root/extension against all nine stored P1 regions\n")
}

# A separately initialized mutation RNG draws independent aggregate Poisson bins.
exposure <- pgForwardExposure(rep(8L, 5L), n = 6L, ancestry_seed = 9912L)
expected <- 2e-4 * 1000 * exposure$branch_generations
mutations <- pgWithSeed(9941L, function()
    replicate(20000L, rpois(length(expected), expected)))
stopifnot(all(abs(rowMeans(mutations) - expected) <
                  6 * sqrt(pmax(expected, 1e-12) / ncol(mutations)) + 1e-10),
          all(abs(apply(mutations, 1, var) - expected) <
                  .1 * pmax(expected, 1)))
zero <- pgForwardSfs(rep(3L, 4L), n = 4L, mu = 0)
stopifnot(length(zero$sfs) == 3L, all(zero$sfs == 0),
          identical(pgForwardSfs(rep(3L, 4L), n = 4L),
                    pgForwardSfs(rep(3L, 4L), n = 4L)))
failure(pgForwardExposure(0L, 2L))
failure(pgForwardExposure(2L, c(2L, 3L)))
failure(pgForwardExposure(2L, 5L))
failure(pgForwardExposure(2L, 2.5))
failure(pgForwardExposure(2L, 2L, sample_indices = c(1L, 1L)))
failure(pgForwardExposure(2L, 2L, sample_indices = c(1L, 5L)))
failure(pgForwardExposure(2L, 2L, max_extension = -1L))
failure(pgForwardExposure(2L, 2L, max_extension = 0L))
failure(pgForwardSfs(3L, n = 2L, mu = -1))
failure(pgForwardSfs(3L, n = 2L, length_bp = 0))
set.seed(426L)
before <- .Random.seed
invisible(pgForwardSfs(rep(3L, 4L), n = 4L))
stopifnot(identical(before, .Random.seed))
rm(".Random.seed", envir = .GlobalEnv)
invisible(pgForwardSfs(rep(3L, 4L), n = 4L))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
cat("PASS: Poisson means/variances, zero mutation, seed replay, invalid input, RNG restoration\n")

elapsed <- system.time({
    library <- lapply(seq_len(1000L), function(i)
        pgForwardSfs(c(80L, rep(20L, 40L), rep(80L, 20L)),
                        ancestry_seed = 1500000L + i,
                        mutation_seed = 1600000L + i))
})[["elapsed"]]
stopifnot(length(library) == 1000L,
          all(vapply(library, function(x) length(x$sfs) == 11L, logical(1))))
cat(sprintf("PASS: 1000 full-parent-map regions in %.3f seconds\n", elapsed))

#!/usr/bin/env Rscript
# Read the original supplied sources; do not install or modify their packages.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript tools/verify_upstream.R POPGENTOOLS")
root <- normalizePath(args[1L], mustWork = TRUE)
library(hpgR)
checks <- 0L
check <- function(ok, label) {
    if (!isTRUE(ok)) stop("FAIL: ", label)
    checks <<- checks + 1L
}

loadFunctions <- function(path, wanted) {
    env <- new.env(parent = baseenv())
    for (expr in parse(path)) {
        if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
            is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% wanted) {
            eval(expr, env)
        }
    }
    if (!all(wanted %in% ls(env, all.names = TRUE))) stop("Missing reference function")
    env
}

manifest <- sourceProvenance()
for (i in which(manifest$upstream != "simuPOP")) {
    path <- file.path(root, manifest$upstream[i], manifest$file[i])
    hash <- strsplit(system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE),
                     "[[:space:]]+")[[1L]][1L]
    check(identical(hash, manifest$sha256[i]), paste("source hash", manifest$id[i]))
}

wanted <- c(".w_matrix", ".expm1d", ".etjj_constant", ".truncated_sfs")
packageNames <- c("wMatrix", "expm1d", "etjjConstant", "truncatedSfs")
renameCalls <- function(x) {
    if (is.symbol(x)) {
        name <- as.character(x)
        return(if (name %in% wanted) as.name(packageNames[match(name, wanted)]) else x)
    }
    if (is.call(x) || is.expression(x) || is.pairlist(x)) {
        for (i in seq_along(x)) {
            if (is.null(x[[i]]) || identical(x[[i]], quote(expr = ))) next
            x[[i]] <- renameCalls(x[[i]])
        }
    }
    x
}
sfs <- loadFunctions(file.path(root, "demographyTools/R/demography_sfs.R"), wanted)
for (nm in wanted) {
    name <- packageNames[match(nm, wanted)]
    check(identical(deparse(body(get(name, asNamespace("hpgR"))), width.cutoff = 500L),
                    deparse(renameCalls(body(sfs[[nm]])), width.cutoff = 500L)),
          paste("identical R calculation after function renaming", nm))
}
for (n in c(2L, 3L, 4L, 8L, 12L, 20L)) {
    reference <- sfs$.truncated_sfs(n, 20000, Inf)[seq.int(2L, n)] * 1e-8 * 10000
    # Keep the public interface's multiplication order in both calculations.
    reference <- (1e-8 * 10000) * sfs$.truncated_sfs(n, 20000, Inf)[seq.int(2L, n)]
    check(identical(expectedSFS(n, 10000, 1e-8, 10000), reference),
          paste("bit-identical SFS", n))
}

post <- loadFunctions(file.path(root, "demographyTools/R/bayesian_utils.R"),
                      ".posterior_summary")
check(identical(body(posteriorSummary), body(post$.posterior_summary)),
      "identical posterior-summary R body")
set.seed(8149)
for (p in c(1L, 3L, 8L)) {
    x <- matrix(rnorm(100L * p), ncol = p)
    for (f in list(NULL, exp)) {
        check(identical(posteriorSummary(x, f), post$.posterior_summary(x, f)),
              paste("bit-identical posterior summary", p, is.null(f)))
    }
}

# Evaluate the original method functions on a minimal S4 reference object.
# The reference method bodies are untouched; no Bioconductor classes are needed.
ref <- new.env(parent = baseenv())
methods::setClass("hpgOracleRateMap", slots = c(position = "numeric", rate = "numeric"))
ref$RateMap <- function(position, rate)
    methods::new("hpgOracleRateMap", position = as.numeric(position), rate = as.numeric(rate))
for (expr in parse(file.path(root, "tskitR/R/msprime-rate-map.R"))) {
    if (is.call(expr) && identical(expr[[1L]], as.name("setMethod")) &&
        expr[[2L]] %in% c("sequenceLength", "getRate", "getCumulativeMass", "meanRate")) {
        ref[[expr[[2L]]]] <- eval(expr[[4L]], ref)
    }
}
for (rates in list(c(.1, .2), c(0, .125), c(1e-8, 2e-8))) {
    actual <- rateMap(c(0, 10, 30), rates)
    original <- ref$RateMap(c(0, 10, 30), rates)
    pos <- c(0, 1, 9.5, 10, 20, 30)
    check(identical(rateAt(actual, pos), ref$getRate(original, pos)), "bit-identical map rates")
    check(identical(cumulativeRate(actual, pos), ref$getCumulativeMass(original, pos)),
          "bit-identical cumulative map")
    check(identical(meanRate(actual), ref$meanRate(original)), "bit-identical mean rate")
}
cat("PASS:", checks, "direct-source R checks; no full upstream package installed.\n")

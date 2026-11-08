#!/usr/bin/env Rscript
# Verify long matches using founder-copy labels replayed through the same meioses.
# Rscript tools/verify_population_ibd.R [OBSERVED_RDS]
main <- function() {
    library(hpgR)
    args <- commandArgs(trailingOnly = TRUE)
    script <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())])
    root <- normalizePath(file.path(dirname(script), '..'), mustWork = TRUE)
    input <- if (length(args)) args[1L] else
        file.path(root, 'notebooks/data/population_history/observed.rds')
    expected <- readRDS(input)
    specification <- populationHistoryDesign()
    design <- specification$design
    candidate <- specification$candidates$recent_founding
    state <- new.env(parent = emptyenv())
    state$labels <- NULL
    state$sample <- NULL
    state$busy <- FALSE
    state$founderSize <- candidate$founder
    state$ancestralRows <- 2L * design$ancestral_size
    traceName <- '.hpgRIbdValidation'
    if (exists(traceName, .GlobalEnv, inherits = FALSE)) stop('Validation state already exists')
    assign(traceName, state, .GlobalEnv)
    ns <- asNamespace('hpgR')
    trace('simupopSample', where = ns, print = FALSE, exit = quote({
        state <- get('.hpgRIbdValidation', .GlobalEnv)
        value <- returnValue()
        if (is.null(state$labels) && size == state$founderSize &&
            nrow(haplotypes) == state$ancestralRows) {
            state$labels <- matrix(rep(seq_len(nrow(value$haplotypes)),
                ncol(value$haplotypes)), nrow = nrow(value$haplotypes))
        } else if (!is.null(state$labels)) {
            rows <- as.vector(rbind(2L * value$indices - 1L, 2L * value$indices))
            state$sample <- state$labels[rows, , drop = FALSE]
        }
    }))
    trace('simupopMeiosis', where = ns, print = FALSE, exit = quote({
        state <- get('.hpgRIbdValidation', .GlobalEnv)
        if (!state$busy && !is.null(state$labels)) {
            state$busy <- TRUE
            state$labels <- hpgR:::simupopMeiosis(state$labels, parents, loci, fractions,
                events = list(starts = starts, switches = switches))
            state$busy <- FALSE
        }
    }))
    on.exit({
        untrace('simupopSample', where = ns)
        untrace('simupopMeiosis', where = ns)
        rm(list = traceName, envir = .GlobalEnv)
    }, add = TRUE)
    simulation <- simulatePopulationHistory(candidate, design, design$observation_seed,
                                             recordEvery = 1L)
    stopifnot(identical(simulation$observed, expected))
    rows <- seq.int(1L, nrow(expected$haplotypes), 2L)
    h <- expected$haplotypes[rows, , drop = FALSE]
    labels <- state$sample[rows, , drop = FALSE]
    p <- colMeans(h)
    use <- which(p >= design$ibd_minimum_maf & p <= 1 - design$ibd_minimum_maf)
    m <- expected$markers[use, ]
    h <- h[, use, drop = FALSE]
    segments <- detectIbdSegments(h, m$genetic_position, m$chromosome,
        minimumCm = design$ibd_length_breaks[1L])
    agreement <- vapply(seq_len(nrow(segments)), function(k) {
        index <- seq.int(use[segments$fromSite[k]], use[segments$toSite[k]])
        mean(labels[segments$hap1[k], index] == labels[segments$hap2[k], index])
    }, numeric(1))
    set.seed(937713L)
    shuffled <- apply(h, 2L, sample)
    stopifnot(identical(colSums(h), colSums(shuffled)))
    null <- detectIbdSegments(shuffled, m$genetic_position, m$chromosome,
                              minimumCm = design$ibd_length_breaks[1L])
    precision <- weighted.mean(agreement, segments$lengthCm)
    stopifnot(precision > .95, nrow(null) == 0L)
    print(data.frame(matches = nrow(segments),
        lengthWeightedFounderAgreement = precision,
        fractionSegmentsAbove90Agreement = mean(agreement >= .9),
        shuffledMatches = nrow(null)), row.names = FALSE)
    cat('PASS: observed DNA reproduced; founder labels and marker-wise control checked.\n')
}
main()

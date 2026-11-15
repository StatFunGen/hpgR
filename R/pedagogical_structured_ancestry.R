#' Coalescence of a pair sampled across two populations
#'
#' Each surviving lineage independently switches parental populations.
#' Lineages in the same population coalesce with probability 1 divided by
#' its number of gene copies. After the split time, both enter the ancestor.
#' @param loci Number of independent regions.
#' @param migration Vector of symmetric migration probabilities in successive
#'   backward generations before the split. Its length is the split time.
#' @param copies,ancestralCopies Number of gene copies per daughter population
#'   and in the common ancestor, respectively.
#' @return One time to coalescence in generations per region.
#' @examples
#' simulateCrossPopulationTimes(5, rep(0.02, 50), 200, 600)
#' @export
simulateCrossPopulationTimes <- function(loci, migration, copies, ancestralCopies) {
    checkInputs(loci = scalarInteger,
        migration = list(numericVector, minimum = 0, maximum = 1, empty = TRUE),
        copies = scalarInteger, ancestralCopies = scalarInteger)
    together <- rep(FALSE, loci)
    times <- rep(NA_integer_, loci)
    for (g in seq_along(migration)) {
        live <- which(is.na(times))
        switch <- runif(length(live)) < 2 * migration[g] * (1 - migration[g])
        together[live] <- xor(together[live], switch)
        eligible <- live[together[live]]
        times[eligible[runif(length(eligible)) < 1 / copies]] <- g
    }
    live <- is.na(times)
    times[live] <- length(migration) + 1L + stats::rgeom(sum(live), 1 / ancestralCopies)
    times
}

#' ABBA and BABA counts under a split and a pulse
#'
#' Trace one copy per population at independent regions in the tree
#' ((P1, P2), P3). A fraction of the P2 lineages enter P3 at the pulse.
#' Mutations on the internal branch are Poisson conditional on its length.
#' @param loci Number of independent regions.
#' @param ancestry Fraction of P2 ancestry from the P3 pulse.
#' @param Ne Common diploid effective size.
#' @param split12,split123,pulseTime Event times in generations before sampling,
#'   with the pulse no older than the P1--P2 split.
#' @param mutationRate Mutation rate per region per branch-generation.
#' @return Named ABBA and BABA mutation counts.
#' @examples
#' simulateAbbaBaba(1000, .2, 1000, 1000, 1600, 300, .0002)
#' @export
simulateAbbaBaba <- function(loci, ancestry, Ne, split12, split123, pulseTime,
                             mutationRate) {
    checkInputs(loci = scalarInteger, ancestry = probability, Ne = positiveScalar,
        split12 = positiveScalar, split123 = positiveScalar,
        pulseTime = list(positiveScalar, zero = TRUE),
        mutationRate = list(positiveScalar, zero = TRUE))
    if (pulseTime > split12 || split12 > split123)
        stop("event times must satisfy pulseTime <= split12 <= split123", call. = FALSE)
    moved <- runif(loci) < ancestry
    window <- split123 - ifelse(moved, pulseTime, split12)
    wait <- rexp(loci, rate = 1 / (2 * Ne))
    early <- wait < window
    pair <- sample(c("12", "23", "13"), loci, replace = TRUE)
    pair[early] <- ifelse(moved[early], "23", "12")
    length <- rexp(loci, rate = 1 / (2 * Ne)) + ifelse(early, window - wait, 0)
    mutations <- rpois(loci, mutationRate * length)
    c(ABBA = sum(mutations[pair == "23"]), BABA = sum(mutations[pair == "13"]))
}

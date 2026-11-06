# From simuPOP: RandomSelection, CloneGenoTransmitter::applyDuringMating and
# ParentsTagger::applyDuringMating.

simupopDrawParents <- function(parentCopies, offspringCopies) {
    sample.int(parentCopies, offspringCopies, replace = TRUE)
}

simupopParentMaps <- function(sizes) {
    lapply(seq_len(length(sizes) - 1L), function(g)
        simupopDrawParents(2 * sizes[g], 2 * sizes[g + 1L]))
}

#' Copy recorded parental genotypes
#' @param genotypes Matrix with haploid individuals in rows and loci in columns.
#' @param parents One-based parental row index for each offspring.
#' @return Integer matrix of offspring genotypes.
#' @export
copyGenotypes <- function(genotypes, parents) {
    genotypes <- genotypeMatrix(genotypes)
    if (!is.numeric(parents) || !length(parents) || any(!is.finite(parents)) ||
        any(parents < 1 | parents > nrow(genotypes) | parents != floor(parents))) {
        stop("parents must index parental genotype rows", call. = FALSE)
    }
    unname(genotypes[as.integer(parents), , drop = FALSE])
}

#' Replay a recorded copying history
#' @param initial Initial haploid genotype matrix.
#' @param parents List of one-based parental-index vectors, one per generation.
#' @return List of integer genotype matrices, including the initial state.
#' @export
replayCopying <- function(initial, parents) {
    initial <- unname(genotypeMatrix(initial))
    if (!is.list(parents)) stop("parents must be a list", call. = FALSE)
    states <- vector("list", length(parents) + 1L)
    states[[1L]] <- initial
    for (g in seq_along(parents)) {
        states[[g + 1L]] <- copyGenotypes(states[[g]], parents[[g]])
    }
    states
}

#' Generate copying records using simuPOP's neutral copying model
#'
#' R translation of RandomSelection, CloneGenoTransmitter and ParentsTagger.
#' Each row is a haploid individual; a diploid Wright-Fisher copy-pool example
#' supplies 2N rows. This model has no mutation, recombination or selection.
#' @param initial Initial haploid genotype matrix.
#' @param sizes Diploid copy-pool sizes including generation zero.
#' @param seed Positive integer R seed. The caller's RNG state is restored.
#' @return List with complete states, one-based parent indices, sizes and seed.
#' @examples
#' initial <- matrix(c(rep(1L, 3), rep(0L, 7)), ncol = 1L)
#' record <- simulateForward(initial, sizes = c(5, 5, 3, 3, 5), seed = 8149L)
#' states <- replayCopying(initial, record$parents)
#' stopifnot(identical(states, record$states))
#' vapply(states, function(G) mean(G[, 1]), numeric(1L))
#' @export
simulateForward <- function(initial, sizes, seed = 8149L) {
    initial <- unname(genotypeMatrix(initial))
    if (!is.numeric(sizes) || !length(sizes) || any(!is.finite(sizes)) ||
        any(sizes < 1 | sizes != floor(sizes)) ||
        any(sizes > .Machine$integer.max / 2) || nrow(initial) != 2 * sizes[1L]) {
        stop("sizes must give diploid copy-pool sizes; initial must have 2*sizes[1] rows",
             call. = FALSE)
    }
    checkInputs(seed = scalarInteger)
    seed <- as.integer(seed)
    simupopWithSeed(seed, function() {
        parents <- simupopParentMaps(sizes)
        states <- replayCopying(initial, parents)
        list(states = states, parents = parents, sizes = as.integer(sizes),
             seed = seed,
             rng = "R Mersenne-Twister",
             model = "R translation of simuPOP haploid RandomSelection; 2N copy pool; no mutation")
    })
}

#' Simulate a diploid population with translated simuPOP operators
#'
#' Each gamete chooses a parental diploid independently and uniformly, with
#' replacement and selfing allowed. Recombination switches between that parent's
#' two homologues to transmit one gamete. This is sex-independent diploid
#' Wright--Fisher reproduction. Mutation, parent choice, transmission and sampling
#' translate the corresponding neutral simuPOP operators into R.
#' All operators use R's recorded random-number stream.
#'
#' Sequences start with allele zero at every site. Symmetric biallelic mutation
#' is the default. Mutation acts on the parental generation before reproduction.
#' Adjacent-marker recombination fractions are obtained from Morgan distances
#' with Haldane's mapping function. Different chromosomes assort independently.
#' A specified burn-in length does not certify equilibrium. Compare its
#' diversity record with a longer ancestral run under the same parameters.
#'
#' @param sizes Diploid sizes, beginning with the ancestral size at generation
#'   zero, followed by the size of each offspring generation.
#' @param loci Number of explicitly represented sites on each chromosome.
#'   Invariant sites are retained.
#' @param mutationRate Probability of a zero-to-one mutation per site per generation.
#' @param sampleSize Number of diploid individuals sampled without replacement
#'   from the final generation.
#' @param geneticPositions List of increasing positions in Morgans, one vector
#'   per chromosome. The default has zero recombination within chromosomes.
#' @param positions List of increasing physical positions, one vector per
#'   chromosome. The default is seq_len(loci) for each chromosome.
#' @param burnin Number of ancestral generations at sizes[1].
#' @param founderSize Optional number of diploid individuals sampled without
#'   replacement after burn-in. Their descendants follow sizes[-1].
#' @param reverseMutationRate Probability of a one-to-zero mutation. Defaults
#'   to mutationRate, giving a finite biallelic mutation--drift equilibrium.
#' @param recordEvery Record allele frequencies every this many generations.
#'   The initial and final generations are always recorded.
#' @param seed Positive integer R seed. The caller's RNG state is restored.
#' @return A list with observed final haplotypes and genotypes, calibration,
#'   generation records and demographic truth.
#'   Haplotype rows are paired homologues of sampled individuals. Marker rows
#'   include all represented sites. Generation records contain complete
#'   population frequencies and must not be passed to final-sample estimators.
#'   Recorded heterozygosity is the site mean of 2*p*(1-p). Multiplying by
#'   2*N/(2*N-1) gives mean pairwise difference between distinct population
#'   copies. The finite biallelic model allows recurrent and reverse mutation;
#'   it is not an infinite-sites simulator.
#'   Public sample identifiers are arbitrary labels. Final population indices
#'   are stored only in truth. truth$ancestralSize is the size before founder
#'   sampling; truth$reproductiveSizes gives the actual generation-zero size
#'   after sampling followed by every offspring size.
#' @export
simulatePopulation <- function(sizes, loci, mutationRate, sampleSize,
                               geneticPositions = NULL, positions = NULL,
                               burnin = 0L, founderSize = NULL,
                               reverseMutationRate = mutationRate,
                               recordEvery = 1L, seed = 8149L) {
    if (!is.numeric(sizes) || !length(sizes) || any(!is.finite(sizes)) ||
        any(sizes < 1 | sizes != floor(sizes) | sizes > .Machine$integer.max / 2)) {
        stop("sizes must give positive diploid population sizes", call. = FALSE)
    }
    if (!is.numeric(loci) || !length(loci) || any(!is.finite(loci)) ||
        any(loci < 1 | loci != floor(loci)) || sum(loci) > .Machine$integer.max) {
        stop("loci must give the number of sites on each chromosome", call. = FALSE)
    }
    sizes <- as.integer(sizes)
    loci <- as.integer(loci)
    checkInputs(mutationRate = list(positiveScalar, zero = TRUE),
                reverseMutationRate = list(positiveScalar, zero = TRUE))
    if (mutationRate > 1 || reverseMutationRate > 1) {
        stop("mutation probabilities must be between zero and one", call. = FALSE)
    }
    checkInputs(sampleSize = scalarInteger, burnin = list(scalarInteger, minimum = 0L),
                recordEvery = scalarInteger, seed = scalarInteger)
    sampleSize <- as.integer(sampleSize)
    burnin <- as.integer(burnin)
    recordEvery <- as.integer(recordEvery)
    seed <- as.integer(seed)
    if (!is.null(founderSize)) {
        founderSize <- scalarInteger(founderSize, "founderSize")
        if (founderSize > sizes[1L]) {
            stop("founderSize cannot exceed the ancestral population", call. = FALSE)
        }
    }
    finalSize <- if (length(sizes) == 1L && !is.null(founderSize)) founderSize else
        tail(sizes, 1L)
    if (sampleSize > finalSize) {
        stop("sampleSize cannot exceed the final diploid population", call. = FALSE)
    }
    checkPositions <- function(x, name, genetic = FALSE) {
        if (is.null(x)) return(if (genetic) lapply(loci, function(n) rep(0, n)) else
            lapply(loci, seq_len))
        if (!is.list(x) || length(x) != length(loci) ||
            any(lengths(x) != loci) || any(vapply(x, function(z)
                !is.numeric(z) || any(!is.finite(z)) || any(z < 0) ||
                any(diff(z) < if (genetic) 0 else .Machine$double.eps), logical(1)))) {
            stop(name, " must contain ordered positions for each chromosome",
                 call. = FALSE)
        }
        lapply(x, as.numeric)
    }
    positions <- checkPositions(positions, "positions")
    geneticPositions <- checkPositions(geneticPositions, "geneticPositions", TRUE)
    fractions <- lapply(geneticPositions, function(z)
        c(-expm1(-2 * diff(z)) / 2, 0))
    markers <- data.frame(chromosome = rep(seq_along(loci), loci),
        position = unlist(positions, use.names = FALSE),
        genetic_position = unlist(geneticPositions, use.names = FALSE))
    markers$site <- seq_len(nrow(markers))
    simupopPopulation(sizes, loci, mutationRate, reverseMutationRate,
        sampleSize, burnin, founderSize, recordEvery, seed, markers, fractions)
}

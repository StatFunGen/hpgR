#' Small genealogy and diversity experiments
#'
#' These recipes call shared ancestry and mutation functions with explicit experimental designs. Simulation truth is stored separately from observations. The terminal-extension model adds a fixed duration to terminal branches and fits its spectrum by a Poisson composite likelihood.
#'
#' @name models_genealogy_controls
#' @param seed Nonnegative integer random seed.
#' @param spectra Region-by-derived-count matrix with n-1 columns of nonnegative spectrum counts.
#' @param mutationRate Positive mutation rate per base per generation.
#' @param sequenceLength Positive analyzed length per independent region in bases.
#' @param N Positive integer diploid population size.
#' @param mu Nonnegative mutation rate per base per generation.
#' @param replicates Positive integer number of independent simulated populations or genealogies.
#' @param base_seed Nonnegative integer starting seed for region or replicate assignments.
#' @param tree Compressed genealogy returned by the tree constructors.
#' @param selected Distinct sample indices retained in the subtree.
#' @param Ne Positive diploid effective size for a constant-size Kingman model.
#' @param n Integer number of sampled haploid copies, at least two.
#' @param extension Nonnegative terminal extension in generations before the sample.
#' @param regions Integer number of independent regions, at least two.
#' @param bootstrap Integer number of region-resampling bootstrap replicates, at least two.
#' @param generations Nonnegative integer number of forward reproductive updates.
#' @param G Binary haplotype matrix: sampled copies in rows, derived variants in columns.
#' @returns matchedDraws returns replicate results and spectra. recordedSubset returns a nested tree and ages. sfs returns counts for derived-copy classes 1 through n-1. founderDiversityExample returns replicate nonIBD and diversityRatio arrays, their summary and truth. terminalSpectrumFit returns estimate, expected, constantNe, constantExpected and observed. terminalExtensionExample returns observed regions, spectra, fit, bootstrap intervals and truth. specifiedLocalTreeExample returns observed, calibration and truth with two specified local genealogies.
#' @examples
#' founderDiversityExample(N=4, generations=2, replicates=2, seed=42)$summary
#' sfs(rbind(c(0,1),c(1,1),c(0,0)))
NULL

# Sources: Pritchard 2.2; canonical Notes 8.10--8.12.

#' @rdname models_genealogy_controls
#' @export
matchedDraws <- function(N, mu, replicates = 600L, base_seed) {
    rows <- list(); spectra <- matrix(0, 2 * replicates, 11)
    for (model in 1:2) for (i in seq_len(replicates)) {
        a <- pgRegionSeed(base_seed + model * 10000L, i)
        m <- pgRegionSeed(base_seed + model * 10000L + 500000L, i, mutation = TRUE)
        genealogy <- if (model == 1L) pgDiscreteGenealogy(rep(N, 61L), ancestry_seed = a)
                     else pgKingmanTree(N, ancestry_seed = a)
        summary <- pgTreeSummaries(genealogy)
        mutations <- pgMutateTree(genealogy, 1000, mu, mutation_seed = m)
        row <- (model - 1L) * replicates + i
        spectra[row, ] <- mutations$sfs
        weights <- 2 * (1:11) * (12 - (1:11)) / (12 * 11)
        internal_ages <- genealogy$nodes$time[is.na(genealogy$nodes$sample_id)]
        first_age <- min(internal_ages)
        first_parents <- genealogy$nodes$id[genealogy$nodes$time == first_age]
        first_loss <- sum(table(genealogy$edges$parent)[as.character(first_parents)] - 1L)
        rows[[row]] <- data.frame(N = N, model = c("Discrete WF", "Kingman")[model],
            replicate = i, pair_time = summary$pair_time, TMRCA = summary$TMRCA,
            branch_length = summary$total_branch_length, pi_region = sum(mutations$sfs * weights),
            S = sum(mutations$sfs),
            pair_differences = sum(mutations$haplotypes[1L, ] != mutations$haplotypes[2L, ]),
            first_loss = first_loss,
            multifurcation = summary$max_children > 2L, ancestry_seed = a, mutation_seed = m)
    }
    list(results = do.call(rbind, rows), sfs = spectra)
}

# Founder-label and allele summaries of the same source-derived copying records.
#' @rdname models_genealogy_controls
#' @export
founderDiversityExample <- function(N = 40L, generations = 60L,
                                    replicates = 200L, seed = 814906L) {
    N <- scalarInteger(N, "N")
    generations <- scalarInteger(generations, "generations", 0L)
    replicates <- scalarInteger(replicates, "replicates", 2L)
    seed <- scalarInteger(seed, "seed")
    copies <- 2L * N
    initial <- cbind(seq_len(copies), rep(0:1, each = N))
    H0 <- 2 * N^2 / (copies * (copies - 1))
    nonIBD <- ratio <- matrix(0, replicates, generations + 1L)
    for (i in seq_len(replicates)) {
        record <- simulateForward(initial, rep(N, generations + 1L), seed + i - 1L)
        for (t in seq_along(record$states)) {
            state <- record$states[[t]]
            counts <- tabulate(state[, 1L], copies)
            nonIBD[i, t] <- 1 - sum(counts * (counts - 1)) / (copies * (copies - 1))
            k <- sum(state[, 2L])
            ratio[i, t] <- 2 * k * (copies - k) / (copies * (copies - 1)) / H0
        }
    }
    list(summary = data.frame(generation = 0:generations,
        nonIBD = colMeans(nonIBD), diversityRatio = colMeans(ratio),
        expectation = (1 - 1 / copies)^(0:generations)),
        nonIBD = nonIBD, diversityRatio = ratio,
        truth = list(N = N, H0 = H0, replicates = replicates,
                     seeds = seed + seq_len(replicates) - 1L))
}

# Complete small example. Existing source-derived Kingman and mutation
# functions supply trees and sequences; pgMakeTree adds the stated tip time.
#' @rdname models_genealogy_controls
#' @export
terminalExtensionExample <- function(Ne = 80, n = 12L, mutationRate = 2e-5,
                                     sequenceLength = 1000, extension = 40,
                                     regions = 48L, bootstrap = 120L, seed = 814903L) {
    n <- scalarInteger(n, "n", 3L)
    regions <- scalarInteger(regions, "regions", 2L)
    bootstrap <- scalarInteger(bootstrap, "bootstrap", 2L)
    seed <- scalarInteger(seed, "seed")
    extension <- positiveScalar(extension, "extension", zero = TRUE)
    trees <- lapply(seq_len(regions), function(j) {
        tree <- pgKingmanTree(Ne, n, pgRegionSeed(seed, j))
        nodes <- tree$nodes[c("id", "time")]
        internal <- !nodes$id %in% tree$sample_nodes
        nodes$time[internal] <- nodes$time[internal] + extension
        pgMakeTree(nodes, tree$edges[c("parent", "child")], tree$sample_nodes,
                     tree$sample_ids, "specified terminal extension")
    })
    observed <- lapply(seq_along(trees), function(j) pgMutateTree(trees[[j]],
        sequenceLength, mutationRate, pgRegionSeed(seed + 100000L, j, TRUE)))
    spectra <- t(vapply(observed, function(x) sfs(x$haplotypes), integer(n - 1L)))
    fit <- terminalSpectrumFit(spectra, mutationRate, sequenceLength)
    draws <- pgWithSeed(seed + 1L, function() t(replicate(bootstrap, {
        rows <- sample.int(regions, regions, replace = TRUE)
        terminalSpectrumFit(spectra[rows, , drop = FALSE], mutationRate, sequenceLength)$estimate
    })))
    interval <- apply(draws, 2L, stats::quantile, c(.025, .975))
    list(observed = observed, spectra = spectra, fit = fit, bootstrap = draws,
        interval = interval, truth = list(Ne = Ne, extension = extension,
            mutationRate = mutationRate, sequenceLength = sequenceLength, n = n,
            regions = regions, seed = seed, trees = trees))
}

# Two local genealogies specified before mutation. Node and edge construction
# is a small classroom recipe; mutation placement uses the shared source kernel.
#' @rdname models_genealogy_controls
#' @export
specifiedLocalTreeExample <- function(seed = 814905L) {
    seed <- scalarInteger(seed, "seed")
    make <- function(pair, time) {
        other <- setdiff(1:3, pair)
        pgMakeTree(data.frame(id = 1:7, time = c(rep(0, 4), time, 400, 800)),
            data.frame(parent = c(5L, 5L, 6L, 7L, 6L, 7L),
                       child = c(pair, other, 4L, 5L, 6L)),
            1:4, paste0("S", 1:4), "specified local genealogy")
    }
    trees <- list(left = make(1:2, 120), right = make(c(1L, 3L), 180))
    data <- lapply(seq_along(trees), function(i)
        pgMutateTree(trees[[i]], 50000, 1e-6, mutation_seed = seed + i - 1L))
    position <- c(data[[1L]]$positions_bp, data[[2L]]$positions_bp + 50000)
    order <- order(position)
    observed <- list(haplotypes = cbind(data[[1L]]$haplotypes,
        data[[2L]]$haplotypes)[, order, drop = FALSE],
        positions = position[order], sequenceLength = 100000)
    list(observed = observed, calibration = list(mutationRate = 1e-6),
        truth = list(breakpoint = 50000, TMRCA = 800, trees = trees,
                     mutationSeeds = seed + 0:1))
}

#' A recorded Kingman genealogy with branch mutations
#'
#' Event records for the genealogy animation, with diploid Ne = 1000, mutation rate 1e-8 per base per generation and length 100000 bases. The recorded exponential mergers and branch mutations can be replayed through the shared ancestry and mutation kernels.
#'
#' @name models_coalescent_example
#' @returns List containing seed, n, Ne, mu, L, nodes, edges, merger events and mutations with carrier sets.
#' @examples
#' example <- coalescentExample()
#' example$events[[1]]$wait
NULL

# Source model: Pritchard 2.2, pp.91--94 and 102--103, Fig.2.35;
# Laruson and Reed (2021), Sec.9.1, pp.178--180.
# Stored realization used by the animation. Source-kernel replay tests verify
# the node/edge records and mutation carrier sets.
#' @rdname models_coalescent_example
#' @export
coalescentExample <- function() {
    readRDS(system.file("examples", "coalescent.rds", package = "hpgR", mustWork = TRUE))
}

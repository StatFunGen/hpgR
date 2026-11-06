#' Discrete-copy branch exposure and frequency spectra
#'
#' Generate complete forward parent maps, trace a fixed final sample and count branch copy-generations by descendant count. Independent conditional Poisson mutation counts give the infinite-sites spectrum. This preserves discrete mergers.
#'
#' @name models_forward_sfs
#' @param sizes Nonempty positive integer diploid sizes in forward generation order, beginning with the ancestral population.
#' @param n Integer number of sampled haploid copies, at least two.
#' @param ancestry_seed Nonnegative integer seed for ancestry draws.
#' @param max_extension Nonnegative integer limit for generations added before the stated history to complete sample ancestry.
#' @param sample_indices Distinct one-based indices of fixed terminal haploid copies.
#' @param length_bp Positive analyzed region length in bases.
#' @param mu Nonnegative mutation rate per base per generation.
#' @param mutation_seed Nonnegative integer seed for conditional mutation draws.
#' @returns pgForwardExposure returns branch_generations, root_age, extension_generations, n and ancestry_seed. pgForwardSfs adds sfs, expected_sfs, exposure, analyzed_length_bp, mu and mutation_seed.
#' @examples
#' x <- pgForwardSfs(rep(4L,3), n=4, length_bp=100, mu=1e-3, ancestry_seed=42)
#' x$sfs
NULL

# The exposure-only path keeps parent maps and descendant counts for repeated
# no-recombination predictions, avoiding full node/edge table construction.
# xi_i | maps ~ Poisson(mu * length_bp * branch_generations[i]).
# Aggregating edge counts preserves the law; pgRecordAncestry also stores
# individual edges and permits recombination.

#' @rdname models_forward_sfs
#' @export
pgForwardExposure <- function(sizes, n = 12L, ancestry_seed = 814901L,
                                max_extension = 20000L,
                                sample_indices = seq_len(n)) {
    checkInputs(sizes = list(numericVector, minimum = 1, integer = TRUE,
                            maximum = .Machine$integer.max / 2),
                n = list(numericVector, length = 1L, minimum = 2, integer = TRUE,
                         maximum = 2 * tail(sizes, 1L)),
                max_extension = list(numericVector, length = 1L, minimum = 0, integer = TRUE),
                sample_indices = list(numericVector, length = n, minimum = 1,
                    maximum = 2 * tail(sizes, 1L), integer = TRUE, distinct = TRUE))
    pgWithSeed(ancestry_seed, function() {
        # Generate even maps older than an early sample MRCA: this preserves the
        # exact forward RNG sequence used by pgRecordAncestry() for this seed.
        parent_maps <- simupopParentMaps(sizes)
        copies <- sample_indices
        branch_generations <- numeric(n - 1L)
        root_age <- 0L
        addSampleBranches <- function() {
            # tabulate(copies) counts sampled descendants of each current copy;
            # the outer tabulate counts branches carrying each descendant count.
            branch_generations <<- branch_generations +
                tabulate(tabulate(copies), nbins = n - 1L)
            root_age <<- root_age + 1L
        }
        for (i in rev(seq_along(parent_maps))) {
            if (length(unique(copies)) == 1L) break
            addSampleBranches()
            copies <- parent_maps[[i]][copies]
        }
        extension <- 0L
        while (length(unique(copies)) > 1L) {
            pgCheck(extension < max_extension,
                     paste("root completion limit reached; retain failed seed", ancestry_seed))
            extension <- extension + 1L
            addSampleBranches()
            older_map <- simupopDrawParents(2 * sizes[1L], 2 * sizes[1L])
            copies <- older_map[copies]
        }
        list(branch_generations = branch_generations, root_age = root_age,
             extension_generations = extension, n = n,
             ancestry_seed = ancestry_seed)
    })
}

#' @rdname models_forward_sfs
#' @export
pgForwardSfs <- function(sizes, length_bp = 1000, n = 12L, mu = 2e-5,
                           ancestry_seed = 814901L, mutation_seed = 914901L,
                           max_extension = 20000L,
                           sample_indices = seq_len(n)) {
    checkInputs(length_bp = positiveScalar, mu = list(positiveScalar, zero = TRUE))
    exposure <- pgForwardExposure(sizes, n, ancestry_seed, max_extension,
                                     sample_indices)
    expected <- mu * length_bp * exposure$branch_generations
    pgCheck(all(is.finite(expected)), "mutation exposure is too large")
    sfs <- pgWithSeed(mutation_seed, function() rpois(n - 1L, expected))
    list(sfs = sfs, expected_sfs = expected, exposure = exposure,
         analyzed_length_bp = length_bp, mu = mu, mutation_seed = mutation_seed)
}

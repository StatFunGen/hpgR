# From msprime's msp_run_coalescent and mutgen_place_mutations,
# supplied in Daniel Nachun's tskitR.

#' Simulate constant-size Hudson ancestry in R
#'
#' Translate the neutral single-population Hudson event kernels from the
#' msprime C sources included in Daniel Nachun's tskitR. R segment lists
#' represent ancestral material and event records. Random draws use R's
#' generator with the supplied seed.
#' @param samples Number of present-day diploid sampled individuals.
#' @param populationSize Constant diploid effective size.
#' @param sequenceLength Positive integer sequence length in bases.
#' @param recombinationRate Uniform rate per base per generation.
#'   Ancestral breakpoints fall between integer positions.
#' @param seed Positive integer R seed. Caller RNG state is restored.
#' @param recordFullArg Retain recombination and noncoalescing ancestry records.
#' @return An hpgR_ancestry R list with node and edge tables, zero-based sample
#'   IDs, sequence length, event records and source metadata. This restricted
#'   interface does not implement structured/changing populations, migration,
#'   ancient samples, selection, gene conversion, DTWF or the full msprime API.
#' @export
simulateAncestry <- function(samples, populationSize, sequenceLength,
                             recombinationRate = 0, seed = 8149L,
                             recordFullArg = FALSE) {
    checkInputs(samples = scalarInteger, populationSize = positiveScalar,
                sequenceLength = positiveScalar)
    samples <- as.integer(samples)
    if (sequenceLength != floor(sequenceLength))
        stop("sequenceLength must be an integer for discrete ancestral breakpoints",
             call. = FALSE)
    checkInputs(recombinationRate = list(positiveScalar, zero = TRUE), seed = scalarInteger,
                recordFullArg = logicalScalar)
    seed <- as.integer(seed)
    msprimeAncestry(samples, populationSize, sequenceLength,
        recombinationRate, seed, recordFullArg)
}

#' Add continuous infinite-sites binary mutations in R
#'
#' The edgewise Poisson count and uniform position/time operations translate
#' mutgen_place_mutations in the vendored msprime source. Distinct positions
#' and binary root state zero give one 0-to-1 mutation at each new site.
#' @param ancestry An hpgR_ancestry R object.
#' @param mutationRate Rate per base per generation.
#' @param seed Positive integer R seed. Caller RNG state is restored.
#' @return A new hpgR_ancestry object with R site and mutation tables. Existing
#'   mutations are replaced. The input ancestry is unchanged.
#' @export
mutateAncestry <- function(ancestry, mutationRate, seed = 8150L) {
    checkInputs(mutationRate = list(positiveScalar, zero = TRUE), seed = scalarInteger)
    seed <- as.integer(seed)
    msprimeMutations(ancestry, mutationRate, seed)
}

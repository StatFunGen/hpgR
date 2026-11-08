# From msprime's msp_run_coalescent and mutgen_place_mutations,
# supplied in Daniel Nachun's tskitR.

#' Simulate Hudson ancestry in R
#'
#' Translate the neutral single-population Hudson event kernels from the
#' msprime C sources included in Daniel Nachun's tskitR. R segment lists
#' represent ancestral material and event records. Random draws use R's
#' generator with the supplied seed.
#' @param samples Number of present-day diploid sampled individuals.
#' @param populationSize Diploid effective size before the recorded recent history.
#' @param populationSizes Optional diploid sizes in successive one-generation
#'   epochs backward from sampling. Older epochs use populationSize.
#' @param sequenceLength Positive integer sequence length in bases.
#' @param recombinationRate Uniform rate per base per generation.
#'   Ancestral breakpoints fall between integer positions.
#' @param seed Positive integer R seed. Caller RNG state is restored.
#' @param recordFullArg Retain recombination and noncoalescing ancestry records.
#' @return An hpgR_ancestry R list with node and edge tables, zero-based sample
#'   IDs, sequence length, event records and source metadata. This restricted
#'   interface follows one population with piecewise-constant sizes.
#' @export
simulateAncestry <- function(samples, populationSize, sequenceLength,
                             recombinationRate = 0, seed = 8149L,
                             recordFullArg = FALSE, populationSizes = NULL) {
    checkInputs(samples = scalarInteger, populationSize = positiveScalar,
                sequenceLength = positiveScalar)
    samples <- as.integer(samples)
    if (sequenceLength != floor(sequenceLength))
        stop("sequenceLength must be an integer for discrete ancestral breakpoints",
             call. = FALSE)
    checkInputs(recombinationRate = list(positiveScalar, zero = TRUE), seed = scalarInteger,
                recordFullArg = logicalScalar)
    seed <- as.integer(seed)
    msprimeAncestryCopies(2L * samples, populationSize, sequenceLength,
        recombinationRate, seed, recordFullArg, populationSizes = populationSizes)
}

#' Add mutations to an ancestry record in R
#'
#' Draw Poisson mutation counts on branches, then sample positions and times
#' uniformly along each branch. The ancestral state is 0. Under the binary
#' model, each integer position inside an edge span has the same mutation rate.
#' @param ancestry An hpgR_ancestry R object.
#' @param model Mutation model: infiniteSites uses unique continuous positions;
#'   binary uses integer sites and symmetric recurrent 0/1 mutations.
#' @param mutationRate Rate per base per generation.
#' @param seed Positive integer R seed. Caller RNG state is restored.
#' @return A new hpgR_ancestry object with R site and mutation tables. Existing
#'   mutations are replaced. The input ancestry is unchanged.
#' @export
mutateAncestry <- function(ancestry, mutationRate, seed = 8150L,
                             model = c("infiniteSites", "binary")) {
    checkInputs(mutationRate = list(positiveScalar, zero = TRUE), seed = scalarInteger)
    seed <- as.integer(seed)
    msprimeMutations(ancestry, mutationRate, seed, model = match.arg(model))
}

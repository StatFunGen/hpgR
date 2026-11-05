#' Allele-frequency change under genotype-dependent fitness
#'
#' Random mating supplies genotype probabilities p squared, 2p(1-p), and (1-p) squared. Reproduction weights these by fitnesses w, ordered as two, one and zero copies of the tracked allele.
#'
#' @name models_selection
#' @param p0 Scalar starting allele frequency in [0,1].
#' @param w Three nonnegative finite fitnesses, ordered as two, one and zero tracked copies, with positive mean fitness. balancedFrequency requires heterozygote advantage.
#' @param generations Nonnegative integer number of reproductive updates.
#' @param p Scalar current allele frequency in [0,1].
#' @param Ne Positive integer diploid effective size.
#' @param initial Scalar frequency before reproduction.
#' @param selected Scalar expected frequency after selection.
#' @returns selectedFrequency returns the frequency after fitness weighting. selectionFrequencies returns a vector including generation zero. balancedFrequency returns the internal equilibrium under heterozygote advantage. selectionDeclineProbability returns the probability that a Binomial(2Ne,selected) draw falls below the initial frequency.
#' @examples
#' selectedFrequency(.3, c(1.1,1.05,1))
#' balancedFrequency(c(.25,1,.92))
NULL

# Fitness is ordered as two, one and zero copies of the tracked allele.
#' @rdname models_selection
#' @export
selectionFrequencies <- function(p0, w, generations) {
    .probability(p0, "p0")
    .fitnessWeights(w)
    .scalarInteger(generations, "generations", 0L)
    p <- numeric(generations + 1)
    p[1] <- p0
    for (t in seq_len(generations)) p[t + 1] <- selectedFrequency(p[t], w)
    p
}

#' @rdname models_selection
#' @export
selectionDeclineProbability <- function(Ne, initial, selected) {
    .scalarInteger(Ne, "Ne")
    .probability(initial, "initial")
    .probability(selected, "selected")
    pbinom(ceiling(2 * Ne * initial) - 1, size = 2 * Ne, prob = selected)
}

#' @rdname models_selection
#' @export
balancedFrequency <- function(w) {
    .fitnessWeights(w)
    if (w[2] <= max(w[c(1, 3)]))
        stop("w must give the heterozygote the greatest fitness", call. = FALSE)
    # Equilibrium of the allele in the first homozygote.
    (w[2] - w[3]) / ((w[2] - w[1]) + (w[2] - w[3]))
}

#' @rdname models_selection
#' @export
selectedFrequency <- function(p, w) {
  .probability(p, "p")
  .fitnessWeights(w)
  q <- 1 - p
  mean_fitness <- p^2 * w[1] + 2 * p * q * w[2] + q^2 * w[3]
  if (mean_fitness == 0) stop("the current population must have positive mean fitness", call. = FALSE)
  (p^2 * w[1] + p * q * w[2]) / mean_fitness
}

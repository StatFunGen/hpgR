#' Allele-frequency pooling, migration and drift
#'
#' Pooling combines source Hardy--Weinberg genotype frequencies without remating. Migration exchanges a fraction m symmetrically between two populations; the drift simulator then draws 2Ne allele copies independently within each population.
#'
#' @name models_demographic_processes
#' @param p1 Frequency or frequency vector for population one.
#' @param p2 Frequency or equal-length frequency vector for population two.
#' @param m Scalar symmetric migration fraction in [0,1].
#' @param generations Nonnegative integer number of reproductive updates.
#' @param n Vector of positive integer trial counts.
#' @param confidence Scalar confidence level strictly between zero and one.
#' @param p Source frequency vector in [0,1]; exactly two values for frequencyDifferentiation.
#' @param weights Nonnegative source weights summing to one, one per frequency.
#' @param start Numeric vector of two starting allele frequencies.
#' @param Ne Positive integer diploid population size per population.
#' @returns poolAlleleFrequencies returns source and pooled genotype probabilities, pooled_p, H_S, H_T and F_ST. Migration functions return generation-by-two-population matrices including generation zero. pairwiseFst returns an elementwise vector; frequencyDifferentiation returns p_bar, H_T, H_S and F_ST. F_ST is NA when pooled heterozygosity is zero. zeroCountUpperBound returns one binomial probability bound per n.
#' @examples
#' poolAlleleFrequencies(c(.2,.8))
#' migrationFrequencies(.2,.8,.05,3)
#' zeroCountUpperBound(c(10,100))
NULL

#' @rdname models_demographic_processes
#' @export
poolAlleleFrequencies <- function(p, weights = rep(1 / length(p), length(p))) {
    .probability(p, "p", length = NULL)
    .probability(weights, "weights", length(p))
    if (abs(sum(weights) - 1) >= 1e-10) stop("weights must sum to 1", call. = FALSE)
    genotype <- rbind(AA = p^2, Aa = 2 * p * (1 - p), aa = (1 - p)^2)
    pooled <- as.vector(genotype %*% weights)
    names(pooled) <- rownames(genotype)
    pooled_p <- sum(weights * p)
    H_S <- sum(weights * genotype['Aa', ])
    H_T <- 2 * pooled_p * (1 - pooled_p)
    list(genotype = genotype, pooled = pooled, pooled_p = pooled_p,
        H_S = H_S, H_T = H_T, F_ST = if (H_T > 0) (H_T - H_S) / H_T else NA_real_)
}

#' @rdname models_demographic_processes
#' @export
simulateMigrationDrift <- function(start, Ne, m, generations) {
  .probability(start, "start", 2L)
  .scalarInteger(Ne, "Ne")
  .probability(m, "m")
  .scalarInteger(generations, "generations", 0L)
  frequency <- matrix(NA_real_, nrow = generations + 1, ncol = 2)
  frequency[1, ] <- start
  for (t in seq_len(generations)) {
    mixed <- c((1 - m) * frequency[t, 1] + m * frequency[t, 2],
               (1 - m) * frequency[t, 2] + m * frequency[t, 1])
    frequency[t + 1, ] <- rbinom(2, size = 2 * Ne, prob = mixed) / (2 * Ne)
  }
  frequency
}

#' @rdname models_demographic_processes
#' @export
migrationFrequencies <- function(p1, p2, m, generations) {
  .probability(p1, "p1")
  .probability(p2, "p2")
  .probability(m, "m")
  .scalarInteger(generations, "generations", 0L)
  frequency <- matrix(NA_real_, nrow = generations + 1, ncol = 2)
  frequency[1, ] <- c(p1, p2)
  for (t in seq_len(generations)) {
    frequency[t + 1, ] <- c(
      (1 - m) * frequency[t, 1] + m * frequency[t, 2],
      (1 - m) * frequency[t, 2] + m * frequency[t, 1]
    )
  }
  colnames(frequency) <- c("population 1", "population 2")
  frequency
}

#' @rdname models_demographic_processes
#' @export
pairwiseFst <- function(p1, p2) {
  .probability(p1, "p1", length = NULL)
  .probability(p2, "p2", length(p1))
  pooled <- (p1 + p2) / 2
  h_total <- 2 * pooled * (1 - pooled)
  h_within <- p1 * (1 - p1) + p2 * (1 - p2)
  ifelse(h_total > 0, (h_total - h_within) / h_total, NA_real_)
}

#' @rdname models_demographic_processes
#' @export
frequencyDifferentiation <- function(p) {
  .probability(p, "p", 2L)
  p_bar <- mean(p)
  H_T <- 2 * p_bar * (1 - p_bar)
  H_S <- mean(2 * p * (1 - p))
  c(p_bar = p_bar, H_T = H_T, H_S = H_S,
    F_ST = if (H_T > 0) (H_T - H_S) / H_T else NA_real_)
}

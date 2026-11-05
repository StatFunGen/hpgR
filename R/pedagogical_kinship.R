#' Estimate pairwise kinship from genotype sharing
#'
#' These complete-genotype calculations compare identity by state with reference-frequency expectations or pairwise heterozygosity. They are small moment calculations for the stated model; the PLINK-labelled calculation solves the three IBS equations without constraining the resulting IBD weights.
#'
#' @name models_kinship
#' @param pair Complete numeric matrix with markers in rows and exactly two people in columns.
#' @param frequency One reference allele frequency per marker; frequency-based estimation requires polymorphic markers.
#' @param genotype_i Complete numeric genotype-count vector (0,1,2) for the first person.
#' @param genotype_j Complete numeric genotype-count vector of equal length for the second person.
#' @param p Nonempty vector of marker allele frequencies in [0,1].
#' @returns ibsExpectation returns an IBS-category by IBD-category matrix. plinkKinship returns Z0, Z1, Z2 and kinship. frequencyKinship and kingKinship return scalar estimates; KING returns NA when neither person has a heterozygous marker. compareKinship returns all estimates and their denominator ratio. Pair simulators return a marker-by-two-people genotype matrix.
#' @examples
#' set.seed(42)
#' f <- rep(.3, 100)
#' pair <- simulateSiblingPair(f)
#' compareKinship(pair, f)
NULL

#' @rdname models_kinship
#' @export
ibsExpectation <- function(p) {
  .probability(p, "p", length = NULL)
  r <- 1 - p
  rbind(
    IBS0 = c(mean(2 * p^2 * r^2), 0, 0),
    IBS1 = c(mean(4 * p * r * (p^2 + r^2)), mean(2 * p * r), 0),
    IBS2 = c(mean(p^4 + r^4 + 4 * p^2 * r^2), mean(p^2 + r^2), 1)
  )
}

#' @rdname models_kinship
#' @export
plinkKinship <- function(genotype_i, genotype_j, frequency) {
  .dosageVector(genotype_i, "genotype_i")
  .dosageVector(genotype_j, "genotype_j")
  .sameLength(genotype_i, genotype_j, "genotype_i and genotype_j")
  .probability(frequency, "frequency", length(genotype_i))
  if (!any(frequency > 0 & frequency < 1))
    stop("frequency must include polymorphic markers", call. = FALSE)
  ibs <- ifelse(genotype_i == genotype_j, 2L,
    ifelse(abs(genotype_i - genotype_j) == 2L, 0L, 1L))
  observed <- tabulate(ibs + 1L, nbins = 3L) / length(ibs)
  z <- drop(solve(ibsExpectation(frequency), observed))
  names(z) <- c("Z0", "Z1", "Z2")
  c(z, kinship = unname(z["Z1"] / 4 + z["Z2"] / 2))
}

#' @rdname models_kinship
#' @export
frequencyKinship <- function(genotype_i, genotype_j, frequency) {
  .dosageVector(genotype_i, "genotype_i")
  .dosageVector(genotype_j, "genotype_j")
  .sameLength(genotype_i, genotype_j, "genotype_i and genotype_j")
  .probability(frequency, "frequency", length(genotype_i))
  if (sum(frequency * (1 - frequency)) == 0)
    stop("frequency must include polymorphic markers", call. = FALSE)
  0.5 - sum((genotype_i - genotype_j)^2) /
    (8 * sum(frequency * (1 - frequency)))
}

#' @rdname models_kinship
#' @export
kingKinship <- function(genotype_i, genotype_j) {
  .dosageVector(genotype_i, "genotype_i")
  .dosageVector(genotype_j, "genotype_j")
  .sameLength(genotype_i, genotype_j, "genotype_i and genotype_j")
  heterozygotes <- sum(genotype_i == 1L) + sum(genotype_j == 1L)
  if (heterozygotes == 0) return(NA_real_)
  0.5 - sum((genotype_i - genotype_j)^2) / (2 * heterozygotes)
}

#' @rdname models_kinship
#' @export
compareKinship <- function(pair, frequency) {
  .dosageMatrix(pair, "pair")
  if (ncol(pair) != 2L) stop("pair must have two columns", call. = FALSE)
  genotype_i <- pair[, 1]
  genotype_j <- pair[, 2]
  plink <- plinkKinship(genotype_i, genotype_j, frequency)
  pair_scale <- 2 * (sum(genotype_i == 1L) + sum(genotype_j == 1L))
  frequency_scale <- 8 * sum(frequency * (1 - frequency))
  c(
    PLINK = unname(plink["kinship"]),
    frequency_moment = frequencyKinship(genotype_i, genotype_j, frequency),
    KING = kingKinship(genotype_i, genotype_j),
    KING_to_frequency_scale = pair_scale / frequency_scale
  )
}

#' @rdname models_kinship
#' @export
simulateUnrelatedPair <- function(frequency) {
  .probability(frequency, "frequency", length = NULL)
  cbind(
    rbinom(length(frequency), 2, frequency),
    rbinom(length(frequency), 2, frequency)
  )
}

#' @rdname models_kinship
#' @export
simulateSiblingPair <- function(frequency) {
  .probability(frequency, "frequency", length = NULL)
  n <- length(frequency)
  father <- cbind(rbinom(n, 1, frequency), rbinom(n, 1, frequency))
  mother <- cbind(rbinom(n, 1, frequency), rbinom(n, 1, frequency))

  child_1 <- father[cbind(seq_len(n), sample.int(2, n, replace = TRUE))] +
    mother[cbind(seq_len(n), sample.int(2, n, replace = TRUE))]
  child_2 <- father[cbind(seq_len(n), sample.int(2, n, replace = TRUE))] +
    mother[cbind(seq_len(n), sample.int(2, n, replace = TRUE))]

  cbind(child_1, child_2)
}

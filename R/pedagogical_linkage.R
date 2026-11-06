#' Linkage evidence from recombination and sibling sharing
#'
#' Linkage LOD compares the binomial recombination likelihood with theta = 0.5. Log probabilities keep this comparison stable for large samples. The sibling-sharing calculation assumes independent affected pairs under a Binomial(2, 0.5) sharing null.
#'
#' @name models_linkage
#' @param nonrecombinant Nonnegative integer number of nonrecombinant transmissions.
#' @param recombinant Nonnegative integer number of recombinant transmissions.
#' @param counts Three nonnegative integer affected-pair counts, ordered by 0, 1 and 2 shared copies; their sum must be positive.
#' @param theta Scalar recombination fraction between zero and one half.
#' @returns linkageLOD returns a scalar base-10 log likelihood ratio. fitLinkage returns a one-row data frame with counts, theta_hat constrained to [0,0.5], and max_LOD. siblingSharingTest returns total shared copies, mean sharing and an upper-tail null probability.
#' @examples
#' fitLinkage(9, 1)
#' linkageLOD(.1, 9000, 1000)
#' siblingSharingTest(c(2, 8, 10))
NULL

#' @rdname models_linkage
#' @export
siblingSharingTest <- function(counts) {
    stopifnot(length(counts) == 3L, all(is.finite(counts)),
        all(counts >= 0), all(counts == floor(counts)), sum(counts) > 0)
    n <- sum(counts)
    total <- sum((0:2) * counts)
    c(total = total, mean = total / n,
        tail = pbinom(total - 1, size = 2 * n, prob = .5, lower.tail = FALSE))
}

#' @rdname models_linkage
#' @export
linkageLOD <- function(theta, nonrecombinant, recombinant) {
  checkInputs(theta = probability)
  if (theta > 0.5) stop("theta must be between 0 and 0.5", call. = FALSE)
  checkInputs(nonrecombinant = list(scalarInteger, minimum = 0L),
              recombinant = list(scalarInteger, minimum = 0L))
  nonrecombinant <- as.integer(nonrecombinant)
  recombinant <- as.integer(recombinant)
  n <- nonrecombinant + recombinant
  (dbinom(recombinant, n, theta, log = TRUE) -
     dbinom(recombinant, n, 0.5, log = TRUE)) / log(10)
}

#' @rdname models_linkage
#' @export
fitLinkage <- function(nonrecombinant, recombinant) {
  checkInputs(nonrecombinant = list(scalarInteger, minimum = 0L),
              recombinant = list(scalarInteger, minimum = 0L))
  nonrecombinant <- as.integer(nonrecombinant)
  recombinant <- as.integer(recombinant)
  n <- nonrecombinant + recombinant
  if (n == 0L) stop("at least one informative meiosis is required", call. = FALSE)
  theta_hat <- min(recombinant / n, 0.5)
  data.frame(
    nonrecombinant = nonrecombinant,
    recombinant = recombinant,
    theta_hat = theta_hat,
    max_LOD = linkageLOD(theta_hat, nonrecombinant, recombinant)
  )
}

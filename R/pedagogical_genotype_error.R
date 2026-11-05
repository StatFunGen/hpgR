#' Estimate allele frequency with genotype-call error
#'
#' True genotypes follow Hardy--Weinberg probabilities for counts 0, 1 and 2. An error event redraws the reported count uniformly from these three values, including the true value. Thus an error-event rate epsilon gives an incorrect-call probability 2 epsilon / 3.
#'
#' @name models_genotype_error
#' @param d Nonempty numeric vector of observed genotype counts 0, 1 or 2.
#' @param f Scalar allele frequency; the EM default is the observed mean dosage divided by two.
#' @param epsilon Probability of a uniform genotype redraw, between zero and one.
#' @param n Positive integer number of people.
#' @param frequency Scalar allele frequency between zero and one.
#' @param error Probability of a uniform genotype redraw.
#' @param tolerance Positive maximum absolute parameter change for convergence.
#' @param max_iterations Positive integer iteration limit.
#' @param g Nonempty numeric vector of genotype counts 0, 1 or 2.
#' @returns hweProbability returns three named genotype probabilities. genotypeLogLik and genotypeErrorLogLik return scalar log likelihoods. genotypeErrorMatrix is observed-call by true-genotype; Joint and Posterior are person by true-genotype matrices. genotypeErrorUpdate returns updated f and epsilon. fitGenotypeError returns f, epsilon, posterior, history, converged and iterations, warning at the iteration limit. simulateGenotypeCalls returns truth, observed and error_event vectors.
#' @examples
#' hweProbability(.3)
#' set.seed(42)
#' calls <- simulateGenotypeCalls(100, .3, .05)
#' genotypeErrorPosterior(calls$observed[1:4], .3, .05)
NULL

#' @rdname models_genotype_error
#' @export
simulateGenotypeCalls <- function(n, frequency, error) {
    n <- scalarInteger(n, "n")
    probability(frequency, "frequency")
    probability(error, "error")
    truth <- sample(0:2, n, replace = TRUE, prob = hweProbability(frequency))
    observed <- truth
    error_event <- runif(n) < error
    observed[error_event] <- sample(0:2, sum(error_event), replace = TRUE)
    list(truth = truth, observed = observed, error_event = error_event)
}

#' @rdname models_genotype_error
#' @export
hweProbability <- function(f) {
  probability(f, "f")
  c(G0 = (1 - f)^2, G1 = 2 * f * (1 - f), G2 = f^2)
}

#' @rdname models_genotype_error
#' @export
genotypeLogLik <- function(f, g) {
  dosageVector(g, "g")
  sum(log(hweProbability(f)[g + 1]))
}

#' @rdname models_genotype_error
#' @export
genotypeErrorMatrix <- function(epsilon) {
  probability(epsilon, "epsilon")
  out <- matrix(epsilon / 3, 3, 3,
                dimnames = list(observed = 0:2, true = 0:2))
  diag(out) <- 1 - 2 * epsilon / 3
  out
}

#' @rdname models_genotype_error
#' @export
genotypeErrorJoint <- function(d, f, epsilon) {
  dosageVector(d, "d")
  genotypeErrorMatrix(epsilon)[d + 1, , drop = FALSE] *
    rep(hweProbability(f), each = length(d))
}

#' @rdname models_genotype_error
#' @export
genotypeErrorPosterior <- function(d, f, epsilon) {
  joint <- genotypeErrorJoint(d, f, epsilon)
  if (any(rowSums(joint) == 0))
    stop("observed calls have zero probability under this model", call. = FALSE)
  joint / rowSums(joint)
}

#' @rdname models_genotype_error
#' @export
genotypeErrorLogLik <- function(d, f, epsilon) {
  sum(log(rowSums(genotypeErrorJoint(d, f, epsilon))))
}

#' @rdname models_genotype_error
#' @export
genotypeErrorUpdate <- function(d, f, epsilon) {
  weight <- genotypeErrorPosterior(d, f, epsilon)
  expected_dosage <- as.numeric(weight %*% (0:2))
  posterior_match <- weight[cbind(seq_along(d), d + 1)]
  c(
    f = mean(expected_dosage) / 2,
    epsilon = min(1 - 1e-10, 1.5 * mean(1 - posterior_match))
  )
}

#' @rdname models_genotype_error
#' @export
fitGenotypeError <- function(d, f = mean(d) / 2, epsilon = 0.05,
                   tolerance = 1e-10, max_iterations = 10000) {
  dosageVector(d, "d")
  probability(f, "f")
  probability(epsilon, "epsilon")
  positiveScalar(tolerance, "tolerance")
  max_iterations <- scalarInteger(max_iterations, "max_iterations")
  converged <- FALSE
  history <- data.frame(
    iteration = 0, f = f, epsilon = epsilon,
    loglik = genotypeErrorLogLik(d, f, epsilon)
  )

  for (iteration in seq_len(max_iterations)) {
    updated <- genotypeErrorUpdate(d, f, epsilon)
    f_new <- updated[["f"]]
    epsilon_new <- updated[["epsilon"]]
    history <- rbind(history, data.frame(
      iteration = iteration, f = f_new, epsilon = epsilon_new,
      loglik = genotypeErrorLogLik(d, f_new, epsilon_new)
    ))

    change <- max(abs(c(f_new - f, epsilon_new - epsilon)))
    f <- f_new
    epsilon <- epsilon_new
    if (change < tolerance) {
      converged <- TRUE
      break
    }
  }

  if (!converged) warning("Genotype EM reached max_iterations before convergence",
                          call. = FALSE)
  list(
    f = f, epsilon = epsilon,
    posterior = genotypeErrorPosterior(d, f, epsilon),
    history = history, converged = converged, iterations = iteration
  )
}

#' Fit a two-group normal mixture
#'
#' The fitted mixture has two means and common variance one. Responsibilities are posterior probabilities of group one. The simulator fixes equal group counts, then randomly orders the observations.
#'
#' @name models_em
#' @param observed Finite numeric observations.
#' @param pi Mixing probability of group one, between zero and one.
#' @param mu Numeric vector of two means, ordered as groups zero and one.
#' @param start List with an interior mixing probability pi and two means mu.
#' @param tolerance Positive maximum absolute parameter change for convergence.
#' @param max_iterations Positive integer iteration limit.
#' @param n Even integer number of observations, at least two.
#' @param means Numeric vector of two generating means.
#' @param sd Positive common generating standard deviation; fitting assumes sd = 1.
#' @returns simulateGaussianMixture returns observed values and binary group labels. gaussianResponsibilities returns a probability vector; gaussianMixtureLogLik returns a scalar log likelihood. gaussianMixtureUpdate returns pi and mu. fitGaussianMixture also returns history, converged and iterations; history includes iteration zero. An exhausted iteration limit raises a warning.
#' @examples
#' set.seed(42)
#' x <- simulateGaussianMixture(40, c(-2, 2))
#' fit <- fitGaussianMixture(x$observed, list(pi=.5, mu=c(-1, 1)))
#' fit$mu
NULL

#' @rdname models_em
#' @export
simulateGaussianMixture <- function(n, means, sd = 1) {
    checkInputs(n = list(scalarInteger, minimum = 2L))
    n <- as.integer(n)
    if (n %% 2L != 0L) stop("n must be even for equal group sizes", call. = FALSE)
    checkInputs(means = list(numericVector, length = 2L), sd = positiveScalar)
    group <- sample(rep(c(0, 1), each = n / 2))
    observed <- rnorm(n, mean = means[group + 1], sd = sd)
    list(observed = observed, group = group)
}

#' @rdname models_em
#' @export
gaussianResponsibilities <- function(observed, pi, mu) {
  checkInputs(observed = numericVector, pi = probability, mu = list(numericVector, length = 2L))
  log_odds <- log(pi) - log1p(-pi) + dnorm(observed, mu[2], 1, log = TRUE) -
    dnorm(observed, mu[1], 1, log = TRUE)
  plogis(log_odds)
}

#' @rdname models_em
#' @export
gaussianMixtureLogLik <- function(observed, pi, mu) {
  checkInputs(observed = numericVector, pi = probability, mu = list(numericVector, length = 2L))
  log_0 <- log1p(-pi) + dnorm(observed, mu[1], 1, log = TRUE)
  log_1 <- log(pi) + dnorm(observed, mu[2], 1, log = TRUE)
  largest <- pmax(log_0, log_1)
  sum(largest + log(exp(log_0 - largest) + exp(log_1 - largest)))
}

#' @rdname models_em
#' @export
gaussianMixtureUpdate <- function(observed, pi, mu) {
  r <- gaussianResponsibilities(observed, pi, mu)
  weight_1 <- sum(r)
  weight_0 <- sum(1 - r)
  if (weight_1 == 0 || weight_0 == 0)
    stop("both mixture groups need positive posterior weight", call. = FALSE)
  list(
    pi = mean(r),
    mu = c(
      sum((1 - r) * observed) / weight_0,
      sum(r * observed) / weight_1
    )
  )
}

#' @rdname models_em
#' @export
fitGaussianMixture <- function(observed, start, tolerance = 1e-8,
                               max_iterations = 500L) {
  checkInputs(observed = numericVector, tolerance = positiveScalar,
              max_iterations = scalarInteger)
  max_iterations <- as.integer(max_iterations)
  if (!is.list(start) || !all(c("pi", "mu") %in% names(start)))
    stop("start must contain pi and two group means in mu", call. = FALSE)
  probability(start$pi, "start$pi")
  numericVector(start$mu, "start$mu", length = 2L)
  if (start$pi == 0 || start$pi == 1)
    stop("start$pi must be strictly between zero and one", call. = FALSE)
  converged <- FALSE
  parameter <- start
  history <- data.frame(
    iteration = 0L, pi = parameter$pi,
    mu_0 = parameter$mu[1], mu_1 = parameter$mu[2],
    loglik = gaussianMixtureLogLik(observed, parameter$pi, parameter$mu)
  )

  for (iteration in seq_len(max_iterations)) {
    updated <- gaussianMixtureUpdate(observed, parameter$pi, parameter$mu)
    change <- max(abs(
      c(updated$pi, updated$mu) - c(parameter$pi, parameter$mu)
    ))
    parameter <- updated
    history <- rbind(
      history,
      data.frame(
        iteration = iteration, pi = parameter$pi,
        mu_0 = parameter$mu[1], mu_1 = parameter$mu[2],
        loglik = gaussianMixtureLogLik(
          observed, parameter$pi, parameter$mu
        )
      )
    )
    if (change < tolerance) {
      converged <- TRUE
      break
    }
  }

  parameter$history <- history
  parameter$converged <- converged
  parameter$iterations <- iteration
  if (!converged) warning("Gaussian EM reached max_iterations before convergence",
                          call. = FALSE)
  parameter
}

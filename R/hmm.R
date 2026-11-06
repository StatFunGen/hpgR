# From Daniel Nachun's haplotypeTools: HMMProbSpace::forward,
# HMMProbSpace::backward and HMMProbSpace::posterior (GPL >= 3).

#' State probabilities in a finite hidden Markov model
#'
#' @param emission Nonnegative observation-by-state likelihood matrix.
#' @param transition Row-stochastic state transition matrix, or an array with
#'   dimensions states, states, observations minus one.
#' @param initial Initial state probabilities. Defaults to equal weights.
#' @return A list with alpha, beta, posterior and log_likelihood.
#' @export
hmmPosterior <- function(emission, transition, initial = NULL) {
    checkInputs(emission = list(numericMatrix, minimum = 0),
                transition = list(transitionProbabilities, states = ncol(emission),
                                  observations = nrow(emission)))
    n <- nrow(emission)
    k <- ncol(emission)
    varying <- length(dim(transition)) == 3L
    if (is.null(initial)) initial <- rep(1 / k, k)
    checkInputs(initial = list(probability, length = k, sumToOne = TRUE))
    alpha <- beta <- matrix(0, n, k, dimnames = dimnames(emission))
    scale <- numeric(n)
    alpha[1L, ] <- initial * emission[1L, ]
    scale[1L] <- sum(alpha[1L, ])
    if (scale[1L] <= 0) stop('observations have zero probability')
    alpha[1L, ] <- alpha[1L, ] / scale[1L]
    if (n > 1L) for (j in 2:n) {
        tr <- if (varying) matrix(transition[, , j - 1L], k, k) else transition
        alpha[j, ] <- as.numeric(alpha[j - 1L, ] %*% tr) * emission[j, ]
        scale[j] <- sum(alpha[j, ])
        if (scale[j] <= 0) stop('observations have zero probability')
        alpha[j, ] <- alpha[j, ] / scale[j]
    }
    beta[n, ] <- 1 / k
    if (n > 1L) for (j in rev(seq_len(n - 1L))) {
        tr <- if (varying) matrix(transition[, , j], k, k) else transition
        beta[j, ] <- as.numeric(tr %*% (emission[j + 1L, ] * beta[j + 1L, ]))
        beta[j, ] <- beta[j, ] / sum(beta[j, ])
    }
    posterior <- alpha * beta
    posterior <- posterior / rowSums(posterior)
    list(alpha = alpha, beta = beta, posterior = posterior,
         log_likelihood = sum(log(scale)))
}

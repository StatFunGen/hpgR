# Daniel Nachun, haplotypeTools, GPL >= 3.
# R translation of HMMProbSpace::forward, backward and posterior in
# inst/include/haplotypeTools/hmm_engine_prob.h. See source-map.tsv.
# Dense transition steps use matrix multiplication. Backward messages use
# the source's self-normalized mode. Rows are observations, columns states.

#' State probabilities in a finite hidden Markov model
#'
#' @param emission Nonnegative observation-by-state likelihood matrix.
#' @param transition Row-stochastic state transition matrix, or an array with
#'   dimensions states, states, observations minus one.
#' @param initial Initial state probabilities. Defaults to equal weights.
#' @return A list with alpha, beta, posterior and log_likelihood.
#' @export
hmmPosterior <- function(emission, transition, initial = NULL) {
    if (!is.matrix(emission) || !is.numeric(emission) ||
        !nrow(emission) || !ncol(emission) || any(!is.finite(emission)) ||
        any(emission < 0)) stop('emission must be a nonnegative likelihood matrix')
    n <- nrow(emission)
    k <- ncol(emission)
    varying <- length(dim(transition)) == 3L
    expected <- if (varying) c(k, k, n - 1L) else c(k, k)
    if (!identical(as.integer(dim(transition)), as.integer(expected)) ||
        any(!is.finite(transition)) || any(transition < 0))
        stop('transition dimensions or probabilities are invalid')
    if (varying && n > 1L) {
        for (j in seq_len(n - 1L))
            if (any(abs(rowSums(matrix(transition[, , j], k, k)) - 1) > 1e-10))
                stop('transition rows must sum to one')
    } else if (!varying && any(abs(rowSums(transition) - 1) > 1e-10)) {
        stop('transition rows must sum to one')
    }
    if (is.null(initial)) initial <- rep(1 / k, k)
    if (length(initial) != k || any(!is.finite(initial)) || any(initial < 0) ||
        abs(sum(initial) - 1) > 1e-10) stop('initial probabilities must sum to one')
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

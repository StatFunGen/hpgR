# From Daniel Nachun's demographyTools: .w_matrix, .expm1d,
# .etjj_constant and .truncated_sfs.

#' Predict the spectrum under constant population size
#'
#' Compute expected segregating-site counts in the neutral, constant-size,
#' infinite-sites coalescent model. Ne counts diploid individuals. A count of
#' i denotes a derived allele carried by i of n sampled haploid sequences.
#' @param n Number of sampled haploid sequences, at least two.
#' @param Ne Diploid effective population size.
#' @param mu Mutation rate per base per generation.
#' @param length_bp Analyzed sequence length, including invariant bases.
#' @return Numeric vector for derived counts 1 through n - 1.
#' @examples
#' expectedSFS(n = 6, Ne = 10000, mu = 1e-8, length_bp = 10000)
#' n <- 12L
#' Ne <- 10000
#' mu <- 1e-8
#' L <- 10000
#' prediction <- expectedSFS(n, Ne, mu, L)
#' stopifnot(isTRUE(all.equal(prediction, 4 * Ne * mu * L / seq_len(n - 1L))))
#' @export
expectedSFS <- function(n, Ne, mu, length_bp) {
    checkInputs(n = list(scalarInteger, minimum = 2L), Ne = positiveScalar,
                mu = list(positiveScalar, zero = TRUE), length_bp = positiveScalar)
    n <- as.integer(n)
    expected <- truncatedSfs(n, 2 * Ne, Inf)
    mu * length_bp * expected[seq.int(2L, n)]
}

wMatrix <- function(n) {
    if (n <= 1L) return(matrix(0, 0L, 0L))
    W <- matrix(0, n - 1L, n - 1L)
    bb <- seq_len(n - 1L)
    W[, 1L] <- 6 / (n + 1)
    if (n > 2L) {
        W[, 2L] <- 30 * (n - 2 * bb) / ((n + 1) * (n + 2))
        if (n > 3L) {
            for (j in 2L:(n - 2L)) {
                jcol <- j + 1L
                W[, jcol] <- W[, jcol - 2L] *
                    (-(1 + j) * (3 + 2 * j) * (n - j) /
                     (j * (2 * j - 1) * (n + j + 1))) +
                    W[, jcol - 1L] *
                    ((3 + 2 * j) * (n - 2 * bb) /
                     (j * (n + j + 1)))
            }
        }
    }
    t(W)
}

expm1d <- function(x) {
    ifelse(abs(x) < 1e-6, {
        c_n <- 1; ret <- 1
        for (n in 2:10) { c_n <- c_n * x / n; ret <- ret + c_n }
        ret
    }, expm1(x) / x)
}

etjjConstant <- function(n, Ne, tau) {
    j <- seq(2L, n)
    coal_rate <- j * (j - 1) / (2 * Ne)
    if (is.infinite(tau)) {
        return(1 / coal_rate)
    }
    scaled_time <- coal_rate * tau
    tau * expm1d(-scaled_time)
}

truncatedSfs <- function(n, Ne, tau) {
    if (n == 0L) return(0)
    etjj <- etjjConstant(n, Ne, tau)
    W <- wMatrix(n)
    sfs_mid <- as.numeric(colSums(etjj * W))
    k <- seq_len(n - 1L)
    before_tmrca <- tau - sum(sfs_mid * k / n)
    if (is.infinite(tau)) before_tmrca <- 0
    c(0, sfs_mid, before_tmrca)
}

# From demographyTools::.sfs_expected_single_pop.
sfsExpectedSingle <- function(epochs, n, mu) {
    if (n < 2L) return(numeric(0))
    pop <- "population"
    n_epochs <- nrow(epochs)

    j <- seq(2L, n)
    coal_rate <- j * (j - 1) / 2
    etjj_total <- numeric(n - 1L)
    no_coal_prob <- rep(1, n - 1L)

    for (k in seq_len(n_epochs)) {
        Ne <- epochs[[pop]][k]
        if (k < n_epochs) {
            tau <- epochs$time_generations[k + 1L] -
                   epochs$time_generations[k]
        } else {
            tau <- Inf
        }
        etjj_epoch <- etjjConstant(n, Ne, tau)
        etjj_total <- etjj_total + no_coal_prob * etjj_epoch
        if (is.finite(tau)) {
            scaled_time <- coal_rate / Ne * tau
            no_coal_prob <- no_coal_prob * exp(-scaled_time)
        }
    }

    W <- wMatrix(n)
    sfs <- as.numeric(colSums(etjj_total * W)) * mu
    sfs
}


#' Predict the SFS along a specified size history
#'
#' From demographyTools::.sfs_expected_single_pop. The neutral infinite-sites
#' model accumulates coalescent waiting-time components over the size epochs.
#' The oldest population is at mutation-drift equilibrium.
#' @inheritParams expectedSFS
#' @param sizes Diploid sizes in forward order. The first is the ancestral
#'   size; each subsequent entry is the parental size for a recent generation.
#' @return Expected segregating-site counts for derived counts 1 to n - 1.
#' @export
expectedHistorySFS <- function(n, sizes, mu, length_bp) {
    checkInputs(n = list(scalarInteger, minimum = 2L),
        sizes = list(numericVector, minimum = .Machine$double.xmin),
        mu = list(positiveScalar, zero = TRUE), length_bp = positiveScalar)
    epochs <- data.frame(time_generations = seq_along(sizes) - 1L,
                         population = 2 * c(rev(sizes[-1L]), sizes[1L]))
    sfsExpectedSingle(epochs, as.integer(n), mu * length_bp)
}

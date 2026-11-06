#' Construct binary and donor-copying HMM factors
#'
#' Donor redraws sample from the supplied weights and can select the current donor again. Target alleles equal the donor allele with probability 1 - error. Each missing target allele has emission likelihood one for every donor.
#'
#' @name models_hmm
#' @param observed Binary target-allele vector; copyingEmission also accepts NA for an untyped marker.
#' @param switch_probability Scalar probability of changing binary state at each interval.
#' @param emission_error Scalar probability of a binary observation error.
#' @param frequency Nonempty vector of marker allele frequencies in [0,1].
#' @param copy_probability Scalar probability of copying the preceding allele rather than drawing from the current marker frequency.
#' @param reference Binary numeric matrix with markers in rows and reference donors in columns.
#' @param error Scalar probability of flipping the copied allele.
#' @param redraw Scalar probability of redrawing a donor at each interval.
#' @param stay Scalar probability of retaining the donor before a possible redraw.
#' @param weights Nonempty donor probability vector summing to one.
#' @returns binaryHMM returns posterior probabilities of state one, in observation order. copyingTransition returns a row-stochastic donor-by-donor matrix. copyingEmission returns a marker-by-donor likelihood matrix. simulateCopyingPath returns one-based donor indices in path and binary copied alleles in haplotype. simulateReferenceHaplotype returns a binary vector.
#' @examples
#' binaryHMM(c(0, 0, 1, 0, 1), .1, .05)
#' tr <- copyingTransition(.9, c(.4, .6))
#' rowSums(tr)
#' copyingEmission(cbind(c(0,1),c(1,1)), c(NA,1), .02)
NULL

# Observation probabilities for binary states and a symmetric calling error.
#' @rdname models_hmm
#' @export
binaryHMM <- function(observed, switch_probability, emission_error) {
    checkInputs(observed = binaryVector, switch_probability = probability,
                emission_error = probability)
    transition <- matrix(c(1 - switch_probability, switch_probability,
        switch_probability, 1 - switch_probability), 2, byrow = TRUE)
    emission <- cbind(ifelse(observed == 0, 1 - emission_error, emission_error),
        ifelse(observed == 1, 1 - emission_error, emission_error))
    hmmPosterior(emission, transition)$posterior[, 2L]
}

#' Infer local ancestry from source allele frequencies
#'
#' A haploid allele has Bernoulli likelihood given its ancestry source.
#' Across genetic distance d, the state persists with probability exp(-redrawRate*d);
#' otherwise it is redrawn from initial. Factors use copyingTransition and
#' hmmPosterior. Missing alleles have likelihood 1.
#' @param observed Numeric binary allele vector, with NA for untyped markers.
#' @param frequency Marker-by-source matrix of allele-1 frequencies.
#' @param geneticPositions Nondecreasing marker positions in Morgans.
#' @param redrawRate Nonnegative ancestry-redraw rate per Morgan.
#' @param initial Source probabilities. Defaults to equal weights.
#' @return The hmmPosterior list of alpha, beta, posterior and log_likelihood;
#'   posterior rows are markers and columns follow the source order in frequency.
#' @examples
#' F <- cbind(A=c(.2,.2), B=c(.8,.8))
#' ancestryHMM(c(1,1), F, c(0,1), -log(.8))$posterior
#' @export
ancestryHMM <- function(observed, frequency, geneticPositions,
                        redrawRate = 8, initial = NULL) {
    checkInputs(frequency = list(numericMatrix, minimum = 0, maximum = 1),
                observed = list(binaryVector, length = nrow(frequency), missing = TRUE, logical = TRUE),
                geneticPositions = list(numericVector, length = nrow(frequency), nondecreasing = TRUE),
                redrawRate = list(positiveScalar, zero = TRUE))
    k <- ncol(frequency)
    if (is.null(initial)) initial <- rep(1 / k, k)
    checkInputs(initial = list(probability, length = k))
    if (abs(sum(initial) - 1) > 1e-10)
        stop("initial must sum to one", call. = FALSE)
    emission <- matrix(1, nrow(frequency), k, dimnames = dimnames(frequency))
    one <- which(!is.na(observed) & observed == 1)
    zero <- which(!is.na(observed) & observed == 0)
    emission[one, ] <- frequency[one, , drop = FALSE]
    emission[zero, ] <- 1 - frequency[zero, , drop = FALSE]
    stay <- exp(-redrawRate * diff(geneticPositions))
    transition <- array(0, c(k, k, length(stay)))
    for (j in seq_along(stay)) transition[, , j] <- copyingTransition(stay[j], initial)
    hmmPosterior(emission, transition, initial)
}

# Redraw a donor from its population weights; a redraw can keep the same donor.
#' @rdname models_hmm
#' @export
copyingTransition <- function(stay, weights) {
    checkInputs(stay = probability, weights = list(probability, length = NULL, sumToOne = TRUE))
    k <- length(weights)
    stay * diag(k) + (1 - stay) * matrix(rep(weights, each = k), nrow = k)
}

#' @rdname models_hmm
#' @export
copyingEmission <- function(reference, observed, error) {
    checkInputs(reference = binaryMatrix,
                observed = list(binaryVector, length = nrow(reference), missing = TRUE),
                error = probability)
    out <- matrix(1, nrow(reference), ncol(reference),
        dimnames = list(marker = rownames(reference), state = colnames(reference)))
    typed <- which(!is.na(observed))
    for (j in typed) out[j, ] <- ifelse(reference[j, ] == observed[j], 1 - error, error)
    out
}

#' @rdname models_hmm
#' @export
simulateCopyingPath <- function(reference, redraw, error) {
    checkInputs(reference = binaryMatrix, redraw = probability, error = probability)
    M <- nrow(reference)
    K <- ncol(reference)
    path <- target <- integer(M)
    current_donor <- sample.int(K, 1)
    for (s in seq_len(M)) {
        if (s > 1 && runif(1) < redraw) current_donor <- sample.int(K, 1)
        path[s] <- current_donor
        target[s] <- reference[s, current_donor]
        if (runif(1) < error) target[s] <- 1 - target[s]
    }
    list(path = path, haplotype = target)
}

#' @rdname models_hmm
#' @export
simulateReferenceHaplotype <- function(frequency, copy_probability = 0.70) {
  checkInputs(frequency = list(probability, length = NULL), copy_probability = probability)
  haplotype <- integer(length(frequency))
  haplotype[1] <- rbinom(1, 1, frequency[1])
  if (length(haplotype) > 1L) for (s in 2:length(haplotype)) {
    haplotype[s] <- if (runif(1) < copy_probability) {
      haplotype[s - 1]
    } else {
      rbinom(1, 1, frequency[s])
    }
  }
  haplotype
}

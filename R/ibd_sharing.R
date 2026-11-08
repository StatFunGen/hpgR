# From demographyTools: compute_prob_tmrca and compute_exp_ibd
# in src/ibd_ne_interface.cpp (Daniel Nachun; GPL >= 3).

#' Expected length in long IBD segments
#'
#' A lineage pair first shares an ancestor in generation g with probability
#' prod(1 - 1/(2 Ne[k]), k < g) / (2 Ne[g]). The expected shared length
#' includes chromosome ends.
#' @param effectiveSizes Diploid effective sizes for generations 1, 2, ... .
#' @param minimumCm Minimum segment length in centimorgans.
#' @param chromosomeLengthsCm Analyzed chromosome lengths in centimorgans.
#' @param haplotypePairs Number of haplotype pairs compared.
#' @return Expected total shared cM by generation. Sum to obtain total sharing
#'   above the threshold. Extend effectiveSizes until the truncated tail is negligible.
#' @export
expectedIbdLength <- function(effectiveSizes, minimumCm, chromosomeLengthsCm,
                               haplotypePairs = 1) {
    checkInputs(effectiveSizes = list(numericVector, minimum = .5),
        minimumCm = list(positiveScalar, zero = TRUE),
        chromosomeLengthsCm = list(numericVector, minimum = 0),
        haplotypePairs = positiveScalar)
    beta <- .5 / effectiveSizes
    probability <- c(1, head(cumprod(1 - beta), -1L)) * beta
    generation <- seq_along(beta)
    cg50 <- minimumCm * generation / 50
    lengths <- chromosomeLengthsCm[chromosomeLengthsCm >= minimumCm]
    haplotypePairs * probability * exp(-cg50) *
        ((cg50 + 1) * sum(lengths) - minimumCm * cg50 * length(lengths))
}

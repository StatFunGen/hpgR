# From tskit's trees.c, supplied in Daniel Nachun's tskitR:
# diversity_summary_func, D_summary_func, r2_summary_func,
# tsk_treeseq_update_site_afs and fold.

tskitDiversity <- function(count, n) {
    # Sum the per-allele contribution over the two alleles.
    2 * count * (n - count) / (n * (n - 1))
}

tskitD <- function(pAB, pA, pB) pAB - pA * pB

tskitR2 <- function(pAB, pA, pB) {
    D <- tskitD(pAB, pA, pB)
    D * D / (pA * pB * (1 - pA) * (1 - pB))
}

tskitSfs <- function(count, n, polarised = TRUE) {
    # Polarised sites contribute the derived count. In the folded spectrum,
    # each allele contributes 1/2 to the same minor-count bin.
    count <- count[count > 0 & count < n]
    if (polarised) tabulate(count, n - 1L) else
        tabulate(pmin(count, n - count), floor(n / 2))
}

#' Four-population allele-frequency statistic
#'
#' Average the biallelic site contributions to f4. The two allele
#' orientations sum to (pA - pB) * (pC - pD) at each site.
#' @param pA,pB,pC,pD Frequencies of the same allele at the same sites in
#'   four disjoint samples, in the order f4(A, B; C, D).
#' @return The mean f4 contribution per supplied site.
#' @source Tskit's f4_summary_func in trees.c, supplied in Daniel Nachun's tskitR.
#' @examples
#' frequencyF4(c(.8, .7), c(.2, .3), c(.9, .6), c(.1, .5))
#' @export
frequencyF4 <- function(pA, pB, pC, pD) {
    checkInputs(pA = list(probability, length = NULL),
        pB = list(probability, length = length(pA)),
        pC = list(probability, length = length(pA)),
        pD = list(probability, length = length(pA)))
    difference <- pC - pD
    mean(pA * (1 - pB) * difference - (1 - pA) * pB * difference)
}

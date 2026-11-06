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

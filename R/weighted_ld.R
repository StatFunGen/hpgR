# From Daniel Nachun's ancestryTools: assign_cm_bins,
# aggregate_weighted_ld_bins and compute_pair_count_histogram.

#' Assign markers to genetic-distance bins
#' @param cM Marker positions in centimorgans.
#' @param cm_min Left edge of the first bin.
#' @param bin_width_cm Bin width in centimorgans.
#' @param n_bins Number of bins.
#' @return Zero-based integer bin indices; -1 marks excluded positions.
#' @export
assignCMBins <- function(cM, cm_min, bin_width_cm, n_bins) {
    checkInputs(n_bins = scalarInteger)
    n_bins <- as.integer(n_bins)
    if (!is.numeric(cM) || any(!is.finite(cM)) ||
        length(cm_min) != 1L || !is.finite(cm_min) ||
        length(bin_width_cm) != 1L || !is.finite(bin_width_cm)) {
        stop("positions and bin limits must be finite", call. = FALSE)
    }
    out <- rep(-1L, length(cM))
    if (bin_width_cm <= 0) return(out)
    for (i in seq_along(cM)) {
        d <- cM[i] - cm_min
        if (d < 0) next
        b <- floor(d / bin_width_cm)
        out[i] <- if (b >= 0 && b < n_bins) as.integer(b) else -1L
    }
    out
}

#' Count marker pairs at each bin separation
#' @param snps_per_bin Nonnegative integer counts per bin.
#' @return Numeric vector of counts at lags zero through B - 1.
#' @export
pairCountHistogram <- function(snps_per_bin) {
    if (!is.numeric(snps_per_bin) || !length(snps_per_bin) ||
        any(!is.finite(snps_per_bin)) || any(snps_per_bin < 0) ||
        any(snps_per_bin != floor(snps_per_bin))) {
        stop("snps_per_bin must contain nonnegative integer counts", call. = FALSE)
    }
    B <- length(snps_per_bin)
    hist <- numeric(B)
    for (k in seq.int(0L, B - 1L)) {
        h <- 0
        for (b in seq_len(B - k)) {
            h <- h + as.double(snps_per_bin[b]) * as.double(snps_per_bin[b + k])
        }
        hist[k + 1L] <- h
    }
    hist
}

#' Accumulate weighted allele counts within genetic-distance bins
#' @param genotype Numeric matrix with individuals in rows and markers in columns.
#' @param p_A,p_B Source-population allele frequencies, one per marker.
#' @param p_hat Cohort allele frequencies, one per marker.
#' @param bin_id Zero-based bin indices, or -1 to exclude a marker.
#' @param n_bins Number of bins.
#' @return List containing weighted sums, bin counts, squared sums and pair counts.
#' @export
weightedLDBins <- function(genotype, p_A, p_B, p_hat, bin_id, n_bins) {
    checkInputs(n_bins = scalarInteger)
    n_bins <- as.integer(n_bins)
    if (!is.matrix(genotype) || !is.numeric(genotype) ||
        any(!is.finite(genotype)) || nrow(genotype) < 1L || ncol(genotype) < 1L) {
        stop("genotype must be a finite nonempty numeric matrix", call. = FALSE)
    }
    N <- nrow(genotype)
    L <- ncol(genotype)
    if (any(lengths(list(p_A, p_B, p_hat, bin_id)) != L) ||
        any(!is.finite(c(p_A, p_B, p_hat, bin_id))) ||
        any(c(p_A, p_B, p_hat) < 0 | c(p_A, p_B, p_hat) > 1) ||
        any(bin_id != floor(bin_id))) {
        stop("use one valid frequency and integer bin index per marker", call. = FALSE)
    }
    S <- matrix(0, N, n_bins)
    snps_per_bin <- integer(n_bins)
    SS <- S_col_sum <- numeric(n_bins)
    contrast <- p_A - p_B
    for (i in seq_len(L)) {
        b <- bin_id[i]
        if (b < 0 || b >= n_bins) next
        b <- b + 1L
        snps_per_bin[b] <- snps_per_bin[b] + 1L
        cc <- contrast[i]
        mm <- 2 * p_hat[i]
        if (cc == 0) next
        for (s in seq_len(N)) {
            y <- cc * (genotype[s, i] - mm)
            S[s, b] <- S[s, b] + y
        }
    }
    for (b in seq_len(n_bins)) {
        col_sum <- col_sq <- 0
        for (s in seq_len(N)) {
            v <- S[s, b]
            col_sum <- col_sum + v
            col_sq <- col_sq + v * v
        }
        S_col_sum[b] <- col_sum
        SS[b] <- col_sq
    }
    list(S = S, snps_per_bin = snps_per_bin, SS = SS, S_col_sum = S_col_sum,
         pair_count = pairCountHistogram(snps_per_bin))
}

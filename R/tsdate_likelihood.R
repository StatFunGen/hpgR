# R translations of tskitR/src/tsdate_discrete.cpp and
# tsdate_rescaling.cpp (Daniel Nachun, GPL >= 3). Source function names and
# hashes are recorded in inst/provenance/source-map.tsv. Base R dpois replaces
# the identical Poisson log-PMF; triangular indexing remains zero based at APIs.

tdTsdateLogsumexp <- function(x) {
    a <- max(x)
    if (!is.finite(a)) return(a)
    a + log(sum(exp(x - a)))
}
tdTsdateMutationLikFixed <- function(muts, span, timediff, mutation_rate,
                                         log_space, standardize) {
    y <- stats::dpois(muts, timediff * mutation_rate * span, log = log_space)
    if (standardize && is.finite(max(y)))
        y <- if (log_space) y - max(y) else y / max(y)
    y
}
tdTsdatePrecalcMutationLiks <- function(mut_edges, edge_spans, is_fixed_child,
                                            timediff_lower_tri, mutation_rate,
                                            log_space, standardize) {
    index <- rep(-1L, length(mut_edges))
    keys <- data.frame(muts = integer(), span = numeric())
    values <- list()
    for (i in which(!is_fixed_child)) {
        k <- which(keys$muts == mut_edges[i] & keys$span == edge_spans[i])
        if (!length(k)) {
            keys <- rbind(keys, data.frame(muts = mut_edges[i], span = edge_spans[i]))
            k <- nrow(keys)
            values[[k]] <- tdTsdateMutationLikFixed(mut_edges[i], edge_spans[i],
                timediff_lower_tri, mutation_rate, log_space, standardize)
        }
        index[i] <- k - 1L
    }
    cache <- if (length(values)) do.call(rbind, values) else
        matrix(numeric(), 0L, length(timediff_lower_tri))
    list(cache = cache, edge_cache_idx = index, n_unique = as.integer(nrow(keys)))
}
tdTriangularSum <- function(input, starts, log_space) {
    ends <- c(starts[-1L], length(input))
    vapply(seq_along(starts), function(i) {
        value <- input[seq.int(starts[i] + 1L, ends[i])]
        if (log_space) tdTsdateLogsumexp(value) else sum(value)
    }, 0)
}
tdTsdateRowsumLowerTriLog <- function(input, row_starts)
    tdTriangularSum(input, row_starts, TRUE)
tdTsdateRowsumLowerTriLin <- function(input, row_starts)
    tdTriangularSum(input, row_starts, FALSE)
tdTsdateRowsumUpperTriLog <- function(input, col_starts)
    tdTriangularSum(input, col_starts, TRUE)
tdTsdateRowsumUpperTriLin <- function(input, col_starts)
    tdTriangularSum(input, col_starts, FALSE)
tdConstrainAges <- function(time, fixed, parent, child, epsilon, iterations = 0L) {
    # tsdate_constrain_ages, including the final rootward constraint sweep.
    parent <- parent + 1L
    child <- child + 1L
    cavity <- matrix(0, length(parent), 2L)
    if (iterations > 0L) for (iteration in seq_len(iterations)) {
        if (all(time[parent] - time[child] > epsilon)) break
        for (e in seq_along(parent)) {
            p <- parent[e]; c <- child[e]
            time[c] <- time[c] - cavity[e, 1L]
            time[p] <- time[p] - cavity[e, 2L]
            adjustment <- time[c] - time[p]
            cavity[e, ] <- 0
            if (adjustment > 0) {
                if (!fixed[c] && !fixed[p]) cavity[e, ] <- c(-adjustment, adjustment) / 2
                else if (fixed[c] && !fixed[p]) cavity[e, 2L] <- adjustment
                else if (!fixed[c] && fixed[p]) cavity[e, 1L] <- -adjustment
            }
            time[c] <- time[c] + cavity[e, 1L]
            time[p] <- time[p] + cavity[e, 2L]
        }
    }
    for (e in seq_along(parent)) if (time[child[e]] + epsilon >= time[parent[e]])
        time[parent[e]] <- time[child[e]] + epsilon
    time
}

# Validation scope is recorded in inst/VALIDATION.md.
# Adapted from Daniel Nachun, tskitR/R/msprime-rate-map.R.
# Numeric calculations are retained; S4 slots become a small S3 list.
# See inst/provenance/source-map.tsv for source hashes and validation scope.

#' Define rates along a sequence
#'
#' Store a constant rate within each interval. Physical positions are in bases;
#' a recombination rate is per base per generation. Its integral is in Morgans.
#' @param position Increasing breakpoints beginning at zero.
#' @param rate Nonnegative finite rates, one per interval.
#' @return A rate map represented by a list of positions and rates.
#' @export
rateMap <- function(position, rate) {
    position <- as.numeric(position)
    rate <- as.numeric(rate)
    if (length(position) < 2L || any(!is.finite(position)) ||
        position[1L] != 0 || any(diff(position) <= 0) ||
        length(rate) != length(position) - 1L ||
        any(!is.finite(rate)) || any(rate < 0)) {
        stop("use increasing positions from zero and one nonnegative rate per interval",
             call. = FALSE)
    }
    structure(list(position = position, rate = rate), class = "hpgRateMap")
}

#' @rdname rateMap
#' @param x A map from rateMap().
#' @param pos Positions at which to query the map.
#' @export
rateAt <- function(x, pos) {
    mapPositions(x, pos)
    idx <- findInterval(pos, x$position, rightmost.closed = TRUE)
    idx <- pmax(pmin(idx, length(x$rate)), 1L)
    x$rate[idx]
}

#' @rdname rateMap
#' @export
cumulativeRate <- function(x, pos) {
    mapPositions(x, pos)
    widths <- diff(x$position)
    cumulative <- c(0, cumsum(widths * x$rate))
    vapply(pos, function(xi) {
        i <- findInterval(xi, x$position, rightmost.closed = TRUE)
        i <- max(min(i, length(x$rate)), 1L)
        cumulative[i] + (xi - x$position[i]) * x$rate[i]
    }, numeric(1L))
}

#' @rdname rateMap
#' @export
meanRate <- function(x) {
    mapPositions(x, numeric())
    widths <- diff(x$position)
    sum(widths * x$rate) / utils::tail(x$position, 1L)
}

mapPositions <- function(x, pos) {
    if (!inherits(x, "hpgRateMap") || !is.numeric(pos) || any(!is.finite(pos)) ||
        any(pos < 0 | pos > utils::tail(x$position, 1L))) {
        stop("pos must lie within the supplied rate map", call. = FALSE)
    }
    invisible(NULL)
}

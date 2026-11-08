# From demographyTools: ibd_hashing::ExtendHash and SeedHash (Daniel Nachun).
# GPL (>= 3). Fixed-word matching, with skip_ratio = max_seeds = 0.

#' Find long matching segments in phased haplotypes
#'
#' Extend pairs with identical marker words and close a match after more than
#' `gap` unmatched words. Marker selection is supplied by the caller. Long
#' identity by state supports IBD only when marker density and error rates make
#' chance matches sufficiently rare.
#' @param haplotypes Complete 0/1 matrix, haplotypes in rows and markers in columns.
#' @param geneticPositions Marker positions in Morgans, increasing within chromosomes.
#' @param chromosome Chromosome label for each marker.
#' @param wordSize Number of markers in a matching word.
#' @param minimumCm Minimum reported segment length in centimorgans.
#' @param gap Number of unmatched words allowed inside a segment.
#' @return Data frame with haplotype indices, chromosome, marker endpoints,
#'   genetic endpoints in cM, and length in cM. Marker endpoints index the input.
#' @export
detectIbdSegments <- function(haplotypes, geneticPositions, chromosome = NULL,
                              wordSize = 1L, minimumCm = 8, gap = 0L) {
    x <- populationHaplotypes(haplotypes)
    checkInputs(geneticPositions = list(numericVector, length = ncol(x)),
                wordSize = list(scalarInteger, maximum = 64L),
                minimumCm = positiveScalar,
                gap = list(scalarInteger, minimum = 0L))
    if (is.null(chromosome)) chromosome <- rep("1", ncol(x))
    if (length(chromosome) != ncol(x) || anyNA(chromosome))
        stop("chromosome must label every marker", call. = FALSE)
    pairs <- utils::combn(nrow(x), 2L)
    chunks <- list()
    emit <- function(take, start, last, index, id) {
        if (!length(take)) return(NULL)
        first <- index[(start[take] - 1L) * wordSize + 1L]
        final <- index[pmin(last[take] * wordSize, length(index))]
        lengthCm <- 100 * (geneticPositions[final] - geneticPositions[first])
        keep <- lengthCm >= minimumCm
        data.frame(hap1 = pairs[1L, take[keep]], hap2 = pairs[2L, take[keep]],
            chromosome = rep(as.character(id), sum(keep)), fromSite = first[keep],
            toSite = final[keep], startCm = 100 * geneticPositions[first[keep]],
            endCm = 100 * geneticPositions[final[keep]], lengthCm = lengthCm[keep])
    }
    for (id in unique(chromosome)) {
        index <- which(chromosome == id)
        if (is.unsorted(geneticPositions[index], strictly = TRUE))
            stop("positions must increase within each chromosome", call. = FALSE)
        start <- last <- integer(ncol(pairs))
        for (word in seq_len(ceiling(length(index) / wordSize))) {
            cols <- index[seq.int((word - 1L) * wordSize + 1L,
                                 min(word * wordSize, length(index)))]
            signature <- if (length(cols) == 1L) x[, cols] else
                apply(x[, cols, drop = FALSE], 1L, paste0, collapse = "")
            same <- signature[pairs[1L, ]] == signature[pairs[2L, ]]
            start[same & start == 0L] <- word
            last[same] <- word
            expired <- which(start > 0L & last < word - gap)
            if (length(expired)) {
                chunks[[length(chunks) + 1L]] <- emit(expired, start, last, index, id)
                start[expired] <- last[expired] <- 0L
            }
        }
        chunks[[length(chunks) + 1L]] <- emit(which(start > 0L), start, last, index, id)
    }
    result <- do.call(rbind, chunks)
    if (is.null(result)) result <- data.frame(hap1 = integer(), hap2 = integer(),
        chromosome = character(), fromSite = integer(), toSite = integer(),
        startCm = numeric(), endCm = numeric(), lengthCm = numeric())
    rownames(result) <- NULL
    result
}

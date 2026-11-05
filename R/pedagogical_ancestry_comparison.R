# Classroom comparisons of inferred and recorded R node/edge tables.
# Recorded ancestry supplies the truth for evaluating inference results.

ancestryClades <- function(x, position) {
    parent <- ancestryParent(x, position)
    paths <- lapply(x$sampleIds, ancestryPath, parent = parent)
    nodes <- unique(unlist(paths))
    clades <- vapply(nodes, function(u) {
        descendants <- which(vapply(paths, function(p) u %in% p, logical(1)))
        if (length(descendants) < 2L || length(descendants) == length(paths)) return("")
        paste(descendants, collapse = ",")
    }, "")
    unique(clades[nzchar(clades)])
}

#' Compare inferred and recorded local ancestry
#'
#' Evaluate clade agreement and ancestor-age error over the union of genomic
#' intervals. Sample-copy order must correspond between the two objects.
#' Unary records and internal node identifiers do not affect clade comparisons.
#' @param estimate An inferred or dated hpgR_ancestry object.
#' @param truth Recorded hpgR_ancestry with the same ordered sampled copies.
#' @return Interval-level clade precision/recall, pairwise MRCA and root-age
#'   absolute errors, plus length-weighted summaries. Clades exclude tips and
#'   the full-sample root. Pair errors average over all sampled-copy pairs.
#' @export
compareAncestry <- function(estimate, truth) {
    ancestryCheckTables(estimate)
    ancestryCheckTables(truth)
    if (length(estimate$sampleIds) != length(truth$sampleIds) ||
        estimate$sequenceLength != truth$sequenceLength)
        stop("ancestry objects must describe the same ordered sample and region", call. = FALSE)
    cuts <- sort(unique(c(0, truth$sequenceLength, estimate$edges$left,
        estimate$edges$right, truth$edges$left, truth$edges$right)))
    pairs <- utils::combn(seq_along(truth$sampleIds), 2L)
    ages <- function(x, position) {
        parent <- ancestryParent(x, position)
        paths <- lapply(x$sampleIds, ancestryPath, parent = parent)
        time <- stats::setNames(x$nodes$time, x$nodes$id)
        vapply(seq_len(ncol(pairs)), function(j) {
            a <- paths[[pairs[1L, j]]]; b <- paths[[pairs[2L, j]]]
            common <- a[a %in% b]
            if (length(common)) unname(time[as.character(common[1L])]) else NaN
        }, 0)
    }
    result <- lapply(seq_len(length(cuts) - 1L), function(j) {
        pos <- mean(cuts[c(j, j + 1L)])
        a <- ancestryClades(estimate, pos)
        b <- ancestryClades(truth, pos)
        shared <- length(intersect(a, b))
        timesA <- ages(estimate, pos)
        timesB <- ages(truth, pos)
        data.frame(left = cuts[j], right = cuts[j + 1L],
            cladePrecision = if (length(a)) shared / length(a) else 1,
            cladeRecall = if (length(b)) shared / length(b) else 1,
            pairMAE = mean(abs(timesA - timesB)),
            rootMAE = abs(max(timesA) - max(timesB)))
    })
    intervals <- do.call(rbind, result)
    summaries <- vapply(intervals[c("cladePrecision", "cladeRecall", "pairMAE", "rootMAE")],
        stats::weighted.mean, 0, w = intervals$right - intervals$left)
    list(intervals = intervals, lengthWeighted = summaries)
}

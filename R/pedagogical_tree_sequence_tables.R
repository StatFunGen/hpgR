# Classroom R representation, traversal and display of ancestry tables.
# Half-open edge intervals and zero-based node IDs follow the source table
# conventions. Display trees suppress unary nodes; recorded tables retain
# their mutation-bearing branches. These helpers use explicit R traversals.

.ancestryCheckTables <- function(x) {
    if (!inherits(x, "hpgR_ancestry") || !is.data.frame(x$nodes) ||
        !all(c("id", "time", "flags") %in% names(x$nodes)) ||
        !is.data.frame(x$edges) ||
        !all(c("left", "right", "parent", "child") %in% names(x$edges))) {
        stop("ancestry must be an hpgR_ancestry object with R node and edge tables",
             call. = FALSE)
    }
    nodes <- x$nodes
    edges <- x$edges
    if (!nrow(nodes) || anyDuplicated(nodes$id) ||
        any(!is.finite(nodes$time) | nodes$time < 0) ||
        any(nodes$id != floor(nodes$id) | nodes$id < 0) ||
        !is.numeric(x$sampleIds) || length(x$sampleIds) < 2L ||
        anyDuplicated(x$sampleIds) || any(!x$sampleIds %in% nodes$id) ||
        length(x$sequenceLength) != 1L || !is.finite(x$sequenceLength) ||
        x$sequenceLength <= 0 || any(!is.finite(as.matrix(edges))) ||
        any(edges$left < 0 | edges$right > x$sequenceLength |
            edges$left >= edges$right) ||
        any(!edges$parent %in% nodes$id | !edges$child %in% nodes$id) ||
        any(nodes$time[match(edges$parent, nodes$id)] <=
            nodes$time[match(edges$child, nodes$id)])) {
        stop("invalid ancestry tables, sample IDs, intervals or node ages", call. = FALSE)
    }
    invisible(x)
}

.ancestryParent <- function(ancestry, position) {
    edges <- ancestry$edges
    rows <- which(edges$left <= position & position < edges$right)
    if (anyDuplicated(edges$child[rows])) {
        stop("a child has multiple parents at the same genomic position", call. = FALSE)
    }
    parent <- stats::setNames(rep(-1L, nrow(ancestry$nodes)), ancestry$nodes$id)
    parent[as.character(edges$child[rows])] <- edges$parent[rows]
    parent
}

.ancestryPath <- function(node, parent) {
    path <- integer()
    while (node != -1L) {
        if (node %in% path || !as.character(node) %in% names(parent)) {
            stop("invalid local ancestry path", call. = FALSE)
        }
        path <- c(path, as.integer(node))
        node <- unname(parent[as.character(node)])
    }
    path
}

.ancestryMarginal <- function(ancestry, position) {
    parent <- .ancestryParent(ancestry, position)
    samples <- ancestry$sampleIds
    paths <- lapply(samples, .ancestryPath, parent = parent)
    present <- unique(unlist(paths, use.names = FALSE))
    childCounts <- tabulate(match(unname(parent[as.character(present)]), present),
                            nbins = length(present))
    keep <- present[present %in% samples | childCounts >= 2L]
    if (!length(keep)) stop("no local sample ancestry", call. = FALSE)
    list(parent = parent, paths = paths, keep = keep)
}

.ancestryLocal <- function(ancestry, left, right) {
    marginal <- .ancestryMarginal(ancestry, (left + right) / 2)
    parent <- marginal$parent
    paths <- marginal$paths
    keep <- marginal$keep
    samples <- ancestry$sampleIds
    displayParent <- vapply(keep, function(node) {
        path <- .ancestryPath(node, parent)[-1L]
        above <- path[path %in% keep]
        if (length(above)) above[1L] else -1L
    }, integer(1))
    pairPaths <- paths[1:2]
    shared <- pairPaths[[1L]][pairPaths[[1L]] %in% pairPaths[[2L]]]
    pairNode <- if (length(shared)) shared[1L] else NA_integer_
    roots <- keep[displayParent == -1L]
    ages <- stats::setNames(ancestry$nodes$time, ancestry$nodes$id)
    descends <- function(sample, node) node %in% paths[[match(sample, samples)]]
    # Order tips by a depth-first traversal, so each clade is contiguous.
    tipOrder <- integer()
    visit <- function(node) {
        if (node %in% samples) tipOrder <<- c(tipOrder, node)
        for (child in keep[displayParent == node]) visit(child)
    }
    for (root in roots) visit(root)
    x <- vapply(keep, function(node) {
        children <- tipOrder[vapply(tipOrder, descends, logical(1), node = node)]
        mean(range(match(children, tipOrder)))
    }, numeric(1))
    tree <- data.frame(node = as.integer(keep), parent = displayParent,
        time = unname(ages[as.character(keep)]), x = x,
        is_sample = keep %in% samples, row.names = NULL)
    list(tree = tree, root_node = if (length(roots) == 1L) as.integer(roots) else NA_integer_,
        pair_node = as.integer(pairNode), tmrca = if (length(roots) == 1L)
        unname(ages[as.character(roots)]) else NaN,
        pair_time = if (!is.na(pairNode)) unname(ages[as.character(pairNode)]) else NaN)
}

.ancestryExtract <- function(ancestry) {
    .ancestryCheckTables(ancestry)
    edges <- ancestry$edges
    cuts <- sort(unique(c(0, ancestry$sequenceLength, edges$left, edges$right)))
    local <- lapply(seq_len(length(cuts) - 1L), function(j)
        .ancestryLocal(ancestry, cuts[j], cuts[j + 1L]))
    # Distinct recorded boundaries can describe the same unary-suppressed tree.
    signatures <- vapply(local, function(x)
        paste(x$tree$node, x$tree$parent, sep = ":", collapse = ","), character(1))
    start <- c(1L, which(signatures[-1L] != head(signatures, -1L)) + 1L)
    end <- c(start[-1L] - 1L, length(local))
    intervals <- data.frame(left = cuts[start], right = cuts[end + 1L],
        tmrca = vapply(local[start], `[[`, numeric(1), "tmrca"),
        pair_time = vapply(local[start], `[[`, numeric(1), "pair_time"),
        root_node = vapply(local[start], `[[`, integer(1), "root_node"),
        pair_node = vapply(local[start], `[[`, integer(1), "pair_node"))
    positions <- if (is.null(ancestry$sites)) numeric() else ancestry$sites$position
    haplotypes <- matrix(0L, nrow = length(ancestry$sampleIds), ncol = length(positions))
    if (length(positions)) {
        if (anyDuplicated(positions) || is.null(ancestry$mutations)) {
            stop("extraction requires unique site positions and a mutation table",
                 call. = FALSE)
        }
        for (j in seq_along(positions)) {
            parent <- .ancestryParent(ancestry, positions[j])
            rows <- which(ancestry$mutations$site == ancestry$sites$id[j])
            origin <- ancestry$mutations$node[rows]
            state <- if ("derived_state" %in% names(ancestry$mutations))
                as.integer(ancestry$mutations$derived_state[rows]) else rep(1L, length(rows))
            rootState <- if ("ancestral_state" %in% names(ancestry$sites))
                as.integer(ancestry$sites$ancestral_state[j]) else 0L
            haplotypes[, j] <- vapply(ancestry$sampleIds, function(sample) {
                path <- .ancestryPath(sample, parent)
                affected <- which(origin %in% path)
                if (!length(affected)) return(rootState)
                depth <- match(origin[affected], path)
                nearest <- tail(affected[depth == min(depth)], 1L)
                state[nearest]
            }, integer(1))
        }
    }
    list(nHaplotypes = length(ancestry$sampleIds), nSites = length(positions),
        nTrees = length(start), nRecordedTrees = length(local),
        nNodes = nrow(ancestry$nodes), nEdges = nrow(edges),
        timeUnits = ancestry$timeUnits,
        recombinationNodes = sum(bitwAnd(as.integer(ancestry$nodes$flags),
                                         bitwShiftL(1L, 17L)) != 0L),
        observed = list(haplotypes = haplotypes, positions = positions,
                        sequenceLength = ancestry$sequenceLength),
        intervals = intervals, localTrees = lapply(local[start], `[[`, "tree"),
        nodes = ancestry$nodes, edges = edges)
}

#' Inspect R ancestry tables and observed sequences
#' @param ancestry An hpgR_ancestry R object.
#' @return Observed phased haplotypes, ordered positions and analyzed length;
#'   local root/pair ages; full node/edge tables and unary-suppressed display
#'   trees. Node/sample IDs remain zero based. This is R table traversal, not
#'   a port of the complete tskit simplifier or an ancestry-inference method.
#' @export
ancestryData <- function(ancestry) .ancestryExtract(ancestry)

# Convert copying paths and reduce marginal tables for classroom inspection.
.tsiEdges <- function(match, child, positions, sequenceLength) {
    parent <- match$parent
    a <- c(1L, which(diff(parent) != 0L) + 1L)
    b <- c(a[-1L] - 1L, length(parent))
    indices <- seq.int(match$start, match$end)
    boundary <- c(0, positions[-1L], sequenceLength)
    data.frame(left = boundary[indices[a]], right = boundary[indices[b] + 1L],
               parent = parent[a], child = rep(as.integer(child), length(a)))
}

.tsiSquashEdges <- function(edges) {
    if (!nrow(edges)) return(edges)
    edges <- edges[order(edges$parent, edges$child, edges$left, edges$right), ]
    result <- edges[1L, , drop = FALSE]
    if (nrow(edges) > 1L) for (j in 2:nrow(edges)) {
        k <- nrow(result)
        if (edges$parent[j] == result$parent[k] &&
            edges$child[j] == result$child[k] && edges$left[j] == result$right[k]) {
            result$right[k] <- edges$right[j]
        } else result <- rbind(result, edges[j, , drop = FALSE])
    }
    rownames(result) <- NULL
    result
}

.tsiSimplify <- function(x) {
    # Traverse each marginal tree, retain sample and branching nodes, and move
    # each mutation to its surviving descendant branch. Equivalent tables use
    # the same ancestral node ID across intervals. This adapter targets R lists.
    cuts <- sort(unique(c(0, x$sequenceLength, x$edges$left, x$edges$right)))
    ed <- list()
    mut <- list()
    for (j in seq_len(length(cuts) - 1L)) {
        left <- cuts[j]
        right <- cuts[j + 1L]
        marginal <- .ancestryMarginal(x, (left + right) / 2)
        parent <- marginal$parent
        keep <- marginal$keep
        for (u in keep) {
            above <- .ancestryPath(u, parent)[-1L]
            p <- above[above %in% keep]
            if (length(p)) ed[[length(ed) + 1L]] <- data.frame(
                left = left, right = right, parent = p[1L], child = u)
        }
        at <- which(x$sites$position >= left & x$sites$position < right)
        for (s in at) {
            rows <- which(x$mutations$site == x$sites$id[s])
            # Depth order retains recurrent and back-mutation ordering.
            rows <- rows[order(x$nodes$time[match(x$mutations$node[rows], x$nodes$id)],
                               decreasing = TRUE)]
            for (r in rows) {
                origin <- x$mutations$node[r]
                under <- keep[vapply(keep, function(u)
                    origin %in% .ancestryPath(u, parent), logical(1))]
                if (!length(under)) next
                v <- under[which.max(x$nodes$time[match(under, x$nodes$id)])]
                mut[[length(mut) + 1L]] <- data.frame(site = x$mutations$site[r],
                    node = v, derived_state = x$mutations$derived_state[r])
            }
        }
    }
    x$edges <- .tsiSquashEdges(do.call(rbind, ed))
    keep <- c(x$sampleIds, sort(setdiff(unique(c(x$edges$parent, x$edges$child)),
                                      x$sampleIds)))
    x$nodes <- x$nodes[match(keep, x$nodes$id), , drop = FALSE]
    x$mutations <- if (length(mut)) do.call(rbind, mut) else
        data.frame(site = integer(), node = integer(), derived_state = integer())
    x$edges$parent <- as.integer(match(x$edges$parent, keep) - 1L)
    x$edges$child <- as.integer(match(x$edges$child, keep) - 1L)
    x$mutations$node <- as.integer(match(x$mutations$node, keep) - 1L)
    x$nodes$id <- seq_along(keep) - 1L
    x$sampleIds <- seq_along(x$sampleIds) - 1L
    rownames(x$nodes) <- rownames(x$mutations) <- NULL
    x
}

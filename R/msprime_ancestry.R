# From msprime 1.4.1: msp_run_coalescent, msp_std_common_ancestor_event,
# msp_merge_two_ancestors and msp_recombination_event.
# Copyright (C) 2015-2025 University of Oxford and msprime contributors.
# GPL-3.0-or-later; source supplied in Daniel Nachun's tskitR/src/msprime.

msprimeEmptySegments <- function() {
    matrix(numeric(), nrow = 0L, ncol = 3L,
           dimnames = list(NULL, c("left", "right", "node")))
}

msprimeCompressSegments <- function(x) {
    if (nrow(x) < 2L) return(x)
    out <- x[1L, , drop = FALSE]
    for (j in seq.int(2L, nrow(x))) {
        last <- nrow(out)
        if (out[last, "right"] == x[j, "left"] &&
            out[last, "node"] == x[j, "node"]) {
            out[last, "right"] <- x[j, "right"]
        } else out <- rbind(out, x[j, , drop = FALSE])
    }
    rownames(out) <- NULL
    out
}

msprimeStoreNode <- function(state, flags = 0L) {
    id <- length(state$time)
    state$time <- c(state$time, state$currentTime)
    state$flags <- c(state$flags, as.integer(flags))
    as.integer(id)
}

msprimeStoreEdge <- function(state, left, right, parent, child) {
    if (left >= right || parent == child ||
        state$time[parent + 1L] <= state$time[child + 1L]) {
        stop("Invalid Hudson edge: positive span and increasing parent time required",
             call. = FALSE)
    }
    state$edges[[length(state$edges) + 1L]] <- c(left, right, parent, child)
    invisible(NULL)
}

msprimeStoreArgEdges <- function(state, segments, node) {
    for (j in seq_len(nrow(segments))) {
        if (segments[j, "node"] != node) {
            msprimeStoreEdge(state, segments[j, "left"], segments[j, "right"],
                              node, segments[j, "node"])
            segments[j, "node"] <- node
        }
    }
    msprimeCompressSegments(segments)
}

msprimeSplitOverlap <- function(state, points) {
    old <- state$overlap
    breaks <- sort(unique(c(old[, "left"], tail(old[, "right"], 1L), points)))
    left <- head(breaks, -1L)
    idx <- findInterval(left, old[, "left"])
    state$overlap <- cbind(left = left, right = tail(breaks, -1L),
                          count = old[idx, "count"])
    invisible(NULL)
}

msprimeCompressOverlap <- function(state) {
    x <- state$overlap
    if (nrow(x) < 2L) return(invisible(NULL))
    begins <- c(1L, which(diff(x[, "count"]) != 0L) + 1L)
    ends <- c(begins[-1L] - 1L, nrow(x))
    state$overlap <- cbind(left = x[begins, "left"], right = x[ends, "right"],
                          count = x[begins, "count"])
    invisible(NULL)
}

msprimeCover <- function(segments, left) {
    idx <- findInterval(left, segments[, "left"])
    nodes <- rep(NA_integer_, length(left))
    possible <- which(idx > 0L)
    present <- possible[left[possible] < segments[idx[possible], "right"]]
    nodes[present] <- as.integer(segments[idx[present], "node"])
    nodes
}

# From msprime's msp_merge_two_ancestors and copy_overlap_count.
msprimeMergeAncestors <- function(state, first, second) {
    a <- state$lineages[[as.character(first)]]
    b <- state$lineages[[as.character(second)]]
    if (is.null(a) || is.null(b) || first == second) {
        stop("A common-ancestor event must select two distinct live lineages",
             call. = FALSE)
    }
    msprimeSplitOverlap(state, c(a[, "left"], a[, "right"],
                                   b[, "left"], b[, "right"]))
    lo <- min(a[1L, "left"], b[1L, "left"])
    hi <- max(tail(a[, "right"], 1L), tail(b[, "right"], 1L))
    take <- which(state$overlap[, "left"] >= lo &
                  state$overlap[, "right"] <= hi)
    left <- state$overlap[take, "left"]
    right <- state$overlap[take, "right"]
    nodeA <- msprimeCover(a, left)
    nodeB <- msprimeCover(b, left)
    both <- !is.na(nodeA) & !is.na(nodeB)
    coalescence <- any(both)
    newNode <- if (coalescence) msprimeStoreNode(state, 0L) else
        if (state$fullArg) msprimeStoreNode(state, 2^18) else NA_integer_
    result <- msprimeEmptySegments()
    for (j in seq_along(left)) {
        value <- NA_integer_
        if (both[j]) {
            count <- state$overlap[take[j], "count"]
            if (count < 2L) stop("Inconsistent Hudson overlap count", call. = FALSE)
            msprimeStoreEdge(state, left[j], right[j], newNode, nodeA[j])
            msprimeStoreEdge(state, left[j], right[j], newNode, nodeB[j])
            if (count == 2L) {
                state$overlap[take[j], "count"] <- 0L
            } else {
                state$overlap[take[j], "count"] <- count - 1L
                value <- newNode
            }
        } else if (!is.na(nodeA[j])) value <- nodeA[j] else
            if (!is.na(nodeB[j])) value <- nodeB[j]
        if (!is.na(value)) {
            result <- rbind(result, c(left = left[j], right = right[j], node = value))
        }
    }
    if (state$fullArg && nrow(result)) {
        result <- msprimeStoreArgEdges(state, result, newNode)
    } else result <- msprimeCompressSegments(result)
    msprimeCompressOverlap(state)
    state$lineages[[as.character(first)]] <- NULL
    state$lineages[[as.character(second)]] <- NULL
    id <- state$nextLineage
    state$nextLineage <- id + 1L
    if (nrow(result)) state$lineages[[as.character(id)]] <- result
    state$commonAncestorEvents <- state$commonAncestorEvents + 1L
    if (!coalescence) state$nonoverlapEvents <- state$nonoverlapEvents + 1L
    invisible(id)
}

# msp_recombination_event, including a breakpoint in a gap between two carried
# segments. The left lineage persists; the right chain receives a new lineage.
msprimeRecombination <- function(state, lineage, breakpoint) {
    x <- state$lineages[[as.character(lineage)]]
    if (is.null(x) || length(breakpoint) != 1L || !is.finite(breakpoint) ||
        breakpoint != floor(breakpoint) || breakpoint <= x[1L, "left"] ||
        breakpoint >= tail(x[, "right"], 1L)) {
        stop("A recombination cut must be an interior integer of a live lineage span",
             call. = FALSE)
    }
    left <- x[x[, "left"] < breakpoint, , drop = FALSE]
    right <- x[x[, "right"] > breakpoint, , drop = FALSE]
    left[, "right"] <- pmin(left[, "right"], breakpoint)
    right[, "left"] <- pmax(right[, "left"], breakpoint)
    inSegment <- any(x[, "left"] < breakpoint & x[, "right"] > breakpoint)
    if (state$fullArg) {
        leftNode <- msprimeStoreNode(state, 2^17)
        left <- msprimeStoreArgEdges(state, left, leftNode)
        rightNode <- msprimeStoreNode(state, 2^17)
        right <- msprimeStoreArgEdges(state, right, rightNode)
    }
    id <- state$nextLineage
    state$nextLineage <- id + 1L
    state$lineages[[as.character(lineage)]] <- left
    state$lineages[[as.character(id)]] <- right
    state$recombinationEvents <- state$recombinationEvents + 1L
    if (!inSegment) state$gapRecombinations <- state$gapRecombinations + 1L
    invisible(id)
}

msprimeFinishEdges <- function(state) {
    if (!length(state$edges)) {
        return(data.frame(left = numeric(), right = numeric(),
                          parent = integer(), child = integer()))
    }
    x <- as.data.frame(do.call(rbind, state$edges))
    names(x) <- c("left", "right", "parent", "child")
    x$parent <- as.integer(x$parent)
    x$child <- as.integer(x$child)
    x <- x[order(x$parent, x$child, x$left, x$right), , drop = FALSE]
    out <- x[1L, , drop = FALSE]
    if (nrow(x) > 1L) for (j in seq.int(2L, nrow(x))) {
        k <- nrow(out)
        if (out$parent[k] == x$parent[j] && out$child[k] == x$child[j] &&
            out$right[k] == x$left[j]) {
            out$right[k] <- x$right[j]
        } else out <- rbind(out, x[j, , drop = FALSE])
    }
    out <- out[order(state$time[out$parent + 1L], out$parent,
                     out$child, out$left), , drop = FALSE]
    rownames(out) <- NULL
    out
}

msprimeAncestry <- function(samples, populationSize, sequenceLength,
                              recombinationRate = 0, seed = 8149L,
                              recordFullArg = FALSE, events = NULL) {
    checkInputs(samples = scalarInteger)
    samples <- as.integer(samples)
    if (samples > .Machine$integer.max / 2L) stop("Too many sample copies", call. = FALSE)
    msprimeAncestryCopies(2L * samples, populationSize, sequenceLength,
                           recombinationRate, seed, recordFullArg, events)
}

# The source event kernel operates on sampled copies. Diploid sample grouping
# is an interface convention; an odd number of sampled copies uses this same
# kernel with the same diploid Ne and pair-merger rate 1/(2Ne).
msprimeAncestryCopies <- function(copies, populationSize, sequenceLength,
                                    recombinationRate = 0, seed = 8149L,
                                    recordFullArg = FALSE, events = NULL) {
    checkInputs(copies = list(scalarInteger, minimum = 2L), populationSize = positiveScalar,
                sequenceLength = scalarInteger,
                recombinationRate = list(positiveScalar, zero = TRUE),
                seed = list(scalarInteger, minimum = 0L), recordFullArg = logicalScalar)
    copies <- as.integer(copies)
    sequenceLength <- as.integer(sequenceLength)
    seed <- as.integer(seed)
    if (!is.null(events) && !is.list(events)) {
        stop("events must be a recorded Hudson event list", call. = FALSE)
    }
    simupopWithSeed(seed, function() {
        state <- new.env(parent = emptyenv())
        state$lineages <- lapply(seq_len(copies), function(j)
            matrix(c(0, sequenceLength, j - 1L), nrow = 1L,
                   dimnames = list(NULL, c("left", "right", "node"))))
        names(state$lineages) <- as.character(seq_len(copies) - 1L)
        state$nextLineage <- copies
        state$currentTime <- 0
        state$time <- rep(0, copies)
        state$flags <- rep(1L, copies)
        state$edges <- list()
        state$overlap <- matrix(c(0, sequenceLength, copies), nrow = 1L,
            dimnames = list(NULL, c("left", "right", "count")))
        state$fullArg <- recordFullArg
        state$commonAncestorEvents <- 0L
        state$nonoverlapEvents <- 0L
        state$recombinationEvents <- 0L
        state$gapRecombinations <- 0L
        recorded <- list()
        while (length(state$lineages)) {
            step <- length(recorded) + 1L
            if (step > 1000000L) {
                stop("Restricted R Hudson simulator exceeded one million events",
                     call. = FALSE)
            }
            ids <- as.integer(names(state$lineages))
            n <- length(ids)
            if (is.null(events)) {
                masses <- vapply(state$lineages, function(x)
                    recombinationRate *
                    (tail(x[, "right"], 1L) - x[1L, "left"] - 1), numeric(1L))
                totalRecombination <- sum(masses)
                caRate <- if (n > 1L) n * (n - 1L) / (4 * populationSize) else 0
                if (!is.finite(totalRecombination) || !is.finite(caRate)) {
                    stop("Hudson event rates exceed the finite numeric range",
                         call. = FALSE)
                }
                reWait <- if (totalRecombination > 0)
                    stats::rexp(1L, totalRecombination) else Inf
                caWait <- if (n > 1L)
                    stats::rexp(1L, caRate) else Inf
                wait <- min(reWait, caWait)
                if (!is.finite(wait)) {
                    stop("Hudson process has ancestral material but no possible event",
                         call. = FALSE)
                }
                time <- state$currentTime + wait
                if (time <= state$currentTime) time <- state$currentTime +
                    max(.Machine$double.xmin,
                        abs(state$currentTime) * .Machine$double.eps)
                if (reWait == wait) {
                    which <- sample.int(n, 1L, prob = masses)
                    x <- state$lineages[[which]]
                    cut <- x[1L, "left"] +
                        sample.int(tail(x[, "right"], 1L) - x[1L, "left"] - 1L, 1L)
                    event <- list(type = "recombination", time = time,
                                  lineage = ids[which], breakpoint = as.numeric(cut))
                } else event <- list(type = "common_ancestor", time = time,
                                      lineages = ids[sample.int(n, 2L, replace = FALSE)])
            } else {
                if (step > length(events)) {
                    stop("Hudson event replay ended before all local MRCAs were reached",
                         call. = FALSE)
                }
                event <- events[[step]]
            }
            if (!is.list(event) || !is.character(event$type) ||
                length(event$type) != 1L || !is.numeric(event$time) ||
                length(event$time) != 1L || !is.finite(event$time) ||
                event$time <= state$currentTime) {
                stop("Each Hudson event must have a type and increasing finite time",
                     call. = FALSE)
            }
            state$currentTime <- event$time
            if (event$type == "recombination") {
                if (recombinationRate == 0 || !is.numeric(event$lineage) ||
                    length(event$lineage) != 1L ||
                    !event$lineage %in% ids) {
                    stop("Recombination replay must select a live lineage at positive rate",
                         call. = FALSE)
                }
                msprimeRecombination(state, event$lineage, event$breakpoint)
            } else if (event$type == "common_ancestor") {
                if (!is.numeric(event$lineages) || length(event$lineages) != 2L ||
                    any(!event$lineages %in% ids) ||
                    event$lineages[1L] == event$lineages[2L]) {
                    stop("Common-ancestor replay must select two live lineages",
                         call. = FALSE)
                }
                msprimeMergeAncestors(state, event$lineages[1L], event$lineages[2L])
            } else stop("Unknown Hudson event type", call. = FALSE)
            recorded[[step]] <- event
        }
        if (!is.null(events) && length(events) != length(recorded)) {
            stop("Hudson event replay contains events after process completion",
                 call. = FALSE)
        }
        if (any(state$overlap[, "count"] != 0)) {
            stop("Hudson process completed with unresolved overlap counts", call. = FALSE)
        }
        structure(list(nodes = data.frame(id = seq_along(state$time) - 1L,
            time = state$time, flags = state$flags),
            edges = msprimeFinishEdges(state), sampleIds = seq_len(copies) - 1L,
            sequenceLength = as.numeric(sequenceLength), timeUnits = "generations",
            events = recorded,
            eventCounts = c(common_ancestor = state$commonAncestorEvents,
                recombination = state$recombinationEvents,
                nonoverlap_common_ancestor = state$nonoverlapEvents,
                gap_recombination = state$gapRecombinations),
            source = list(implementation = "restricted pure-R Hudson translation",
                reference = "msprime 1.4.1 vendored in tskitR",
                populationSize = populationSize, recombinationRate = recombinationRate,
                discreteGenome = TRUE, recordFullArg = recordFullArg,
                rng = "R Mersenne-Twister", seed = seed)),
            class = c("hpgR_ancestry", "list"))
    })
}

# Restricted translation of vendored msprime 1.4.1 mutation placement.
# Source in Daniel Nachun's tskitR snapshot:
# src/msprime/mutgen.c::mutgen_place_mutations (1077--1158),
# mutation_matrix_choose_root_state and mutation_matrix_transition.
# Copyright University of Oxford; GPL-3.0-or-later. See provenance records.
# Uniform continuous positions, constant rate, and BinaryMutationModel with
# root state 0 are supported. Unique positions mean each new site has one
# 0 -> 1 transition. Site and mutation records use R data frames. batchEdges
# draws all edge counts before positions and times for long discrete records.
# Both draw orders use the same Poisson exposure and uniform placement laws.

.msprimeMutations <- function(ancestry, mutationRate, seed, events = NULL,
                              batchEdges = FALSE) {
    .ancestryCheckTables(ancestry)
    mutationRate <- .positiveScalar(mutationRate, "mutationRate", zero = TRUE)
    seed <- .scalarInteger(seed, "seed", minimum = 0L)
    if (!is.logical(batchEdges) || length(batchEdges) != 1L || is.na(batchEdges))
        stop("batchEdges must be TRUE or FALSE", call. = FALSE)
    nodes <- ancestry$nodes
    edges <- ancestry$edges
    draw <- function() {
        low <- nodes$time[match(edges$child, nodes$id)]
        high <- nodes$time[match(edges$parent, nodes$id)]
        rates <- mutationRate * (edges$right - edges$left) * (high - low)
        if (any(!is.finite(rates))) stop("mutation exposure is too large", call. = FALSE)
        if (batchEdges && is.null(events)) {
            chosen <- rep(seq_len(nrow(edges)), stats::rpois(nrow(edges), rates))
            positions <- stats::runif(length(chosen), edges$left[chosen], edges$right[chosen])
            times <- low[chosen] + stats::runif(length(chosen)) * (high - low)[chosen]
            while (anyDuplicated(positions)) {
                repeated <- which(duplicated(positions))
                positions[repeated] <- stats::runif(length(repeated),
                    edges$left[chosen[repeated]], edges$right[chosen[repeated]])
            }
            events <- data.frame(edge = chosen, position = positions, time = times)
        }
        mutationRows <- vector("list", nrow(edges))
        used <- numeric()
        if (!is.null(events)) {
            required <- c("edge", "position", "time")
            if (!is.data.frame(events) || !all(required %in% names(events)) ||
                any(!is.finite(as.matrix(events[required]))) ||
                any(events$edge != floor(events$edge) | events$edge < 1L |
                    events$edge > nrow(edges)) || anyDuplicated(events$position)) {
                stop("events must give distinct continuous positions, times and one-based edge rows",
                     call. = FALSE)
            }
            j <- events$edge
            if (any(events$position < edges$left[j] | events$position >= edges$right[j]) ||
                any(events$time < low[j] | events$time >= high[j])) {
                stop("mutation events must lie on the specified edge and branch", call. = FALSE)
            }
            mutationRows <- list(if (length(j)) data.frame(position = events$position,
                node = edges$child[j], time = events$time, edge = j - 1L) else NULL)
        } else for (j in seq_len(nrow(edges))) {
            left <- edges$left[j]
            right <- edges$right[j]
            count <- stats::rpois(1L, rates[j])
            positions <- times <- numeric(count)
            for (k in seq_len(count)) {
                repeat {
                    position <- stats::runif(1L, left, right)
                    if (!position %in% used) break
                }
                used <- c(used, position)
                positions[k] <- position
                times[k] <- stats::runif(1L, low[j], high[j])
            }
            if (length(positions)) {
                mutationRows[[j]] <- data.frame(position = positions,
                    node = edges$child[j], time = times, edge = j - 1L)
            }
        }
        rows <- Filter(Negate(is.null), mutationRows)
        if (length(rows)) {
            mutations <- do.call(rbind, rows)
            mutations <- mutations[order(mutations$position), , drop = FALSE]
            rownames(mutations) <- NULL
            mutations$id <- seq_len(nrow(mutations)) - 1L
            mutations$site <- mutations$id
            mutations$derived_state <- 1L
            sites <- data.frame(id = mutations$site,
                position = mutations$position, ancestral_state = 0L)
        } else {
            mutations <- data.frame(position = numeric(), node = integer(),
                time = numeric(), edge = integer(), id = integer(),
                site = integer(), derived_state = integer())
            sites <- data.frame(id = integer(), position = numeric(),
                                ancestral_state = integer())
        }
        out <- ancestry
        out$sites <- sites
        out$mutations <- mutations
        out$mutationRate <- mutationRate
        out$mutationSeed <- seed
        out$mutationModel <- "continuous infinite-sites binary mutation; root state zero"
        out
    }
    .simupopWithSeed(seed, draw)
}

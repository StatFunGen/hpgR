# From msprime 1.4.1: mutgen_place_mutations, mutation_matrix_choose_root_state
# and mutation_matrix_transition, supplied in Daniel Nachun's tskitR.
# Copyright University of Oxford; GPL-3.0-or-later.

msprimeMutations <- function(ancestry, mutationRate, seed, events = NULL,
                              batchEdges = FALSE, model = "infiniteSites") {
    ancestryCheckTables(ancestry)
    checkInputs(mutationRate = list(positiveScalar, zero = TRUE),
                seed = list(scalarInteger, minimum = 0L), batchEdges = logicalScalar)
    seed <- as.integer(seed)
    model <- match.arg(model, c("infiniteSites", "binary"))
    nodes <- ancestry$nodes
    edges <- ancestry$edges
    draw <- function() {
        low <- nodes$time[match(edges$child, nodes$id)]
        high <- nodes$time[match(edges$parent, nodes$id)]
        left <- if (model == "binary") ceiling(edges$left) else edges$left
        right <- if (model == "binary") ceiling(edges$right) else edges$right
        rates <- mutationRate * (right - left) * (high - low)
        if (any(!is.finite(rates))) stop("mutation exposure is too large", call. = FALSE)
        if (batchEdges && is.null(events) && model == "infiniteSites") {
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
        if (model == "binary") {
            if (!is.null(events)) stop("binary mutation draws do not accept continuous event replay",
                                      call. = FALSE)
            # mutgen_place_mutations assigns mass to integers inside each edge span.
            chosen <- rep(seq_len(nrow(edges)), stats::rpois(nrow(edges), rates))
            position <- floor(stats::runif(length(chosen), left[chosen], right[chosen]))
            time <- low[chosen] + stats::runif(length(chosen)) * (high - low)[chosen]
            mutations <- data.frame(position = position, node = edges$child[chosen],
                time = time, edge = chosen - 1L)
            mutations <- mutations[order(mutations$position, -mutations$time), , drop = FALSE]
            positions <- sort(unique(mutations$position))
            sites <- data.frame(id = seq_along(positions) - 1L,
                position = positions, ancestral_state = rep(0L, length(positions)))
            mutations$id <- seq_len(nrow(mutations)) - 1L
            mutations$site <- match(mutations$position, positions) - 1L
            mutations$derived_state <- integer(nrow(mutations))
            # mutation_matrix_transition for the symmetric binary matrix:
            # each hit flips the state inherited from the older mutation.
            for (j in seq_len(nrow(sites))) {
                rows <- which(mutations$site == sites$id[j])
                parent <- ancestryParent(ancestry, sites$position[j])
                processed <- integer()
                for (row in rows) {
                    path <- ancestryPath(mutations$node[row], parent)
                    old <- processed[mutations$node[processed] %in% path]
                    state <- 0L
                    if (length(old)) {
                        depth <- match(mutations$node[old], path)
                        closest <- tail(old[depth == min(depth)], 1L)
                        state <- mutations$derived_state[closest]
                    }
                    mutations$derived_state[row] <- 1L - state
                    processed <- c(processed, row)
                }
            }
            rownames(mutations) <- NULL
        } else {
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
        }
        out <- ancestry
        out$sites <- sites
        out$mutations <- mutations
        out$mutationRate <- mutationRate
        out$mutationSeed <- seed
        out$mutationModel <- if (model == "binary")
            "discrete-site symmetric binary mutation; root state zero" else
            "continuous infinite-sites binary mutation; root state zero"
        out
    }
    simupopWithSeed(seed, draw)
}

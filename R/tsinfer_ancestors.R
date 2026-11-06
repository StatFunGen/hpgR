# Ancestor construction from tsinfer's ancestor_builder.c, vendored in
# Daniel Nachun's tskitR (revision 7bfa63e). University of Oxford, GPL >= 3.
# Original functions: ancestor_builder_finalise, break_ancestor,
# compute_between_focal_sites, compute_ancestral_states, make_ancestor.

tsiDescriptors <- function(G, times) {
    descriptors <- list()
    for (age in sort(unique(times), decreasing = TRUE)) {
        sites <- which(times == age)
        patterns <- vapply(sites, function(j) paste(G[, j], collapse = ""), "")
        for (key in sort(unique(patterns))) {
            focal <- sites[patterns == key]
            carriers <- which(G[, focal[1L]] == 1L)
            start <- 1L
            if (length(focal) > 1L) for (j in seq_len(length(focal) - 1L)) {
                between <- seq.int(focal[j] + 1L, focal[j + 1L] - 1L)
                between <- between[between > focal[j] & between < focal[j + 1L]]
                split <- any(vapply(between[times[between] > age], function(k) {
                    state <- G[carriers, k]
                    ones <- sum(state == 1L)
                    ones != 0L && ones != sum(state >= 0L)
                }, logical(1)))
                if (split) {
                    descriptors[[length(descriptors) + 1L]] <-
                        list(time = age, focal = focal[start:j])
                    start <- j + 1L
                }
            }
            descriptors[[length(descriptors) + 1L]] <-
                list(time = age, focal = focal[start:length(focal)])
        }
    }
    descriptors
}

tsiMakeAncestor <- function(G, times, focal) {
    ancestor <- rep(-1L, ncol(G))
    carriers <- which(G[, focal[1L]] == 1L)
    age <- times[focal[1L]]
    ancestor[focal] <- 1L
    if (length(focal) > 1L) for (j in seq_len(length(focal) - 1L)) {
        between <- seq.int(focal[j] + 1L, focal[j + 1L] - 1L)
        between <- between[between > focal[j] & between < focal[j + 1L]]
        for (k in between) {
            ancestor[k] <- 0L
            if (times[k] > age) {
                state <- G[carriers, k]
                ones <- sum(state == 1L)
                zeros <- sum(state == 0L)
                ancestor[k] <- if (ones + zeros == 0L) -1L else
                    as.integer(ones >= zeros)
            }
        }
    }
    extend <- function(focus, direction) {
        pool <- which(G[, focus] == 1L)
        minimum <- length(pool) %/% 2L
        disagree <- rep(FALSE, nrow(G))
        last <- focus
        k <- focus + direction
        while (k >= 1L && k <= ncol(G)) {
            ancestor[k] <<- 0L
            last <- k
            state <- G[pool, k]
            ones <- sum(state == 1L)
            zeros <- sum(state == 0L)
            consensus <- as.integer(ones >= zeros)
            differs <- state != consensus & state >= 0L
            remove <- disagree[pool] & differs
            if (times[k] > age) ancestor[k] <<-
                if (ones + zeros == 0L) -1L else consensus
            if (times[k] > age || sum(G[, k] == 1L) > ones)
                disagree[pool[!remove]] <- differs[!remove]
            pool <- pool[!remove]
            if (length(pool) <= minimum) break
            k <- k + direction
        }
        last
    }
    end <- extend(tail(focal, 1L), 1L)
    start <- extend(focal[1L], -1L)
    list(haplotype = ancestor, start = start, end = end,
         focal = focal, time = age)
}

tsiGenerateAncestors <- function(G) {
    times <- colSums(G == 1L) / nrow(G)
    descriptors <- tsiDescriptors(G, times)
    lapply(descriptors, function(d) tsiMakeAncestor(G, times, d$focal))
}

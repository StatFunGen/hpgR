# Copying-HMM recurrence from ancestor_matcher_update_site_likelihood_values
# and traceback in tsinfer/ancestor_matcher.c (Oxford, GPL >= 3), and Dan
# Nachun's tsinfer-matching.R / tsinfer-pipeline.R. Dense state arrays replace
# the C tree-compressed likelihood store. Scores are normalized per site.
# The transition is max(L[u]*(1-rho+rho/n), rho/n), including its tie rule.
# Path compression is an optimization, not an additional inferential model.

tsiMatch <- function(haplotype, reference, nodeIds, recombination, mismatch,
                      numNodes) {
    observed <- which(haplotype >= 0L)
    if (!length(observed)) stop("a haplotype must contain observed alleles", call. = FALSE)
    start <- min(observed)
    end <- max(observed)
    m <- length(haplotype)
    rho <- rep_len(recombination, m)
    mu <- rep_len(mismatch, m)
    values <- rep(1, nrow(reference))
    switch <- matrix(FALSE, nrow(reference), m)
    best <- integer(m)
    scale <- numeric(m)
    for (j in start:end) {
        active <- reference[, j] >= 0L
        stay <- values * (1 - rho[j] + rho[j] / numNodes)
        change <- rho[j] / numNodes
        switch[, j] <- stay <= change
        emit <- ifelse(reference[, j] == haplotype[j] | haplotype[j] < 0L,
                       1 - mu[j], mu[j])
        values <- pmax(stay, change) * emit
        values[!active] <- 0
        best[j] <- which.max(values)
        scale[j] <- max(values)
        if (scale[j] <= 0 || !is.finite(scale[j]))
            stop("copying probabilities cannot explain this haplotype", call. = FALSE)
        values <- values / scale[j]
    }
    state <- integer(m)
    state[end] <- best[end]
    if (end > start) for (j in seq.int(end, start + 1L)) {
        u <- state[j]
        state[j - 1L] <- if (switch[u, j] || reference[u, j - 1L] < 0L)
            best[j - 1L] else u
    }
    matched <- rep(-1L, m)
    matched[start:end] <- reference[cbind(state[start:end], start:end)]
    list(parent = as.integer(nodeIds[state[start:end]]), start = start, end = end,
         matched = matched, logScore = sum(log(scale[start:end])),
         state = state, switch = switch, best = best)
}

#' Infer a tree sequence from phased binary sequences
#'
#' Construct ancestral haplotypes from derived-carrier patterns, match older
#' ancestors with the Li--Stephens copying HMM, then match the sampled copies.
#' Ancestor construction and HMM transitions follow the tsinfer sources in
#' Daniel Nachun's tskitR. Dense arrays make the calculation inspectable in R.
#' @param haplotypes Matrix with phased copies in rows and 0/1 alleles in columns.
#'   Zero denotes the independently established ancestral allele.
#' @param positions Strictly increasing site positions in bases.
#' @param sequenceLength Analyzed length including invariant sites.
#' @param recombinationProbability Scalar or per-site copying-switch parameter.
#'   This is the matcher probability, not a per-base biological rate.
#' @param mismatchProbability Scalar or per-site mismatch probability.
#' @return An hpgR_ancestry object containing inferred node, edge, site and
#'   mutation tables. Initial node times express ancestor ordering. Use
#'   dateAncestry to estimate ages in generations. Ancestor and match records
#'   describe the actual sequence-only computation.
#' @export
inferAncestry <- function(haplotypes, positions, sequenceLength,
                          recombinationProbability = .01,
                          mismatchProbability = 1e-20) {
    G <- as.matrix(haplotypes)
    sequenceLength <- positiveScalar(sequenceLength, "sequenceLength")
    if (!is.numeric(G) || nrow(G) < 2L || ncol(G) < 1L ||
        anyNA(G) || any(!G %in% c(0, 1)) ||
        !is.numeric(positions) || length(positions) != ncol(G) ||
        any(!is.finite(positions)) || any(diff(positions) <= 0) ||
        any(positions < 0 | positions >= sequenceLength))
        stop("supply phased binary copies and ordered positions within the analyzed region",
             call. = FALSE)
    for (p in list(recombinationProbability, mismatchProbability))
        if (!is.numeric(p) || !length(p) %in% c(1L, ncol(G)) ||
            any(!is.finite(p) | p <= 0 | p >= 1))
            stop("matcher probabilities must be in (0,1), scalar or one per site", call. = FALSE)
    if (any(colSums(G) == 0 | colSums(G) == nrow(G)))
        stop("inference positions must be polymorphic; include invariants in sequenceLength",
             call. = FALSE)
    storage.mode(G) <- "integer"
    ancestors <- tsiGenerateAncestors(G)
    oldest <- max(vapply(ancestors, `[[`, 0, "time"))
    nodes <- data.frame(id = 0:1, time = c(oldest + 2, oldest + 1), flags = 0L)
    edges <- data.frame(left = 0, right = sequenceLength, parent = 0L, child = 1L)
    mutations <- data.frame(site = integer(), node = integer(), derived_state = integer())
    reference <- matrix(0L, 2L, ncol(G))
    records <- list()
    addGroup <- function(jobs) {
        ids <- which(nodes$time > jobs[[1L]]$time & nodes$id != 0L)
        fits <- lapply(jobs, function(job) tsiMatch(job$haplotype,
            reference[ids, , drop = FALSE], nodes$id[ids],
            recombinationProbability, mismatchProbability, nrow(nodes)))
        for (i in seq_along(jobs)) {
            job <- jobs[[i]]
            fit <- fits[[i]]
            id <- nrow(nodes)
            nodes <<- rbind(nodes, data.frame(id = id, time = job$time, flags = job$flags))
            edges <<- rbind(edges, tsiEdges(fit, id, positions, sequenceLength))
            changed <- which(job$haplotype >= 0 & job$haplotype != fit$matched)
            if (length(changed)) mutations <<- rbind(mutations,
                data.frame(site = changed - 1L, node = id,
                           derived_state = job$haplotype[changed]))
            reconstructed <- job$haplotype
            within <- fit$start:fit$end
            unknown <- within[reconstructed[within] < 0L]
            reconstructed[unknown] <- fit$matched[unknown]
            reference <<- rbind(reference, reconstructed)
            records[[length(records) + 1L]] <<- list(node = id, time = job$time,
                logScore = fit$logScore, parent = fit$parent,
                start = fit$start, end = fit$end, mutations = changed - 1L)
        }
    }
    times <- vapply(ancestors, `[[`, 0, "time")
    for (age in sort(unique(times), decreasing = TRUE)) {
        jobs <- lapply(ancestors[times == age], function(a) {
            a$flags <- 0L
            a
        })
        addGroup(jobs)
    }
    firstSample <- nrow(nodes)
    addGroup(lapply(seq_len(nrow(G)), function(i)
        list(haplotype = G[i, ], time = 0, flags = 1L)))
    x <- structure(list(nodes = nodes, edges = edges,
        sampleIds = firstSample + seq_len(nrow(G)) - 1L,
        sites = data.frame(id = seq_along(positions) - 1L, position = positions,
                           ancestral_state = 0L),
        mutations = mutations, sequenceLength = sequenceLength,
        timeUnits = "ancestor ordering", inference = list(
            method = "tsinfer ancestor construction and copying HMM",
            ancestors = ancestors, matches = records,
            recombinationProbability = recombinationProbability,
            mismatchProbability = mismatchProbability)),
        class = c("hpgR_ancestry", "list"))
    x <- tsiSimplify(x)
    ancestryCheckTables(x)
    if (!identical(unname(ancestryExtract(x)$observed$haplotypes), unname(G)))
        stop("inferred tables failed to reproduce the observed alleles", call. = FALSE)
    x
}

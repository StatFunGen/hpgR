#' Record and replay finite-copy ancestry
#'
#' Discrete generations sample pooled parental haploid copies uniformly. At linked regions, two pooled donors alternate at Poisson crossovers. Nodes use one-based IDs; edge intervals are left-closed and right-open in bases. Samples are fixed before optional older parent maps complete their ancestry. This is a pooled-copy model, with separate records for observed DNA and simulation truth.
#'
#' @name models_copy_records
#' @param index Child-indexed edge lists from pgEdgeIndex.
#' @param sample_nodes Distinct one-based IDs of sampled nodes.
#' @param left Left endpoint of the local interval in bases, included.
#' @param right Right endpoint of the local interval in bases, excluded.
#' @param edges Edge data frame with parent and child IDs; copy records also use left and right interval coordinates in bases.
#' @param number_nodes Positive integer number of nodes in the edge index.
#' @param record Finite-copy ancestry record returned by pgRecordAncestry.
#' @param mutations Mutation data frame with position in bases, child node, time in generations and derived_state.
#' @param require_complete Logical; require complete sample roots when validating.
#' @param frontier List of local interval records returned by pgLocalTrees or pgAdvanceFrontier.
#' @param sizes Nonempty positive integer diploid sizes in forward generation order, beginning with the ancestral population.
#' @param length_bp Positive analyzed region length in bases.
#' @param n Integer number of sampled haploid copies, at least two.
#' @param map_length Nonnegative region map length in Morgans.
#' @param ancestry_seed Nonnegative integer seed for ancestry draws.
#' @param complete_roots Logical; extend older parent maps until each local sample genealogy has one root.
#' @param max_extension Nonnegative integer limit for generations added before the stated history to complete sample ancestry.
#' @param sample_indices Distinct one-based indices of fixed terminal haploid copies.
#' @param parent_ids Integer node IDs in the parental pool.
#' @param offspring_copies Nonnegative integer number of offspring haploid copies to draw.
#' @param child_ids Integer IDs assigned to offspring copies.
#' @param lengths_bp Positive analyzed lengths in bases, one per independent region.
#' @param mu Nonnegative mutation rate per base per generation.
#' @param map_lengths Nonnegative region map lengths in Morgans, one per independent region.
#' @param mutation_seed Nonnegative integer seed for conditional mutation draws.
#' @param ok Scalar logical assertion.
#' @param message Character error message used when the assertion fails.
#' @param base_seed Nonnegative integer starting seed for region or replicate assignments.
#' @param region Positive integer region identifier.
#' @param mutation Logical; select the mutation rather than ancestry seed assignment.
#' @param x Integer values to validate.
#' @param minimum Minimum accepted integer value.
#' @param name Character argument label used in validation errors.
#' @param copies Vector of ancestral node IDs for the current sampled copies.
#' @param position Position in bases, inside the left-closed, right-open analyzed region.
#' @param observed Observed sample_ids and independent regions, each with analyzed_length_bp, positions_bp and a sample-by-variable-site binary haplotype matrix.
#' @param calibration List containing mu_per_base_per_generation and the declared mutation and sampling models.
#' @param seed Nonnegative integer random seed.
#' @param draw Zero-argument function to evaluate with the supplied seed.
#' @returns pgRecordAncestry returns nodes, edges, sample_nodes, local_trees, size history and root-completion records. pgLocalTrees and pgWalk return local intervals with their used edges and descendant sets. pgMutateAncestry returns binary haplotypes, positions and mutation records; pgReplay independently reconstructs that matrix. pgSimulate returns observed, calibration and truth lists; pgObservedSummaries returns per-region diversity and spectra. Validators return invisible TRUE (pgInteger returns invisible NULL); pgCheck returns invisible NULL on success. Low-level indexing helpers return parent indices, edge lists or local frontier records. pgWithSeed returns the supplied function result while restoring RNG state; pgRegionSeed returns an integer region seed.
#' @examples
#' record <- pgRecordAncestry(rep(4L,3), length_bp=100, n=4, ancestry_seed=42)
#' pgValidateAncestry(record)
#' mutation <- pgMutateAncestry(record, 1e-3, 43)
#' stopifnot(identical(pgReplay(record,mutation$mutations), mutation$haplotypes))
NULL

# g increases forward; t = g_sample - g.
# Sources: Pritchard Sec.2.1, printed76--77; Sec.2.2, printed93--94/Eq2.26;
# canonical Gorroochurn Notes Secs.8.2,8.10--8.11; tskit forward tutorial
# https://tskit-dev.github.io/tutorials/forward_sims.html (read2026-11-03).
# Independent regions have independent copying histories; a shared whole-genome
# diploid pedigree is supplied by simulatePopulation.

#' @rdname models_copy_records
#' @export
pgCheck <- function(ok, message) {
    if (!isTRUE(ok)) stop(message, call. = FALSE)
}

#' @rdname models_copy_records
#' @export
pgInteger <- function(x, minimum = 0, name = deparse(substitute(x))) {
    pgCheck(is.numeric(x) && length(x) > 0L && all(is.finite(x)) &&
                 all(x == floor(x)) && all(x >= minimum),
             paste(name, "must contain finite integers >=", minimum))
}

#' @rdname models_copy_records
#' @export
pgWithSeed <- function(seed, draw) {
    checkInputs(seed = list(scalarInteger, minimum = 0L))
    simupopWithSeed(seed, draw)
}

#' @rdname models_copy_records
#' @export
pgParentDraws <- function(parent_ids, offspring_copies) {
    # Index with sample.int: sample(a_single_large_ID, ...) has another R meaning.
    parent_ids[simupopDrawParents(length(parent_ids), offspring_copies)]
}

#' @rdname models_copy_records
#' @export
pgRegionSeed <- function(base_seed, region, mutation = FALSE) {
    checkInputs(base_seed = list(numericVector, length = 1L, minimum = 0,
                                maximum = .Machine$integer.max, integer = TRUE),
                region = list(numericVector, length = 1L, minimum = 1, integer = TRUE))
    # Avoid adjacent replicate base seeds reusing a previous replicate's region2.
    # These are reproducible distinct seed assignments, not a formal parallel-RNG
    # substream guarantee. Dataset builders also check their complete seed ledger.
    as.integer((104729 * base_seed + 13007 * region + 580007 * mutation) %%
                   2147483646) + 1L
}

#' @rdname models_copy_records
#' @export
pgInheritance <- function(parent_ids, child_ids, length_bp, map_length) {
    first <- pgParentDraws(parent_ids, length(child_ids))
    if (map_length == 0) {
        return(data.frame(left = 0, right = length_bp, parent = first, child = child_ids))
    }
    crossovers <- rpois(length(child_ids), map_length)
    rows <- lapply(seq_along(child_ids), function(i) {
        k <- crossovers[i]
        second <- if (k) pgParentDraws(parent_ids, 1L) else first[i]
        cuts <- c(0, sort(runif(k, 0, length_bp)), length_bp)
        pgCheck(all(diff(cuts) > 0), "crossover positions must be distinct and internal")
        donor <- ifelse(seq_len(k + 1L) %% 2L == 1L, first[i], second)
        data.frame(left = head(cuts, -1L), right = tail(cuts, -1L),
                   parent = donor, child = child_ids[i])
    })
    do.call(rbind, rows)
}

#' @rdname models_copy_records
#' @export
pgEdgeIndex <- function(edges, number_nodes) {
    index <- vector("list", number_nodes)
    parts <- split(edges, edges$child)
    if (length(parts)) index[as.integer(names(parts))] <- unname(parts)
    index
}

#' @rdname models_copy_records
#' @export
pgAtPosition <- function(index, copies, position) {
    unique_copies <- unique(copies)
    rows <- lapply(unique_copies, function(id) {
        part <- index[[id]]
        if (is.null(part) || !nrow(part)) return(NULL)
        selected <- which(part$left <= position & position < part$right)
        pgCheck(length(selected) == 1L, "ancestry has a gap or overlapping parents")
        part[selected, , drop = FALSE]
    })
    if (all(vapply(rows, is.null, logical(1)))) return(NULL)
    pgCheck(!any(vapply(rows, is.null, logical(1))), "mixed recording horizons")
    selected <- do.call(rbind, rows)
    list(copies = selected$parent[match(copies, selected$child)],
         right = min(selected$right), edges = selected)
}

#' @rdname models_copy_records
#' @export
pgWalk <- function(index, sample_nodes, left, right) {
    copies <- sample_nodes
    used <- descendants <- list()
    steps <- 0L
    while (length(unique(copies)) > 1L) {
        step <- pgAtPosition(index, copies, left)
        if (is.null(step)) break
        steps <- steps + 1L
        pgCheck(steps <= length(index), "ancestry cycle")
        right <- min(right, step$right)
        used[[steps]] <- step$edges
        descendants[[steps]] <- lapply(step$edges$child, function(id) which(copies == id))
        copies <- step$copies
    }
    list(left = left, right = right, complete = length(unique(copies)) == 1L,
         copies = copies, edges = if (length(used)) do.call(rbind, used) else NULL,
         descendants = unlist(descendants, recursive = FALSE))
}

#' @rdname models_copy_records
#' @export
pgLocalTrees <- function(record) {
    index <- pgEdgeIndex(record$edges, nrow(record$nodes))
    left <- 0
    trees <- list()
    while (left < record$length_bp) {
        tree <- pgWalk(index, record$sample_nodes, left, record$length_bp)
        pgCheck(tree$right > left, "local ancestry did not advance along the region")
        if (!is.null(tree$edges)) {
            tree$edges$left <- left
            tree$edges$right <- tree$right
        }
        tree$root <- if (tree$complete) unique(tree$copies) else NA_integer_
        tree$root_age <- if (tree$complete) record$nodes$time[tree$root] else NA_real_
        trees[[length(trees) + 1L]] <- tree
        left <- tree$right
    }
    trees
}

#' @rdname models_copy_records
#' @export
pgAdvanceFrontier <- function(index, frontier) {
    result <- list()
    for (block in frontier) {
        if (block$complete) {
            result[[length(result) + 1L]] <- block
            next
        }
        left <- block$left
        while (left < block$right) {
            step <- pgAtPosition(index, block$copies, left)
            pgCheck(!is.null(step), "missing ancestral-completion parent map")
            right <- min(block$right, step$right)
            result[[length(result) + 1L]] <- list(
                left = left, right = right, copies = step$copies,
                complete = length(unique(step$copies)) == 1L)
            left <- right
        }
    }
    result
}

#' @rdname models_copy_records
#' @export
pgRecordAncestry <- function(sizes, length_bp, n, map_length = 0,
                               ancestry_seed = 814901L, complete_roots = TRUE,
                               max_extension = 20000L, sample_indices = seq_len(n)) {
    checkInputs(sizes = list(numericVector, minimum = 1, integer = TRUE),
                n = list(numericVector, length = 1L, minimum = 2, integer = TRUE,
                         maximum = 2 * tail(sizes, 1L)),
                max_extension = list(numericVector, length = 1L, minimum = 0, integer = TRUE),
                length_bp = positiveScalar, map_length = list(positiveScalar, zero = TRUE),
                sample_indices = list(numericVector, length = n, minimum = 1,
                    maximum = 2 * tail(sizes, 1L), integer = TRUE, distinct = TRUE),
                complete_roots = logicalScalar)
    pgWithSeed(ancestry_seed, function() {
        generation <- rep(seq_along(sizes) - 1L, 2L * sizes)
        pools <- split(seq_along(generation), generation)
        edge_blocks <- vector("list", length(sizes) - 1L)
        for (i in seq_len(length(sizes) - 1L)) {
            edge_blocks[[i]] <- pgInheritance(pools[[i]], pools[[i + 1L]],
                                               length_bp, map_length)
        }
        empty_edges <- data.frame(left = numeric(), right = numeric(),
                                  parent = integer(), child = integer())
        edges <- if (length(edge_blocks)) do.call(rbind, edge_blocks) else empty_edges
        nodes <- data.frame(id = seq_along(generation), generation = generation,
                            time = length(sizes) - 1L - generation)
        sample_nodes <- tail(pools, 1L)[[1]][sample_indices]
        record <- list(nodes = nodes, edges = edges, length_bp = length_bp,
                       sample_nodes = sample_nodes)
        frontier <- pgLocalTrees(record)
        index <- pgEdgeIndex(edges, nrow(nodes))
        oldest_pool <- pools[[1]]
        extension <- 0L
        while (complete_roots && !all(vapply(frontier, `[[`, logical(1), "complete"))) {
            pgCheck(extension < max_extension,
                     paste("root completion limit reached; retain failed seed", ancestry_seed))
            extension <- extension + 1L
            older_pool <- length(generation) + seq_len(2L * sizes[1])
            older_edges <- pgInheritance(older_pool, oldest_pool, length_bp, map_length)
            length(index) <- max(older_pool)
            parts <- split(older_edges, older_edges$child)
            index[as.integer(names(parts))] <- unname(parts)
            edge_blocks[[length(edge_blocks) + 1L]] <- older_edges
            generation <- c(generation, rep(-extension, length(older_pool)))
            frontier <- pgAdvanceFrontier(index, frontier)
            oldest_pool <- older_pool
        }
        nodes <- data.frame(id = seq_along(generation), generation = generation,
                            time = length(sizes) - 1L - generation)
        edges <- if (length(edge_blocks)) do.call(rbind, edge_blocks) else empty_edges
        rownames(edges) <- NULL
        record$nodes <- nodes
        record$edges <- edges
        record$local_trees <- pgLocalTrees(record)
        record$sample_indices <- sample_indices
        record$sample_generation <- length(sizes) - 1L
        record$size_history <- data.frame(generation = seq_along(sizes) - 1L, N = sizes)
        record$ancestral_N <- sizes[1]
        record$extension_generations <- extension
        record$root_complete <- all(vapply(record$local_trees, `[[`, logical(1), "complete"))
        record$ancestry_seed <- ancestry_seed
        record$map_length <- map_length
        record$completion_model <- "older full-copy maps; discrete WF at ancestral N"
        record$recording_phase <- ifelse(nodes$generation < 0,
                                         "ancestral_completion", "forward")
        # The sample is fixed before completion. No descendant is added after g_sample.
        record
    })
}

#' @rdname models_copy_records
#' @export
pgValidateAncestry <- function(record, require_complete = TRUE) {
    nodes <- record$nodes
    edges <- record$edges
    pgCheck(identical(nodes$id, seq_len(nrow(nodes))), "node IDs must be consecutive rows")
    pgCheck(all(is.finite(nodes$generation)) &&
                 all(nodes$time == record$sample_generation - nodes$generation),
             "node time/generation inconsistency")
    pgCheck(all(record$sample_nodes %in% nodes$id) && !anyDuplicated(record$sample_nodes),
             "invalid sample node IDs")
    pgCheck(all(nodes$time[record$sample_nodes] == 0), "terminal sample changed")
    pgCheck(all(edges$parent %in% nodes$id) && all(edges$child %in% nodes$id),
             "unknown parent/child ID")
    pgCheck(all(is.finite(edges$left)) && all(is.finite(edges$right)) &&
                 all(edges$left >= 0 & edges$right <= record$length_bp &
                     edges$left < edges$right), "invalid edge intervals")
    pgCheck(all(nodes$generation[edges$parent] == nodes$generation[edges$child] - 1L),
             "parent must be exactly one generation older")
    ordered <- edges[order(edges$child, edges$left), , drop = FALSE]
    first <- !duplicated(ordered$child)
    last <- !duplicated(ordered$child, fromLast = TRUE)
    pgCheck(all(ordered$left[first] == 0) &&
                 all(ordered$right[last] == record$length_bp) &&
                 all(ordered$left[!first] == ordered$right[!last]),
             "child inheritance intervals have a gap or overlap")
    roots <- setdiff(nodes$id, edges$child)
    pgCheck(all(nodes$generation[roots] == min(nodes$generation)),
             "record contains roots at different generations")
    trees <- pgLocalTrees(record)
    complete <- all(vapply(trees, `[[`, logical(1), "complete"))
    pgCheck(identical(record$root_complete, complete), "incorrect root-complete flag")
    if (require_complete) pgCheck(complete, "sample ancestry is an incomplete forest")
    invisible(TRUE)
}

#' @rdname models_copy_records
#' @export
pgMutateAncestry <- function(record, mu, mutation_seed = 814902L) {
    checkInputs(mu = list(positiveScalar, zero = TRUE))
    pgValidateAncestry(record)
    # Retain only sample ancestral edges BELOW each local MRCA. Ancestors above
    # it would add sample-fixed states, not segregating sites or pair differences.
    trees <- pgLocalTrees(record)
    edges <- do.call(rbind, lapply(trees, `[[`, "edges"))
    carriers <- unlist(lapply(trees, `[[`, "descendants"), recursive = FALSE)
    duration <- record$nodes$time[edges$parent] - record$nodes$time[edges$child]
    exposure <- (edges$right - edges$left) * duration
    ancestry <- structure(list(nodes = data.frame(id = record$nodes$id - 1L,
        time = record$nodes$time, flags = as.integer(record$nodes$id %in% record$sample_nodes)),
        edges = data.frame(left = edges$left, right = edges$right,
                           parent = edges$parent - 1L, child = edges$child - 1L),
        sampleIds = record$sample_nodes - 1L, sequenceLength = record$length_bp,
        timeUnits = "generations"), class = c("hpgR_ancestry", "list"))
    mutation <- msprimeMutations(ancestry, mu, mutation_seed, batchEdges = TRUE)$mutations
    G <- matrix(0L, length(record$sample_nodes), nrow(mutation))
    for (j in seq_len(nrow(mutation))) G[carriers[[mutation$edge[j] + 1L]], j] <- 1L
    pgCheck(!ncol(G) || all(colSums(G) > 0 & colSums(G) < nrow(G)),
             "mutation is not on a segregating sample branch")
    list(haplotypes = G, positions_bp = mutation$position,
         mutations = data.frame(position = mutation$position, child = mutation$node + 1L,
                                time = mutation$time, derived_state = mutation$derived_state),
         ancestral_exposure = sum(exposure), mutation_seed = mutation_seed)
}

#' @rdname models_copy_records
#' @export
pgValidateMutations <- function(record, mutations) {
    pgCheck(!anyDuplicated(mutations$position) &&
                 all(diff(mutations$position) > 0), "invalid mutation positions/order")
    pgCheck(all(mutations$child %in% record$nodes$id), "unknown mutation child")
    pgCheck(all(mutations$position >= 0 & mutations$position < record$length_bp),
             "mutation position outside region")
    pgCheck(all(mutations$derived_state == 1L), "invalid infinite-sites derived state")
    index <- pgEdgeIndex(record$edges, nrow(record$nodes))
    for (j in seq_len(nrow(mutations))) {
        mutation <- mutations[j, ]
        part <- pgAtPosition(index, mutation$child, mutation$position)
        pgCheck(!is.null(part), "mutation assigned to a root")
        upper <- record$nodes$time[part$copies]
        lower <- record$nodes$time[mutation$child]
        pgCheck(is.finite(mutation$time) && lower < mutation$time && mutation$time < upper,
                 "mutation time outside its branch")
    }
    invisible(TRUE)
}

#' @rdname models_copy_records
#' @export
pgReplay <- function(record, mutations) {
    # Independent carrier reconstruction: trace the parent map at each site.
    # Do not use cached local trees, descendant lists or the observed matrix.
    pgValidateAncestry(record)
    pgValidateMutations(record, mutations)
    G <- matrix(0L, length(record$sample_nodes), nrow(mutations))
    edges <- record$edges
    for (j in seq_len(nrow(mutations))) {
        x <- mutations$position[j]
        keep <- edges$left <= x & x < edges$right
        parent <- integer(nrow(record$nodes))
        parent[edges$child[keep]] <- edges$parent[keep]
        copies <- record$sample_nodes
        target_generation <- record$nodes$generation[mutations$child[j]]
        while (record$nodes$generation[copies[1]] > target_generation) {
            copies <- parent[copies]
            pgCheck(all(copies > 0), "replay reaches an unrecorded parent")
        }
        G[, j] <- as.integer(copies == mutations$child[j])
    }
    pgCheck(!ncol(G) || all(colSums(G) > 0 & colSums(G) < nrow(G)),
             "replayed mutation is outside segregating sample ancestry")
    G
}

#' @rdname models_copy_records
#' @export
pgSimulate <- function(sizes, lengths_bp = rep(1000, 3), n = 12L, mu = 2e-5,
                         map_lengths = rep(0, length(lengths_bp)),
                         ancestry_seed = 814901L, mutation_seed = 914901L,
                         max_extension = 20000L) {
    pgCheck(length(lengths_bp) > 0L && length(map_lengths) == length(lengths_bp),
             "one map length is needed per region")
    sample_ids <- paste0("S", seq_len(n))
    observed <- list(sample_ids = sample_ids, sampling_time = 0, regions = list())
    truth <- list(regions = list())
    for (i in seq_along(lengths_bp)) {
        record <- pgRecordAncestry(sizes, lengths_bp[i], n, map_lengths[i],
                                    pgRegionSeed(ancestry_seed, i),
                                    max_extension = max_extension)
        pgValidateAncestry(record)
        mutation <- pgMutateAncestry(record, mu,
                                      pgRegionSeed(mutation_seed, i, mutation = TRUE))
        G <- mutation$haplotypes
        variant_ids <- if (ncol(G)) paste0("v", seq_len(ncol(G))) else character()
        dimnames(G) <- list(sample_ids, variant_ids)
        observed$regions[[i]] <- list(region_id = paste0("R", i),
                                     analyzed_length_bp = lengths_bp[i],
                                     positions_bp = mutation$positions_bp, haplotypes = G)
        truth$regions[[i]] <- list(ancestry = record, mutations = mutation$mutations,
                                   ancestral_exposure = mutation$ancestral_exposure,
                                   mutation_seed = mutation$mutation_seed,
                                   sample_node_map = setNames(record$sample_nodes, sample_ids))
    }
    calibration <- list(mu_per_base_per_generation = mu, map_lengths_morgan = map_lengths,
                        mutation_model = "continuous-position infinite sites; root state 0",
                        inheritance_model = "iid uniform pooled parental copies; Poisson switches",
                        regions_independent = TRUE, sample_design = "fixed labelled terminal copies",
                        units = list(positions = "bases", map = "Morgans", time = "generations"))
    list(observed = observed, calibration = calibration, truth = truth)
}

#' @rdname models_copy_records
#' @export
pgObservedSummaries <- function(observed, calibration) {
    # Deliberately accepts no truth argument and performs no ancestor lookup.
    pgCheck(length(observed$sample_ids) >= 2L && !anyDuplicated(observed$sample_ids),
             "invalid observed sample labels")
    pgCheck(is.finite(calibration$mu_per_base_per_generation) &&
                 calibration$mu_per_base_per_generation >= 0, "invalid supplied mutation rate")
    n <- length(observed$sample_ids)
    pgCheck(length(observed$regions) > 0L, "observed dataset has no analyzed regions")
    rows <- spectra <- list()
    for (i in seq_along(observed$regions)) {
        region <- observed$regions[[i]]
        G <- region$haplotypes
        pgCheck(is.matrix(G) && nrow(G) == n && all(G %in% 0:1), "invalid binary matrix")
        pgCheck(identical(rownames(G), observed$sample_ids), "sample rows do not match labels")
        pgCheck(length(region$positions_bp) == ncol(G) &&
                     !anyDuplicated(region$positions_bp) &&
                     all(diff(region$positions_bp) > 0), "invalid observed positions")
        L <- region$analyzed_length_bp
        pgCheck(is.finite(L) && L > 0 &&
                     all(region$positions_bp >= 0 & region$positions_bp < L),
                 "invalid analyzed length or position")
        count <- colSums(G)
        pgCheck(all(count > 0 & count < n), "store variable columns only")
        pi_region <- sum(tskitDiversity(count, n))
        rows[[i]] <- data.frame(region = region$region_id, n = n, length_bp = L,
                                S = ncol(G), pi_region = pi_region, pi_site = pi_region / L)
        spectra[[i]] <- tskitSfs(count, n)
    }
    list(regions = do.call(rbind, rows), sfs = spectra)
}

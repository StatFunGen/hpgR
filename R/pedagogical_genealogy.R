#' Construct and summarize sample genealogies
#'
#' Compressed trees retain sampled leaves and branching ancestors. Node ages and edge lengths are generations before the sample. The discrete genealogy retains simultaneous mergers; the Kingman genealogy uses continuous exponential pair mergers. Tree mutation calls the source-derived mutation kernel.
#'
#' @name models_genealogy
#' @param Ne Positive diploid effective size for a constant-size Kingman model.
#' @param n Integer number of sampled haploid copies, at least two.
#' @param ancestry_seed Nonnegative integer seed for ancestry draws.
#' @param sample_ids Distinct nonempty sample labels in row order.
#' @param tree Compressed node/edge tree returned by the genealogy constructors.
#' @param pair Two distinct one-based sample indices for the pairwise common-ancestor time.
#' @param length_bp Positive analyzed region length in bases.
#' @param mu Nonnegative mutation rate per base per generation.
#' @param mutation_seed Nonnegative integer seed for conditional mutation draws.
#' @param sizes Nonempty positive integer diploid sizes in forward generation order, beginning with the ancestral population.
#' @param max_extension Nonnegative integer limit for generations added before the stated history to complete sample ancestry.
#' @param record Finite-copy ancestry record returned by pgRecordAncestry.
#' @param local One local interval record returned by pgLocalTrees.
#' @param nodes Node data frame with one-based id and nonnegative time before the sample.
#' @param edges Edge data frame with parent and child IDs; copy records also use left and right interval coordinates in bases.
#' @param sample_nodes Distinct one-based IDs of sampled nodes.
#' @param model Character label describing the generating model.
#' @param position Position in bases, inside the left-closed, right-open analyzed region.
#' @returns Tree constructors and compressors return nodes, edges, sample_nodes, sample_ids, root, descendants and model; pgCompressTrees returns a list of local trees. pgTreeSummaries returns pair_time, TMRCA, branch lengths by descendant count and waiting_times. pgMutateTree returns haplotypes, positions_bp, mutations and their spectrum. pgSampleLabels validates labels and returns invisible NULL.
#' @examples
#' tree <- pgKingmanTree(20, n=4, ancestry_seed=42)
#' pgTreeSummaries(tree)$TMRCA
#' dna <- pgMutateTree(tree, length_bp=1000, mu=1e-4, mutation_seed=43)
#' dim(dna$haplotypes)
NULL

# Pritchard Sec.2.2, printed87--94, Eqs.2.14 and2.20--2.26.
# A fixed pair's MRCA age differs from the final two-lineage interval when n > 2.

#' @rdname models_genealogy
#' @export
pgSampleLabels <- function(sample_ids, n) {
    pgCheck(is.character(sample_ids) && length(sample_ids) == n &&
                 !anyNA(sample_ids) && all(nzchar(sample_ids)) &&
                 !anyDuplicated(sample_ids), "invalid sample labels")
}

#' @rdname models_genealogy
#' @export
pgMakeTree <- function(nodes, edges, sample_nodes, sample_ids, model) {
    pgInteger(sample_nodes, 1, "sample_nodes")
    n <- length(sample_nodes)
    pgCheck(n >= 2L && !anyDuplicated(sample_nodes), "invalid tree sample")
    pgSampleLabels(sample_ids, n)
    pgCheck(is.data.frame(nodes) && all(c("id", "time") %in% names(nodes)) &&
                 is.data.frame(edges) && all(c("parent", "child") %in% names(edges)),
             "invalid tree node/edge columns")
    pgInteger(nodes$id, 1, "node IDs")
    pgCheck(!anyDuplicated(nodes$id) && all(sample_nodes %in% nodes$id) &&
                 is.numeric(nodes$time) && all(is.finite(nodes$time)) &&
                 all(nodes$time >= 0), "invalid tree nodes")
    pgCheck(all(edges$parent %in% nodes$id) && all(edges$child %in% nodes$id) &&
                 !anyDuplicated(edges$child) && nrow(edges) == nrow(nodes) - 1L,
             "each nonroot tree node needs one parent")
    roots <- setdiff(nodes$id, edges$child)
    pgCheck(length(roots) == 1L, "tree must have one root")
    age <- setNames(nodes$time, nodes$id)
    edges$length <- unname(age[as.character(edges$parent)] -
                               age[as.character(edges$child)])
    pgCheck(all(is.finite(edges$length)) && all(edges$length > 0) &&
                 all(age[as.character(sample_nodes)] == 0), "invalid tree branch ages")
    children <- split(edges$child, edges$parent)
    pgCheck(!any(sample_nodes %in% edges$parent) &&
                 all(lengths(children) >= 2L) &&
                 setequal(setdiff(nodes$id, edges$parent), sample_nodes),
             "compressed tree leaves must be samples; internal nodes must branch")
    descendants <- setNames(vector("list", nrow(nodes)), nodes$id)
    for (i in seq_along(sample_nodes)) descendants[[as.character(sample_nodes[i])]] <- i
    for (id in nodes$id[order(nodes$time)]) {
        child <- children[[as.character(id)]]
        if (!length(child)) next
        sets <- descendants[as.character(child)]
        joined <- unlist(sets, use.names = FALSE)
        pgCheck(all(lengths(sets) > 0L) && !anyDuplicated(joined),
                 "tree children must have disjoint descendants")
        descendants[[as.character(id)]] <- sort(joined)
    }
    pgCheck(identical(descendants[[as.character(roots)]], seq_len(n)),
             "root must cover every sampled copy")
    nodes$sample_id <- sample_ids[match(nodes$id, sample_nodes)]
    rownames(nodes) <- rownames(edges) <- NULL
    list(nodes = nodes, edges = edges, sample_nodes = sample_nodes,
         sample_ids = sample_ids, n = n, root = roots, descendants = descendants,
         model = model)
}

#' @rdname models_genealogy
#' @export
pgCompressLocal <- function(record, local, sample_ids) {
    # Internal operation on one local tree of an already validated copy record.
    parent <- setNames(local$edges$parent, local$edges$child)
    children <- split(local$edges$child, local$edges$parent)
    retained <- unique(c(record$sample_nodes, local$root,
                         as.integer(names(children)[lengths(children) >= 2L])))
    child <- setdiff(retained, local$root)
    compressed_parent <- vapply(child, function(id) {
        older <- unname(parent[as.character(id)])
        while (!older %in% retained) older <- unname(parent[as.character(older)])
        older
    }, integer(1))
    nodes <- record$nodes[match(retained, record$nodes$id), c("id", "time"), drop = FALSE]
    tree <- pgMakeTree(nodes, data.frame(parent = compressed_parent, child = child),
                         record$sample_nodes, sample_ids,
                         "compressed local tree from the recorded full WF parent maps")
    tree$left <- local$left
    tree$right <- local$right
    tree$source_ancestry_seed <- record$ancestry_seed
    tree
}

#' @rdname models_genealogy
#' @export
pgCompressTree <- function(record, position = 0,
                             sample_ids = paste0("S", seq_along(record$sample_nodes))) {
    pgValidateAncestry(record)
    pgCheck(length(position) == 1L && is.finite(position) && position >= 0 &&
                 position < record$length_bp, "position lies outside the region")
    trees <- pgLocalTrees(record)
    selected <- which(vapply(trees, function(x)
        x$left <= position && position < x$right, logical(1)))
    pgCheck(length(selected) == 1L, "position must select one local tree")
    pgCompressLocal(record, trees[[selected]], sample_ids)
}

#' @rdname models_genealogy
#' @export
pgCompressTrees <- function(record,
                              sample_ids = paste0("S", seq_along(record$sample_nodes))) {
    pgValidateAncestry(record)
    lapply(pgLocalTrees(record), function(local)
        pgCompressLocal(record, local, sample_ids))
}

#' @rdname models_genealogy
#' @export
pgTreeSummaries <- function(tree, pair = 1:2) {
    pgInteger(pair, 1, "pair")
    pgCheck(length(pair) == 2L && !anyDuplicated(pair) && all(pair <= tree$n),
             "pair must name two distinct sample indices")
    nodes <- tree$nodes
    edges <- tree$edges
    descendants <- tree$descendants
    carrier_counts <- lengths(descendants[as.character(edges$child)])
    branch_generations <- numeric(tree$n - 1L)
    for (i in seq_len(tree$n - 1L))
        branch_generations[i] <- sum(edges$length[carrier_counts == i])
    shared <- vapply(descendants, function(x) all(pair %in% x), logical(1))
    pair_age <- min(nodes$time[match(as.integer(names(descendants)[shared]), nodes$id)])
    ages <- sort(unique(nodes$time))
    waiting_times <- setNames(numeric(tree$n - 1L), as.character(2:tree$n))
    parent_age <- nodes$time[match(edges$parent, nodes$id)]
    child_age <- nodes$time[match(edges$child, nodes$id)]
    for (i in seq_len(length(ages) - 1L)) {
        middle <- (ages[i] + ages[i + 1L]) / 2
        k <- sum(child_age < middle & middle < parent_age)
        waiting_times[as.character(k)] <- waiting_times[as.character(k)] +
            ages[i + 1L] - ages[i]
    }
    list(pair_time = pair_age, pair = tree$sample_ids[pair],
         TMRCA = nodes$time[match(tree$root, nodes$id)],
         total_branch_length = sum(edges$length),
         branch_generations = branch_generations, waiting_times = waiting_times,
         last_two_wait = unname(waiting_times["2"]),
         max_children = max(tabulate(match(edges$parent, nodes$id))))
}

#' @rdname models_genealogy
#' @export
pgKingmanTree <- function(Ne, n = 12L, ancestry_seed = 814901L,
                            sample_ids = paste0("S", seq_len(n))) {
    pgCheck(length(Ne) == 1L && is.finite(Ne) && Ne > 0, "Ne must be positive")
    pgInteger(n, 2, "n")
    pgCheck(length(n) == 1L, "n must be scalar")
    pgSampleLabels(sample_ids, n)
    ancestry <- msprimeAncestryCopies(n, Ne, 1L, seed = ancestry_seed)
    tree <- pgMakeTree(
        data.frame(id = ancestry$nodes$id + 1L, time = ancestry$nodes$time),
        data.frame(parent = ancestry$edges$parent + 1L, child = ancestry$edges$child + 1L),
        seq_len(n), sample_ids, "constant-size Kingman exponential pair mergers")
    tree$Ne <- Ne
    tree$ancestry_seed <- ancestry_seed
    tree
}

#' @rdname models_genealogy
#' @export
pgDiscreteGenealogy <- function(sizes, n = 12L, ancestry_seed = 814901L,
                                  max_extension = 20000L,
                                  sample_ids = paste0("S", seq_len(n))) {
    pgInteger(sizes, 1, "sizes")
    pgInteger(n, 2, "n")
    pgInteger(max_extension, 0, "max_extension")
    pgCheck(all(sizes <= .Machine$integer.max / 2) && length(n) == 1L &&
                 n <= 2 * tail(sizes, 1L) && length(max_extension) == 1L,
             "invalid population, sample size or completion limit")
    pgSampleLabels(sample_ids, n)
    pgWithSeed(ancestry_seed, function() {
        nodes <- data.frame(id = seq_len(n), time = rep(0, n))
        edges <- data.frame(parent = integer(), child = integer())
        lineages <- seq_len(n)
        horizon <- length(sizes) - 1L
        parent_sizes <- rev(head(sizes, -1L))
        max_age <- horizon + max_extension
        age <- 0L
        while (length(lineages) > 1L) {
            pgCheck(age < max_age, "discrete genealogy reached its root-completion limit")
            N <- if (age < horizon) parent_sizes[age + 1L] else sizes[1L]
            # Consecutive draws use the same N only up to the next size change.
            remaining <- Inf
            if (age < horizon) {
                tail_sizes <- parent_sizes[seq.int(age + 1L, horizon)]
                change <- which(tail_sizes != N)
                if (length(change)) remaining <- min(remaining, change[1L] - 1L)
                else if (sizes[1L] != N) remaining <- min(remaining, horizon - age)
            }
            k <- length(lineages)
            pool <- 2 * N
            probability <- if (k > pool) 1 else
                -expm1(sum(log1p(-(seq_len(k) - 1L) / pool)))
            batch <- min(remaining, 1024L, max(1L, ceiling(1 / probability)))
            draws <- matrix(simupopDrawParents(pool, k * batch), nrow = k)
            collision <- rep(FALSE, batch)
            for (i in 2:k) for (j in seq_len(i - 1L))
                collision <- collision | draws[i, ] == draws[j, ]
            hit <- which(collision)
            if (!length(hit)) {
                age <- age + batch
                next
            }
            first <- hit[1L]
            age <- age + first
            pgCheck(age <= max_age, "discrete genealogy reached its root-completion limit")
            # Retain every merger in this generation, including separate pair
            # mergers and multifurcations. Batched iid draws change computation,
            # not the discrete parent-copy law; unused later columns are ignored.
            groups <- split(lineages, draws[, first])
            older <- integer(length(groups))
            for (j in seq_along(groups)) {
                group <- groups[[j]]
                if (length(group) == 1L) older[j] <- group
                else {
                    id <- nrow(nodes) + 1L
                    nodes <- rbind(nodes, data.frame(id = id, time = age))
                    edges <- rbind(edges, data.frame(parent = id, child = group))
                    older[j] <- id
                }
            }
            lineages <- older
        }
        tree <- pgMakeTree(nodes, edges, seq_len(n), sample_ids,
                             "exact discrete WF iid ancestral-copy draws")
        tree$ancestry_seed <- ancestry_seed
        tree$size_history <- sizes
        tree$extension_generations <- max(0L, age - horizon)
        tree
    })
}

pgTreeAncestry <- function(tree, length_bp) {
    ids <- tree$nodes$id
    structure(list(nodes = data.frame(id = seq_along(ids) - 1L,
        time = tree$nodes$time, flags = as.integer(ids %in% tree$sample_nodes)),
        edges = data.frame(left = 0, right = length_bp,
            parent = match(tree$edges$parent, ids) - 1L,
            child = match(tree$edges$child, ids) - 1L),
        sampleIds = match(tree$sample_nodes, ids) - 1L,
        sequenceLength = length_bp, timeUnits = "generations"),
        class = c("hpgR_ancestry", "list"))
}

#' @rdname models_genealogy
#' @export
pgMutateTree <- function(tree, length_bp = 1000, mu = 2e-5,
                            mutation_seed = 914901L) {
    pgCheck(length(length_bp) == 1L && is.finite(length_bp) && length_bp > 0,
             "length_bp must be positive")
    pgCheck(length(mu) == 1L && is.finite(mu) && mu >= 0, "mu must be nonnegative")
    exposure <- length_bp * tree$edges$length
    pgCheck(all(is.finite(mu * exposure)), "mutation exposure is too large")
    ancestry <- pgTreeAncestry(tree, length_bp)
    mutated <- msprimeMutations(ancestry, mu, mutation_seed)
    observed <- ancestryData(mutated)$observed
    G <- observed$haplotypes
    variant_ids <- if (ncol(G)) paste0("v", seq_len(ncol(G))) else character()
    dimnames(G) <- list(tree$sample_ids, variant_ids)
    mutations <- mutated$mutations
    list(haplotypes = G, positions_bp = observed$positions,
         mutations = data.frame(position = mutations$position,
             child = tree$nodes$id[mutations$node + 1L],
             time = mutations$time, derived_state = mutations$derived_state),
         sfs = tskitSfs(colSums(G), tree$n),
         ancestral_exposure = sum(exposure), mutation_seed = mutation_seed)
}

#' @rdname models_genealogy_controls
#' @export
recordedSubset <- function(tree, selected = 1:4) {
    sets <- lapply(tree$descendants, intersect, y = selected)
    valid <- lengths(sets) >= 2
    keys <- vapply(sets[valid], ptiKey, character(1))
    ages <- tree$nodes$time[match(as.integer(names(sets)[valid]), tree$nodes$id)]
    ages <- tapply(ages, keys, min)
    ages <- setNames(as.numeric(ages), names(ages))
    root_key <- ptiKey(selected)
    root_age <- ages[root_key]
    ages <- ages[ages <= root_age]
    clades <- lapply(names(ages), function(k) as.integer(strsplit(k, ",", fixed = TRUE)[[1]]))
    clades <- c(lapply(selected, function(i) i), clades)
    list(tree = cladeTree(clades, selected), ages = ages)
}

# Derived-allele counts for the small polarized-sequence examples.
#' @rdname models_genealogy_controls
#' @export
sfs <- function(G) {
    tskitSfs(colSums(G), nrow(G))
}

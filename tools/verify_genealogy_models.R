#!/usr/bin/env Rscript
# Optional --truth-file PRIVATE_P1_TRUTH_RDS checks the saved P1 sample trees.
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
lecture <- Sys.getenv("HPGR_COURSE_DIR")
library(hpgR)
library(hpgR)
library(hpgR)

failure <- function(expr) stopifnot(inherits(try(expr, silent = TRUE), "try-error"))
close_mean <- function(x, expected, label, sigma = 6) {
    difference <- abs(mean(x) - expected)
    se <- sd(x) / sqrt(length(x))
    if (!(difference <= sigma * se + 1e-9))
        stop(sprintf("%s: mean %.7g versus expected %.7g, MCSE %.7g",
                     label, mean(x), expected, se), call. = FALSE)
    invisible(c(mean = mean(x), expected = expected, MCSE = se))
}
compare_record <- function(record) {
    samples <- paste0("S", seq_along(record$sample_nodes))
    local_trees <- pg_local_trees(record)
    compressed <- pg_compress_trees(record, samples)
    for (index in seq_along(local_trees)) {
        local <- local_trees[[index]]
        tree <- compressed[[index]]
        summary <- pg_tree_summaries(tree)
        full_exposure <- tabulate(lengths(local$descendants), nbins = tree$n - 1L)
        stopifnot(identical(summary$branch_generations, as.numeric(full_exposure)),
                  summary$TMRCA == local$root_age,
                  summary$total_branch_length == nrow(local$edges),
                  identical(tree$sample_nodes, record$sample_nodes),
                  identical(tree$sample_ids, samples), tree$root == local$root)
        # Every compressed edge follows the original parent path exactly, with
        # matching original IDs and one generation counted for each full edge.
        parent <- setNames(local$edges$parent, local$edges$child)
        for (i in seq_len(nrow(tree$edges))) {
            edge <- tree$edges[i, ]
            current <- edge$child
            steps <- 0L
            while (current != edge$parent) {
                current <- unname(parent[as.character(current)])
                steps <- steps + 1L
            }
            stopifnot(steps == edge$length)
        }
        # Independent pair trace from the uncompressed maps.
        copies <- record$sample_nodes[1:2]
        pair_time <- 0L
        while (copies[1L] != copies[2L]) {
            copies <- unname(parent[as.character(copies)])
            pair_time <- pair_time + 1L
        }
        stopifnot(summary$pair_time == pair_time,
                  sum(summary$waiting_times) == summary$TMRCA,
                  sum(as.numeric(names(summary$waiting_times)) *
                          summary$waiting_times) == summary$total_branch_length)
    }
    invisible(TRUE)
}

designs <- list(list(sizes = rep(4L, 5L), n = 6L, map = 0),
                list(sizes = c(5L, 2L, 4L), n = 8L, map = 0),
                list(sizes = 3L, n = 3L, map = 0),
                list(sizes = rep(3L, 5L), n = 4L, map = .3))
for (i in seq_along(designs)) {
    d <- designs[[i]]
    record <- pg_record_ancestry(d$sizes, 1000, d$n, d$map, ancestry_seed = 61800L + i)
    compare_record(record)
}
cat("PASS: compression preserves local intervals, original node/sample IDs, ages and exposures\n")

truth_index <- match("--truth-file", args)
if (!is.na(truth_index)) {
    stopifnot(truth_index < length(args), file.exists(args[truth_index + 1L]))
    truth <- readRDS(args[truth_index + 1L])
    compared <- 0L
    for (dataset in truth) for (region in dataset$regions) {
        compare_record(region$ancestry)
        compared <- compared + 1L
    }
    cat(sprintf("PASS: exact sample/tree identities against %d stored P1 regions\n", compared))
}

# A fixed two-copy sample has a geometric MRCA age on {1,2,...}.
pair_N <- 4L
pair_trees <- lapply(seq_len(5000L), function(i)
    pg_discrete_genealogy(pair_N, n = 2L, ancestry_seed = 610000L + i))
pair_times <- vapply(pair_trees, function(x) pg_tree_summaries(x)$pair_time, numeric(1))
close_mean(pair_times, 2 * pair_N, "geometric pair mean")
stopifnot(abs(var(pair_times) - 2 * pair_N * (2 * pair_N - 1)) <
              .15 * 2 * pair_N * (2 * pair_N - 1))
for (t in c(1, 4, 8, 16)) {
    probability <- (1 - 1 / (2 * pair_N))^t
    stopifnot(abs(mean(pair_times > t) - probability) <
                  6 * sqrt(probability * (1 - probability) / length(pair_times)))
}
cat("PASS: discrete pair mean, variance and survival match the geometric law\n")

# Independent Kingman simulations test the exact implemented continuous model.
Ne <- 80
n <- 8L
kingman <- lapply(seq_len(3000L), function(i)
    pg_tree_summaries(pg_kingman_tree(Ne, n, ancestry_seed = 620000L + i)))
kingman_stat <- function(name) vapply(kingman, `[[`, numeric(1), name)
close_mean(kingman_stat("pair_time"), 2 * Ne, "Kingman fixed pair mean")
close_mean(kingman_stat("last_two_wait"), 2 * Ne, "Kingman final two-lineage mean")
close_mean(kingman_stat("TMRCA"), 4 * Ne * (1 - 1 / n), "Kingman MRCA mean")
close_mean(kingman_stat("total_branch_length"), 4 * Ne * sum(1 / seq_len(n - 1L)),
           "Kingman total branch mean")
branches <- vapply(kingman, `[[`, numeric(n - 1L), "branch_generations")
waits <- vapply(kingman, `[[`, numeric(n - 1L), "waiting_times")
for (i in seq_len(n - 1L)) {
    close_mean(branches[i, ], 4 * Ne / i, paste("branch bin", i))
    k <- i + 1L
    close_mean(waits[i, ], 4 * Ne / (k * (k - 1)), paste("waiting stage", k))
}
stopifnot(all(vapply(kingman, `[[`, numeric(1), "max_children") == 2),
          any(kingman_stat("pair_time") != kingman_stat("last_two_wait")))
cat("PASS: Kingman pair, MRCA, branch-SFS and all stage waiting-time expectations\n")

# At finite N, an iid parental draw can merge three copies, or separate pairs,
# in the same generation. A Kingman realization never has these events.
small <- lapply(seq_len(800L), function(i)
    pg_discrete_genealogy(4L, n = 6L, ancestry_seed = 630000L + i))
stopifnot(any(vapply(small, function(x)
    pg_tree_summaries(x)$max_children >= 3L, logical(1))))
simultaneous <- vapply(small, function(x) {
    internal <- x$nodes$time[!x$nodes$id %in% x$sample_nodes]
    anyDuplicated(internal) > 0L
}, logical(1))
stopifnot(any(simultaneous))
at_one <- vapply(small, function(x) {
    if (x$nodes$time[match(x$root, x$nodes$id)] <= 1) return(1)
    parent_age <- x$nodes$time[match(x$edges$parent, x$nodes$id)]
    child_age <- x$nodes$time[match(x$edges$child, x$nodes$id)]
    sum(child_age <= 1 & parent_age > 1)
}, numeric(1))
close_mean(at_one, 8 * (1 - (1 - 1 / 8)^6), "occupied parental copies")
cat("PASS: exact finite-N occupancies preserve multifurcations and simultaneous mergers\n")

# Reversing a changing history must use N of the parent generation. In this
# design, the terminal n=4 copies choose from just two parental copies at age1.
changing <- lapply(seq_len(200L), function(i)
    pg_discrete_genealogy(c(20L, 1L, 20L), n = 4L, ancestry_seed = 640000L + i))
stopifnot(all(vapply(changing, function(x) {
    parent_age <- x$nodes$time[match(x$edges$parent, x$nodes$id)]
    child_age <- x$nodes$time[match(x$edges$child, x$nodes$id)]
    sum(child_age <= 1 & parent_age > 1) <= 2L
}, logical(1))))
cat("PASS: changing histories use parental N with the original generation convention\n")

# Compare two independent computations of the same discrete WF law. This is a
# distributional test; seeds need not reproduce each other's realized parent map.
full <- sampled <- numeric(160L)
for (i in seq_along(full)) {
    record <- pg_record_ancestry(rep(8L, 4L), 1000, 6L, ancestry_seed = 650000L + i)
    full[i] <- pg_tree_summaries(pg_compress_tree(record))$total_branch_length
    sampled[i] <- pg_tree_summaries(pg_discrete_genealogy(rep(8L, 4L), n = 6L,
                                  ancestry_seed = 660000L + i))$total_branch_length
}
stopifnot(abs(mean(full) - mean(sampled)) <
              6 * sqrt(var(full) / length(full) + var(sampled) / length(sampled)))
cat("PASS: full-copy maps and sampled-lineage draws agree in distribution\n")

tree <- pg_kingman_tree(80L, n = 6L, ancestry_seed = 670001L)
summary <- pg_tree_summaries(tree)
expected <- 2e-5 * 1000 * summary$branch_generations
mutations <- lapply(seq_len(3000L), function(i)
    pg_mutate_tree(tree, mu = 2e-5, mutation_seed = 671000L + i))
spectra <- vapply(mutations, `[[`, numeric(5L), "sfs")
for (i in seq_len(5L)) {
    close_mean(spectra[i, ], expected[i], paste("mutation bin", i))
    stopifnot(abs(var(spectra[i, ]) - expected[i]) < .2 * max(expected[i], 1))
}
for (mutation in mutations[1:10]) {
    G <- mutation$haplotypes
    edge <- match(mutation$mutations$child, tree$edges$child)
    upper <- tree$nodes$time[match(tree$edges$parent[edge], tree$nodes$id)]
    lower <- tree$nodes$time[match(mutation$mutations$child, tree$nodes$id)]
    stopifnot(identical(rownames(G), tree$sample_ids),
              all(colSums(G) > 0 & colSums(G) < tree$n),
              all(mutation$mutations$time > lower & mutation$mutations$time < upper),
              all(diff(mutation$positions_bp) > 0),
              mutation$ancestral_exposure == 1000 * summary$total_branch_length)
}
pair_differences <- vapply(mutations, function(x)
    sum(x$haplotypes[1L, ] != x$haplotypes[2L, ]), numeric(1))
close_mean(pair_differences, 2 * 2e-5 * 1000 * summary$pair_time,
           "mutations below the fixed pair MRCA")
zero <- pg_mutate_tree(tree, mu = 0)
stopifnot(ncol(zero$haplotypes) == 0L, all(zero$sfs == 0),
          identical(pg_mutate_tree(tree), pg_mutate_tree(tree)))
cat("PASS: conditional infinite-sites mutation exposure, SFS moments and carrier identities\n")

failure(pg_kingman_tree(0L))
failure(pg_kingman_tree(80L, n = 1L))
failure(pg_kingman_tree(80L, n = c(2L, 3L)))
failure(pg_kingman_tree(80L, n = 3L, sample_ids = rep("a", 3L)))
failure(pg_discrete_genealogy(0L))
failure(pg_discrete_genealogy(3L, n = 7L))
failure(pg_discrete_genealogy(3L, n = 2L, max_extension = 0L))
failure(pg_discrete_genealogy(c(3L, 0L), n = 2L))
failure(pg_mutate_tree(tree, mu = -1))
failure(pg_mutate_tree(tree, length_bp = 0))
failure(pg_tree_summaries(tree, pair = c(1L, 1L)))
failure(pg_tree_summaries(tree, pair = c(1L, 7L)))
record <- pg_record_ancestry(3L, 1000, 3L)
failure(pg_compress_tree(record, position = -1))
failure(pg_compress_tree(record, position = 1000))
set.seed(420L)
before <- .Random.seed
invisible(pg_kingman_tree(80L))
invisible(pg_discrete_genealogy(80L))
invisible(pg_mutate_tree(tree))
stopifnot(identical(before, .Random.seed),
          identical(pg_kingman_tree(80L), pg_kingman_tree(80L)),
          identical(pg_discrete_genealogy(80L), pg_discrete_genealogy(80L)))
rm(".Random.seed", envir = .GlobalEnv)
invisible(pg_kingman_tree(80L))
invisible(pg_discrete_genealogy(80L))
invisible(pg_mutate_tree(tree))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
cat("PASS: seed replay, invalid-input failures and caller RNG restoration\n")

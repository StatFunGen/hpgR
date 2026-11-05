#!/usr/bin/env Rscript
# Independent checks of sequence-only topology constraints and mutation-clock fits.
args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg) != 1L) stop("Run this file with Rscript")
tool_path <- normalizePath(sub("^--file=", "", file_arg))
vignettes <- Sys.getenv("HPGR_VIGNETTE_DIR")
if (!nzchar(vignettes)) stop("Set HPGR_VIGNETTE_DIR to the canonical vignette repository")
library(hpgR)

checks <- character()
check <- function(name, test) {
    if (!isTRUE(test)) stop("FAILED: ", name)
    checks <<- c(checks, name)
    cat("PASS:", name, "\n")
}
near <- function(a, b, tolerance = 1e-6)
    isTRUE(all.equal(unname(a), unname(b), tolerance = tolerance, check.attributes = FALSE))
errors <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
G_from_counts <- function(keys, counts) {
    columns <- list()
    for (i in seq_along(keys)) if (counts[i] > 0L) {
        column <- as.numeric(seq_len(4) %in% as.integer(strsplit(keys[i], ",", fixed = TRUE)[[1]]))
        columns <- c(columns, rep(list(column), counts[i]))
    }
    if (!length(columns)) return(matrix(numeric(), 4, 0))
    do.call(cbind, columns)
}
canonical <- function(tree) paste(sort(vapply(pti_tree_clades(tree), pti_key, character(1))), collapse = ";")
trees <- pti_rooted_trees()
check("15 distinct rooted four-tip binaries", length(trees) == 15L &&
      length(unique(vapply(trees, canonical, character(1)))) == 15L)
ranks <- lapply(trees, pti_rank_designs)
check("18 ranked histories, including all balanced orders", sum(lengths(ranks)) == 18L &&
      sum(lengths(ranks) == 2L) == 3L && sum(lengths(ranks) == 1L) == 12L)
check("four/three/two lineage exposure columns", all(vapply(unlist(ranks, recursive = FALSE),
    function(d) identical(as.numeric(colSums(d$A)), c(4, 3, 2)), logical(1))))

balanced <- list(list(1L, 2L), list(3L, 4L))
cat_tree <- list(list(list(1L, 2L), 3L), 4L)
balanced_design <- Filter(function(d) identical(d$order, c("1,2", "3,4", "1,2,3,4")),
                          pti_rank_designs(balanced))[[1]]
cat_design <- pti_rank_designs(cat_tree)[[1]]
check("balanced branch exposure agrees with direct node subtraction",
      near(setNames(as.vector(balanced_design$A %*% c(20, 40, 40)), balanced_design$edge_keys)
           [c("1", "2", "3", "4", "1,2", "3,4")], c(20, 20, 60, 60, 80, 40)))
check("caterpillar branch exposure agrees with direct node subtraction",
      near(setNames(as.vector(cat_design$A %*% c(20, 30, 50)), cat_design$edge_keys)
           [c("1", "2", "3", "4", "1,2", "1,2,3")], c(20, 20, 50, 100, 30, 50)))

empty <- matrix(numeric(), 4, 0)
singletons <- diag(4)
nested <- G_from_counts(c("1,2", "1,2,3"), c(1, 1))
disjoint <- G_from_counts(c("1,2", "3,4"), c(1, 1))
conflict <- G_from_counts(c("1,2", "1,3"), c(1, 1))
check("empty and singleton-only data preserve all candidates",
      length(pti_compatible_trees(empty)) == 15L && length(pti_compatible_trees(singletons)) == 15L)
check("nested, disjoint, and incompatible carrier constraints",
      length(pti_compatible_trees(nested)) == 1L && length(pti_compatible_trees(disjoint)) == 1L &&
      length(pti_compatible_trees(conflict)) == 0L)
check("monomorphic columns do not add below-root branches",
      length(pti_compatible_trees(cbind(nested, rep(0, 4), rep(1, 4)))) == 1L &&
      ncol(pti_polymorphic(cbind(empty, rep(0, 4), rep(1, 4)))) == 0L)
check("invalid inputs are rejected", errors(pti_compatible_trees(matrix(2, 4, 2))) &&
      errors(pti_compatible_trees(matrix(0, 3, 1))) &&
      errors(pti_linked_constraints(nested, c(1, 1))) && errors(pti_fit_ages(nested, 0, 1000)))

zero_fit <- pti_fit_ages(empty, 1e-5, 1000)
check("zero counts retain every ranking and exact zero ages", length(zero_fit$best_fits) == 18L &&
      all(vapply(zero_fit$best_fits, function(x) all(x$durations == 0) && x$loglik == 0, logical(1))))
# Integer branch counts proportional to known branch lengths have an exact MLE.
for (design in list(balanced_design, cat_design)) {
    durations <- if (identical(design, balanced_design)) c(20, 40, 40) else c(20, 30, 50)
    counts <- as.integer(as.vector(design$A %*% durations))
    G <- G_from_counts(design$edge_keys, counts)
    fit <- pti_fit_ages(G, 1, 1)
    got <- fit$best_fits[[1]]
    expected <- setNames(cumsum(durations), design$order)
    check(paste("exact non-boundary MLE:", paste(design$order, collapse = "/")),
          near(got$ages[names(expected)], expected, 1e-5) && got$kkt_violation < 1e-5)
}

observed <- readRDS(file.path(vignettes, "data", "popgen_forward_p1", "observed.rds"))
calibration <- readRDS(file.path(vignettes, "data", "popgen_forward_p1", "calibration.rds"))
region <- observed$D3$regions[[1]]
G <- region$haplotypes[1:4, , drop = FALSE]
mu <- calibration$D3$mu_per_base_per_generation
L <- region$analyzed_length_bp
fit <- pti_fit_ages(G, mu, L)
check("fixed D3 sample stays unresolved", ncol(pti_polymorphic(G)) == 7L &&
      fit$consensus$compatible_count == 3L &&
      identical(vapply(fit$consensus$supported_clades, pti_key, character(1)), "1,3,4"))
check("D3 exact zero edge and all tied refinements retained", length(fit$best_fits) == 3L &&
      all(vapply(fit$best_fits, function(x) near(x$durations, c(50, 0, 75), 1e-5) &&
          x$boundary && x$kkt_violation < 1e-6, logical(1))))
check("consensus age estimates are 50 and 125 generations",
      near(fit$node_ages$estimate_min, c(50, 125), 1e-5) &&
      near(fit$node_ages$estimate_max, c(50, 125), 1e-5))
rescaled <- pti_fit_ages(G, 2 * mu, L)
length_rescaled <- pti_fit_ages(G, mu, 2 * L)
check("mutation rate and analyzed-length rescaling", near(rescaled$node_ages$estimate_min,
      fit$node_ages$estimate_min / 2) && near(length_rescaled$node_ages$estimate_min,
      fit$node_ages$estimate_min / 2))
permutation <- c(4L, 2L, 1L, 3L)
permuted <- pti_fit_ages(G[permutation, rev(seq_len(ncol(G))), drop = FALSE], mu, L)
old_supported <- sort(vapply(permuted$consensus$supported_clades,
    function(x) pti_key(permutation[x]), character(1)))
check("sample and site permutations preserve inference", near(permuted$best_loglik, fit$best_loglik) &&
      identical(old_supported, "1,3,4") &&
      near(permuted$node_ages$estimate_min, fit$node_ages$estimate_min))

pairs <- pti_pair_intervals(G, mu, L)
q <- .05 / 6
independent_lower <- ifelse(pairs$differences == 0, 0,
    qchisq(q / 2, 2 * pairs$differences) / (4 * mu * L))
independent_upper <- qchisq(1 - q / 2, 2 * (pairs$differences + 1)) / (4 * mu * L)
root_ci <- pti_clade_interval(pairs)
check("Bonferroni Poisson intervals agree with chi-squared calculation",
      near(pairs$lower, independent_lower) && near(pairs$upper, independent_upper) && nrow(pairs) == 6L)
check("MRCA bounds are maxima of simultaneous pair bounds",
      near(root_ci$lower, max(pairs$lower)) && near(root_ci$upper, max(pairs$upper)) &&
      root_ci$simultaneous_coverage_at_least == .95)
zero_pairs <- pti_pair_intervals(empty, mu, L)
check("zero differences retain finite upper uncertainty", all(zero_pairs$lower == 0) &&
      all(is.finite(zero_pairs$upper)) && all(zero_pairs$upper > 0))

linked_region <- observed$D4$regions[[1]]
linked_G <- linked_region$haplotypes[1:4, , drop = FALSE]
linked <- pti_linked_constraints(linked_G, linked_region$positions_bp)
quarter_counts <- vapply(seq_len(4), function(j) {
    in_window <- linked_region$positions_bp >= (j - 1) * L / 4 &
                 linked_region$positions_bp < j * L / 4
    length(pti_compatible_trees(linked_G[, in_window, drop = FALSE]))
}, integer(1))
check("fixed-quarter public linked constraints", identical(quarter_counts, c(1L, 1L, 3L, 3L)) &&
      length(pti_compatible_trees(linked_G)) == 0L && nrow(linked$witnesses) == 8L)
check("linked minimum blocks and all alternative cuts retained", linked$minimum_blocks == 3L &&
      linked$minimum_changes == 2L && length(linked$minimal_cut_sets) == 12L && nrow(linked$cut_gaps) == 7L)
check("two disjoint mutation-pattern witnesses", nrow(linked$disjoint_witnesses) == 2L &&
      near(linked$disjoint_witnesses$left, c(152.726451, 643.411416), 1e-5) &&
      near(linked$disjoint_witnesses$right, c(300.095358, 835.802690), 1e-5))

# Brute-force every cut set independently for small examples; no truth or hidden
# breakpoint is used by either the dynamic program or this exhaustive check.
brute <- function(G) {
    G <- pti_polymorphic(G); m <- ncol(G)
    if (m <= 1L) return(list(integer()))
    possible <- list()
    for (mask in 0:(2^(m - 1L) - 1L)) {
        cuts <- which(as.logical(intToBits(mask)[seq_len(m - 1L)]))
        bounds <- c(0L, cuts, m)
        okay <- all(vapply(seq_len(length(bounds) - 1L), function(i)
            length(pti_compatible_trees(G[, seq.int(bounds[i] + 1L, bounds[i + 1L]), drop = FALSE])) > 0,
            logical(1)))
        if (okay) possible[[length(possible) + 1L]] <- cuts
    }
    possible[lengths(possible) == min(lengths(possible))]
}
set_key <- function(x) paste(sort(vapply(x, function(a) paste(a, collapse = ","), character(1))), collapse = ";")
dp_cases <- list(empty, nested, conflict, G_from_counts(c("1,2", "1,3", "2,4", "1,2,3", "1,4"), rep(1, 5)),
                 linked_G[, 2:9, drop = FALSE])
check("dynamic programming equals independent exhaustive partitions", all(vapply(dp_cases, function(g) {
    result <- pti_linked_constraints(g, seq_len(ncol(g)))
    identical(set_key(result$minimal_cut_sets), set_key(brute(g)))
}, logical(1))))

# Predeclared independent seeds, neither filtered nor selected using a result.
# Generate branch mutations outside the estimator, then provide only G, mu, L.
saved_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
seeds <- 624001:624200
covered <- logical(length(seeds)); all_pair_covered <- covered
rank_kkt <- numeric(20); root_estimates <- numeric(20)
for (s in seq_along(seeds)) {
    set.seed(seeds[s])
    design <- if (s %% 2L) cat_design else balanced_design
    durations <- if (s %% 2L) c(30, 40, 70) else c(30, 40, 70)
    known_ages <- setNames(cumsum(durations), design$order)
    branch_lengths <- as.vector(design$A %*% durations)
    rate <- .2
    counts <- rpois(length(branch_lengths), rate * branch_lengths)
    sim_G <- G_from_counts(design$edge_keys, counts)
    intervals <- pti_pair_intervals(sim_G, rate, 1)
    true_pair_ages <- vapply(seq_len(nrow(intervals)), function(i) {
        pair <- c(intervals$i[i], intervals$j[i])
        clades <- pti_tree_clades(design$tree, include_tips = FALSE)
        ancestor <- Filter(function(x) all(pair %in% x), clades)
        unname(known_ages[pti_key(ancestor[[which.min(lengths(ancestor))]])])
    }, numeric(1))
    ci <- pti_clade_interval(intervals)
    covered[s] <- ci$lower <= 140 && ci$upper >= 140
    all_pair_covered[s] <- all(intervals$lower <= true_pair_ages & intervals$upper >= true_pair_ages)
    if (s <= 20L) {
        inferred <- pti_fit_ages(sim_G, rate, 1)
        root_estimates[s] <- inferred$node_ages$estimate_min[inferred$node_ages$clade == "1,2,3,4"]
        rank_kkt[s] <- max(vapply(inferred$best_fits, `[[`, numeric(1), "kkt_violation"))
        stopifnot(inferred$consensus$compatible_count >= 1L)
    }
}
if (is.null(saved_seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
} else assign(".Random.seed", saved_seed, envir = .GlobalEnv)
check("unfiltered simulated sequences retain nominal simultaneous uncertainty",
      mean(covered) >= .9 && mean(all_pair_covered) >= .9)
check("convex optimum certificate across independent simulated inputs", max(rank_kkt) < 1e-4)
cat(sprintf("SIMULATION: %d fixed seeds; root coverage %.3f; all-pair coverage %.3f; max KKT %.3g\n",
            length(seeds), mean(covered), mean(all_pair_covered), max(rank_kkt)))
cat("SIMULATION: first 20 root point estimates:", paste(round(root_estimates, 3), collapse = ", "), "\n")
cat(sprintf("All %d sequence-only tree-inference checks passed.\n", length(checks)))

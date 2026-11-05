test_that("founder IBD and diversity have the same expected copying decay", {
    copies <- 4L
    initial <- cbind(founder = seq_len(copies), allele = rep(0:1, each = 2L))
    parents <- as.matrix(expand.grid(rep(list(seq_len(copies)), copies)))
    pair <- utils::combn(copies, 2L)
    initialH <- mean(initial[pair[1L, ], "allele"] != initial[pair[2L, ], "allele"])
    values <- t(apply(parents, 1L, function(index) {
        child <- copyGenotypes(initial, index)
        c(F = mean(child[pair[1L, ], 1L] == child[pair[2L, ], 1L]),
          H = mean(child[pair[1L, ], 2L] != child[pair[2L, ], 2L]))
    }))
    expect_equal(nrow(values), copies^copies)
    expect_equal(mean(values[, "F"]), 1 / copies)
    expect_equal(mean(1 - values[, "F"]), 1 - 1 / copies)
    expect_equal(mean(values[, "H"]) / initialH, 1 - 1 / copies)
    expect_true(any(abs(1 - values[, "F"] - values[, "H"] / initialH) > .01))
})

test_that("a supplied terminal extension adds singleton branch exposure", {
    # Pritchard's rapid-growth approximation adds time before any sample merger.
    tree <- pg_kingman_tree(80, n = 12L, ancestry_seed = 814903L)
    nodes <- tree$nodes[c("id", "time")]
    extension <- 40
    internal <- !nodes$id %in% tree$sample_nodes
    nodes$time[internal] <- nodes$time[internal] + extension
    extended <- pg_make_tree(nodes, tree$edges[c("parent", "child")],
        tree$sample_nodes, tree$sample_ids, "specified terminal extension")
    base <- pg_tree_summaries(tree)
    altered <- pg_tree_summaries(extended)
    added <- c(tree$n * extension, rep(0, tree$n - 2L))
    expect_equal(altered$branch_generations - base$branch_generations, added,
                 tolerance = 1e-11)
    expect_equal(altered$TMRCA - base$TMRCA, extension)
    expect_equal(altered$pair_time - base$pair_time, extension)
    expect_equal(altered$total_branch_length - base$total_branch_length,
                 tree$n * extension)
})

test_that("diversity and segregating sites recover the same neutral expectation", {
    n <- 12L
    Ne <- 80
    mu <- 2e-5
    L <- 1000
    spectrum <- expectedSFS(n, Ne, mu, L)
    weights <- 2 * seq_len(n - 1L) * (n - seq_len(n - 1L)) / (n * (n - 1L))
    expect_equal(sum(spectrum * weights) / (4 * mu * L), Ne)
    expect_equal(sum(spectrum) / (sum(1 / seq_len(n - 1L)) * 4 * mu * L), Ne)

    # Integer counts give an independent check of the old specified-growth fit.
    B <- 48
    counts <- round(B * (spectrum + c(n * 40 * mu * L, rep(0, n - 2L))))
    ancestral <- sum(counts[-1L]) / (B * 4 * mu * L * sum(1 / seq.int(2, n - 1L)))
    generations <- (counts[1L] - B * 4 * mu * L * ancestral) / (B * mu * L * n)
    objective <- function(par) {
        expected <- B * mu * L * (4 * par[1L] / seq_len(n - 1L) +
            c(n * par[2L], rep(0, n - 2L)))
        -sum(stats::dpois(counts, expected, log = TRUE))
    }
    fit <- stats::optim(c(20, 100), objective, method = "L-BFGS-B",
        lower = c(.01, 0), upper = c(1000, 1000))
    estimate <- terminalSpectrumFit(matrix(counts, nrow = B, ncol = n - 1L,
        byrow = TRUE) / B, mu, L)
    expect_equal(unname(estimate$estimate), c(ancestral, generations))
    expect_identical(fit$convergence, 0L)
    expect_equal(fit$par, c(ancestral, generations), tolerance = .001)
})

test_that("founder summaries reproduce their source copying records", {
    x <- founderDiversityExample(N = 2L, generations = 3L, replicates = 4L, seed = 91L)
    initial <- cbind(1:4, c(0L, 0L, 1L, 1L))
    record <- simulateForward(initial, rep(2L, 4L), seed = 91L)
    direct <- vapply(record$states, function(state) {
        counts <- tabulate(state[, 1L], 4L)
        1 - sum(counts * (counts - 1)) / 12
    }, 0)
    expect_equal(x$nonIBD[1L, ], direct)
    expect_equal(x$summary$expectation, .75^(0:3))
    expect_identical(x$truth$seeds, 91:94)
})

test_that("terminal-extension examples reuse source trees and mutation draws", {
    x <- terminalExtensionExample(regions = 4L, bootstrap = 3L)
    original <- pg_kingman_tree(80, 12L, pg_region_seed(814903L, 1L))
    expect_equal(pg_tree_summaries(x$truth$trees[[1L]])$TMRCA,
                 pg_tree_summaries(original)$TMRCA + 40)
    reference <- pg_mutate_tree(x$truth$trees[[1L]], 1000, 2e-5,
        pg_region_seed(914903L, 1L, TRUE))
    expect_identical(x$observed[[1L]], reference)
    expect_equal(x$spectra[1L, ], sfs(reference$haplotypes), ignore_attr = TRUE)
    expect_identical(x$fit, terminalSpectrumFit(x$spectra, 2e-5, 1000))
    expect_identical(dim(x$bootstrap), c(3L, 2L))
    boundary <- terminalSpectrumFit(matrix(0, 2L, 3L), 1e-5, 1000)
    expect_equal(unname(boundary$estimate), c(0, 0))
})

test_that("specified local trees generate compatible windows and one change", {
    x <- specifiedLocalTreeExample()
    expect_identical(x, specifiedLocalTreeExample())
    observed <- x$observed
    constraints <- pti_linked_constraints(observed$haplotypes, observed$positions)
    expect_identical(constraints$minimum_changes, 1L)
    expect_lt(max(constraints$witnesses$left), x$truth$breakpoint)
    expect_gt(min(constraints$witnesses$right), x$truth$breakpoint)
    for (j in 1:2) {
        limits <- if (j == 1) c(0, 25000) else c(75000, 100000)
        take <- observed$positions >= limits[1L] & observed$positions < limits[2L]
        G <- observed$haplotypes[, take, drop = FALSE]
        fit <- pti_fit_ages(G, x$calibration$mutationRate, diff(limits))
        expect_identical(fit$consensus$compatible_count, 1L)
        expect_true(is.finite(fit$best_loglik))
    }
})

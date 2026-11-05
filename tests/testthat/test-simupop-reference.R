test_that("source-aligned RNG and copying match the original simuPOP fixtures", {
    fixture <- readRDS(test_path("fixtures", "simupop-reference.rds"))
    expect_identical(fixture$rng$engine$version, "1.1.18")
    set.seed(29301L)
    before <- .Random.seed
    kind <- RNGkind()
    out <- simupopTestRng(fixture$rng)
    expect_identical(.Random.seed, before)
    expect_identical(RNGkind(), kind)
    expect_true(all(out$uniform_equal))
    expect_true(all(out$integer_equal))
    expect_true(all(out$geometric_equal))
    expect_true(all(out$copy_equal))
    # The controls demonstrate that a seed number is not a complete RNG state.
    expect_false(any(out$default_uniform_equal))
    expect_false(any(out$r_poisson_equal))
    expect_false(any(out$r_geometric_equal))
})

test_that("parental mutation timing has the correct complete sampling law", {
    fixture <- readRDS(test_path("fixtures", "simupop-reference.rds"))
    for (d in fixture$timing$results) {
        n <- 2 * d$N
        native <- as.numeric(unlist(d$values)) / n
        # This short regression complements the full 4096-run reference check.
        seeds <- as.integer(unlist(d$seeds))[seq_len(256L)] + 99991L
        observed <- vapply(seeds, function(seed)
            mean(simulatePopulation(c(d$N, d$N), 1L, d$u, d$N,
                reverseMutationRate = 0, seed = seed)$observed$haplotypes), 0)
        probability <- vapply(0:n, function(k)
            sum(dbinom(0:n, n, d$u) * dbinom(k, n, (0:n) / n)), 0)
        expectedVariance <- d$u * (1 - d$u) / n * (2 - 1 / n)
        wrongVariance <- d$u * (1 - d$u) / n
        expect_gt(expectedVariance, wrongVariance)
        expect_lt(abs(var(native) - expectedVariance), .004)
        bound <- sqrt(log(4000000) / (2 * length(observed)))
        empirical <- cumsum(tabulate(round(observed * n) + 1L, n + 1L) / length(observed))
        expect_lt(max(abs(empirical - cumsum(probability))), bound)
    }
})

test_that("ancestor construction reproduces the original C implementation", {
    cases <- readRDS(test_path("fixtures", "tsinfer-ancestors.rds"))
    for (fixture in cases) {
        actual <- hpgR:::.tsiGenerateAncestors(fixture$G)
        expected <- fixture$ancestors
        for (j in seq_along(expected))
            expect_equal(actual[[j]][names(expected[[j]])], expected[[j]], tolerance = 1e-14)
    }
})

test_that("dense copying HMM retains original C optimal path scores", {
    cases <- readRDS(test_path("fixtures", "tsinfer-matching.rds"))
    for (fixture in cases) {
        with(fixture, {
            fit <- hpgR:::.tsiMatch(haplotype, reference, nodeIds, rho, mu, numNodes)
            use <- fit$start:fit$end
            rows <- match(fit$parent, nodeIds)
            agree <- reference[cbind(rows, use)] == haplotype[use]
            score <- sum(ifelse(agree, log1p(-mu[use]), log(mu[use]))) +
                log(1 - rho[use[1L]] + rho[use[1L]] / numNodes)
            if (length(use) > 1L) score <- score +
                sum(ifelse(diff(fit$parent) == 0, log(1 - rho[use[-1L]] + rho[use[-1L]] / numNodes),
                           log(rho[use[-1L]] / numNodes)))
            expect_equal(score, c_score, tolerance = 1e-10)
        })
    }
})

test_that("recombining sequence inference reconstructs every observed allele", {
    for (seed in 1:6) {
        a <- simulateAncestry(2L, 50, 1000, .00005, seed = 4400L + seed)
        observed <- ancestryData(mutateAncestry(a, .0002, 5400L + seed))$observed
        inferred <- inferAncestry(observed$haplotypes, observed$positions,
                                  observed$sequenceLength)
        expect_identical(ancestryData(inferred)$observed$haplotypes, observed$haplotypes)
        expect_identical(ancestryData(inferred)$observed$positions, observed$positions)
        expect_equal(inferred$timeUnits, "ancestor ordering")
        expect_true(all(inferred$nodes$time[inferred$edges$parent + 1L] >
                        inferred$nodes$time[inferred$edges$child + 1L]))
        dated <- dateAncestry(inferred, .0002, 50)
        expect_true(all(is.finite(dated$nodes$time)))
        expect_equal(rowSums(dated$dating$posterior), rep(1, nrow(dated$nodes)))
        expect_identical(ancestryData(dated)$observed$haplotypes, observed$haplotypes)
    }
})

test_that("dating uses topological order, not initial numeric ancestor ages", {
    G <- rbind(c(1L, 1L, 0L, 0L), c(1L, 0L, 1L, 0L),
               c(0L, 0L, 0L, 1L), c(0L, 0L, 0L, 0L))
    x <- inferAncestry(G, c(100, 200, 300, 400), 1000)
    a <- dateAncestry(x, 1e-5, 100)
    x$nodes$time <- 10000 * x$nodes$time^2
    b <- dateAncestry(x, 1e-5, 100)
    expect_identical(a$dating, b$dating)
    expect_identical(a$nodes$time, b$nodes$time)
    large <- dateAncestry(x, 1e-5, 1000)
    expect_false(isTRUE(all.equal(a$dating$mean, large$dating$mean)))
})

test_that("coalescent prior moments include every ancestral waiting interval", {
    for (n in 2:20) {
        variance <- hpgR:::.td_tsdate_conditional_coalescent_variance(n)
        expected <- sum(4 / ((2:n)^2 * (1:(n - 1L))^2))
        expect_equal(variance[n + 1L], expected, tolerance = 1e-12)
    }
    expect_equal(hpgR:::.td_tsdate_conditional_coalescent_variance(3)[4], 10/9)
})

test_that("two-copy inside-outside posterior agrees with direct enumeration", {
    G <- rbind(c(1L, 1L, 0L), c(0L, 0L, 1L))
    x <- inferAncestry(G, c(10, 20, 30), 1000)
    dated <- dateAncestry(x, .00001, 100, timepoints = 80L)
    prior <- hpgR:::.td_build_prior_grid(x, 100, timepoints = 80L)
    root <- setdiff(x$nodes$id, x$sampleIds)
    logWeight <- log(hpgR:::.td_ntv_get(prior, root)) +
        dpois(3L, 2 * .00001 * 1000 * (prior$timepoints + 1e-8), log = TRUE)
    weight <- exp(logWeight - max(logWeight))
    weight <- weight / sum(weight)
    expect_equal(dated$dating$posterior[root + 1L, ], weight, tolerance = 1e-10)
})

test_that("inference and dating validate their observed inputs", {
    expect_error(inferAncestry(matrix(0L, 2L, 2L), c(1, 2), 3), "polymorphic")
    expect_error(inferAncestry(matrix(2L, 2L, 2L), c(1, 2), 3), "binary")
    expect_error(inferAncestry(matrix(c(0L,1L), 2L, 2L), c(2, 1), 3), "ordered")
    x <- inferAncestry(matrix(c(0L, 1L), 2L, 1L), 1, 3)
    expect_error(dateAncestry(x, 0, 100), "mutationRate")
    expect_error(dateAncestry(x, .0001, 100, 2L), "timepoints")
})

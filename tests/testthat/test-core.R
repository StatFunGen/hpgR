test_that("rate-map calculations retain interval boundaries and units", {
    x <- rateMap(c(0, 10, 30), c(0.01, 0.02))
    expect_identical(rateAt(x, c(0, 5, 10, 30)), c(0.01, 0.01, 0.02, 0.02))
    expect_equal(cumulativeRate(x, c(0, 5, 10, 20, 30)), c(0, .05, .1, .3, .5))
    expect_equal(meanRate(x), .5 / 30)
    expect_error(rateMap(c(0, 0), 1))
    expect_error(rateAt(x, -1))
    expect_error(cumulativeRate(x, NA_real_))
})

test_that("the expected SFS uses diploid Ne and haploid sample counts", {
    for (n in c(2L, 3L, 4L, 12L, 20L)) {
        expect_equal(expectedSFS(n, 10000, 1e-8, 10000),
                     4 / seq_len(n - 1L), tolerance = 1e-11)
        expect_equal(expectedSFS(n, 10000, 0, 10000), rep(0, n - 1L))
    }
    expect_error(expectedSFS(1, 10000, 1e-8, 10000))
    expect_error(expectedSFS(4, 0, 1e-8, 10000))
})

test_that("zero-based bins and empty bins follow the original header", {
    expect_identical(assignCMBins(c(-.1, 0, .25, .5, .75, 1), 0, .25, 4),
                     c(-1L, 0L, 1L, 2L, 3L, -1L))
    expect_identical(assignCMBins(c(0, 1), 0, 0, 2), c(-1L, -1L))
    expect_identical(pairCountHistogram(c(2L, 0L, 3L)), c(13, 0, 6))
    G <- matrix(c(0, 1, 2, 0), nrow = 2L)
    out <- weightedLDBins(G, c(.75, .75), c(.25, .25),
                          c(.25, .5), c(0L, 2L), 3L)
    expect_identical(out$snps_per_bin, c(1L, 0L, 1L))
    expect_identical(out$pair_count, c(2, 0, 1))
    expect_identical(out$S[, 2L], c(0, 0))
})

test_that("parental copying preserves exact alleles and dimensions", {
    initial <- matrix(0:5, nrow = 3L, byrow = TRUE)
    p <- list(c(3L, 1L, 1L), c(2L, 1L))
    states <- replayCopying(initial, p)
    expect_identical(states[[1]], initial)
    expect_identical(states[[2]], initial[p[[1]], , drop = FALSE])
    expect_identical(states[[3]], initial[p[[1]][p[[2]]], , drop = FALSE])
    expect_identical(dim(copyGenotypes(matrix(1L, 1, 1), 1L)), c(1L, 1L))
    expect_error(copyGenotypes(initial, 0L))
    expect_error(copyGenotypes(initial, 1.5))
    expect_error(copyGenotypes(matrix(.5, 1L), 1L))
})

test_that("posterior summaries retain upstream calculations", {
    x <- matrix(1:12, ncol = 3L)
    actual <- posteriorSummary(x)
    expect_equal(actual$mean, colMeans(x))
    expect_equal(actual$median, apply(x, 2L, median))
    expect_equal(actual$lower, apply(x, 2L, quantile, .025))
})

test_that("source records ship with the installed package", {
    sources <- sourceProvenance()
    expect_true(all(c("rate-map", "sfs", "posterior", "weighted-ld",
                      "parent-choice", "transmission", "parent-records",
                      "native-mating", "native-mutation", "native-sampling") %in% sources$id))
    expect_false(anyDuplicated(sources$id) > 0L)
    expect_true(all(nchar(sources$sha256) == 64L))
})

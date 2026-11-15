test_that("deme rows select the parental population before finite copying", {
    swap <- matrix(c(0, 1, 1, 0), nrow = 2)
    p <- matrix(c(0, 1, 1, 0), nrow = 2)
    expect_identical(simulateDemeFrequencies(p, 50, swap, 0), p)
    expect_identical(simulateDemeFrequencies(p, 50, swap, 1), p[2:1, ])
    expect_error(simulateDemeFrequencies(p, 50, swap / 2, 1), "row sums")

    set.seed(814902)
    m <- matrix(c(.75, .25, .25, .75), nrow = 2, byrow = TRUE)
    result <- simulateDemeFrequencies(p[, rep(1, 20000)], 50, m, 1)
    expect_equal(rowMeans(result), c(.25, .75), tolerance = .003)
})

test_that("ancestral pairs must occupy the same deme before they can merge", {
    set.seed(814903)
    expect_equal(simulateCrossPopulationTimes(10, rep(0, 8), 1, 1), rep(9, 10))
    expect_equal(simulateCrossPopulationTimes(10, numeric(0), 1, 1), rep(1, 10))
    times <- simulateCrossPopulationTimes(30000, .5, 1, 1)
    expect_lt(abs(mean(times == 1) - .5), .02)
})

test_that("the split gives equal discordant counts and the pulse favors ABBA", {
    set.seed(814904)
    split <- simulateAbbaBaba(100000, 0, 1000, 1000, 1600, 300, .0002)
    pulse <- simulateAbbaBaba(100000, .2, 1000, 1000, 1600, 300, .0002)
    expect_lt(abs(unname(diff(split))) / sum(split), .05)
    expect_gt(unname(pulse["ABBA"] / pulse["BABA"]), 1.35)
})

test_that("f4 uses the same allele orientation in all four samples", {
    p <- cbind(c(.8, .6), c(.3, .9), c(.6, .2), c(.2, .7))
    observed <- function(p) do.call(frequencyF4, unname(as.list(as.data.frame(p))))
    f4 <- mean((p[, 1] - p[, 2]) * (p[, 3] - p[, 4]))
    expect_equal(observed(p), f4, tolerance = 1e-14)
    expect_equal(observed(1 - p), f4, tolerance = 1e-14)
    expect_equal(observed(p[, c(2, 1, 3, 4)]), -f4, tolerance = 1e-14)
})

test_that("forward diversity agrees with the exact allele-count process", {
    # Enumerate parental mutations, then every binomial offspring draw.
    sizes <- c(2, 3, 5)
    mu <- .03
    count <- 0:(2 * sizes[1L])
    distribution <- dbinom(count, 2 * sizes[1L], .3)
    diversity <- function(probabilities) {
        n <- length(probabilities) - 1L
        k <- 0:n
        sum(probabilities * 2 * k * (n - k) / (n * (n - 1)))
    }
    exact <- diversity(distribution)
    for (t in 1:2) {
        n <- 2 * sizes[t]
        child <- 2 * sizes[t + 1L]
        nextDistribution <- numeric(child + 1L)
        for (k in 0:n) for (u in 0:k) for (v in 0:(n - k)) {
            weight <- distribution[k + 1L] * dbinom(u, k, 1 - mu) *
                dbinom(v, n - k, mu)
            nextDistribution <- nextDistribution + weight *
                dbinom(0:child, child, (u + v) / n)
        }
        distribution <- nextDistribution
        expect_equal(sum(distribution), 1, tolerance = 1e-12)
        exact <- c(exact, diversity(distribution))
    }
    expect_equal(expectedDiversity(sizes, mu, initialDiversity = .42),
                 exact, tolerance = 1e-12)
    expect_equal(expectedDiversity(sizes, 0, .42),
                 .42 * c(1, cumprod(1 - 1 / (2 * head(sizes, -1L)))))
})

test_that("ancestral size inversion uses the declared recent and ancestral histories", {
    mu <- 2.5e-5
    sizes <- c(rep(20, 20), 20 * 50^((0:40) / 40))
    for (N in c(40, 100, 500)) {
        a <- 2 * mu * (1 - mu)
        q <- (1 - 1 / (2 * N)) * (1 - 2 * mu)^2
        equilibrium <- (1 - 1 / (2 * N)) * a / (1 - q)
        final <- tail(expectedDiversity(sizes, mu, equilibrium), 1L)
        fit <- estimateAncestralSize(final, sizes, mu)
        expect_equal(fit$Ne, N, tolerance = 1e-9)
        expect_equal(fit$mutationInput + fit$ancestralWeight * equilibrium,
                     final, tolerance = 1e-12)
        burn <- tail(expectedDiversity(rep(N, 1001), mu), 1L)
        finite <- tail(expectedDiversity(sizes, mu, burn), 1L)
        expect_equal(estimateAncestralSize(finite, sizes, mu,
                     ancestralGenerations = 1000)$Ne, N, tolerance = 1e-8)
    }
    expect_error(estimateAncestralSize(0, sizes, mu), "outside the range")
    expect_error(estimateAncestralSize(.9, sizes, mu), "outside the range")
    expect_error(estimateAncestralSize(.01, sizes, 0), "positive mutation")
    expect_error(expectedDiversity(sizes, .5), "smaller than")
})

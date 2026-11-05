test_that("final-sample diversity agrees with explicit sequence-pair enumeration", {
    x <- rbind(c(0, 0, 1, 0), c(0, 1, 1, 0),
               c(1, 1, 1, 0), c(1, 1, 0, 0))
    pairs <- utils::combn(nrow(x), 2L)
    direct <- mean(apply(pairs, 2L, function(pair) sum(x[pair[1L], ] != x[pair[2L], ])))
    summary <- summarizePopulationSample(x)
    expect_equal(summary$pairwiseDifferences, direct)
    expect_equal(summary$piPerSite, direct / 4)
    expect_identical(unname(summary$foldedSfs), c(2L, 1L))
    expect_equal(sum(summary$groupedSfs), 3)
    expect_equal(summary$segregatingSites, 3)
    expect_equal(estimateDiversitySize(x, 1e-5)$Ne, direct / (4 * 4e-5))
    # All-one columns and omitted invariant columns make no differences.
    expect_equal(summarizePopulationSample(x, sequenceLength = 8)$piPerSite, direct / 8)
    expect_equal(summarizePopulationSample(1L - x)$foldedSfs, summary$foldedSfs)
    expect_error(summarizePopulationSample(matrix(2L, 3, 4)), "0 and 1")
})

test_that("chromosome summaries retain lengths and missing invariant bases", {
    x <- rbind(c(0, 0, 1, 0), c(0, 1, 1, 0),
               c(1, 1, 1, 0), c(1, 1, 0, 0))
    chromosome <- c("a", "a", "b", "b")
    summary <- summarizePopulationSample(x, chromosome, c(a = 10L, b = 20L))
    expect_equal(summary$sequenceLength, 30)
    expect_equal(sum(vapply(summary$byChromosome, `[[`, numeric(1),
                            "pairwiseDifferences")), summary$pairwiseDifferences)
    expect_equal(summary$byChromosome$a$sequenceLength, 10)
    expect_equal(summarizePopulationSample(x, chromosome)$sequenceLength, 4)
    expect_error(summarizePopulationSample(x, chromosome, 30), "named analyzed")
    expect_error(summarizePopulationSample(x, chromosome, c(a = 1, b = 2)), "smaller")
    expect_error(summarizePopulationSample(x, sfsGroups = list(a = 1, b = 1)), "nonoverlapping")
    expect_equal(summarizePopulationSample(x[1:2, , drop = FALSE])$segregatingSites, 1)
})

test_that("the symmetric parental mutation inverse solves its probability recursion", {
    # Exact expectations are rational when u=1/4 and N=2. Construct a
    # sequence sample with pi=H=9/26, then verify the inverse independently.
    u <- .25
    a <- 2 * u * (1 - u)
    b <- (1 - 2 * u)^2
    Ne <- 2
    c <- 1 / (2 * Ne)
    H <- (1 - c) * a / (1 - (1 - c) * b)
    expect_equal(H, 9 / 26)
    # n=2 gives one difference per differing site, so pi=9/26 exactly.
    x <- rbind(rep(0L, 26), c(rep(1L, 9), rep(0L, 17)))
    result <- estimateDiversitySize(x, u, mutationModel = "symmetric_parental")
    expect_equal(result$piPerSite, H)
    expect_equal(result$Ne, Ne)
    expect_equal((1 - c) * (a + b * H), H)
    expect_equal(estimateDiversitySize(matrix(0L, 4, 5), 1e-5)$Ne, 0)
    expect_error(estimateDiversitySize(x, .5), "smaller")
    saturated <- rbind(rep(0L, 4), rep(1L, 4))
    expect_identical(estimateDiversitySize(saturated, u,
        mutationModel = "symmetric_parental")$Ne, Inf)
})

test_that("parental mutation diversity matches an independently enumerated WF chain", {
    # Build the allele-count Markov chain without using the diversity recursion.
    # First flip each of the M parental copies, then binomially copy the
    # realized mutated parental pool into M offspring positions.
    M <- 4L
    u <- .25
    states <- 0:M
    transition <- matrix(0, M + 1L, M + 1L)
    for (k in states) {
        mutation <- numeric(M + 1L)
        for (retained in 0:k) for (new in 0:(M - k)) {
            j <- retained + new
            mutation[j + 1L] <- mutation[j + 1L] +
                stats::dbinom(retained, k, 1 - u) * stats::dbinom(new, M - k, u)
        }
        for (j in states) transition[k + 1L, ] <- transition[k + 1L, ] +
            mutation[j + 1L] * stats::dbinom(states, M, j / M)
    }
    expect_equal(rowSums(transition), rep(1, M + 1L))
    constraint <- t(transition) - diag(M + 1L)
    constraint[M + 1L, ] <- 1
    stationary <- solve(constraint, c(rep(0, M), 1))
    H <- sum(stationary * 2 * states * (M - states) / (M * (M - 1)))
    a <- 2 * u * (1 - u)
    b <- (1 - 2 * u)^2
    expected <- (1 - 1 / M) * a / (1 - (1 - 1 / M) * b)
    expect_equal(H, expected)
    expect_equal((a + b * H) / (2 * a * (1 - 2 * H)), M / 2)
    rawH <- sum(stationary * 2 * states * (M - states) / M^2)
    expect_equal(rawH, (1 - 1 / M) * H)
})

test_that("LD inversion retains each pair's map fraction and direct correlation", {
    x <- rbind(c(0, 0, 0), c(0, 0, 1), c(0, 1, 0), c(1, 1, 1),
               c(1, 1, 1), c(1, 1, 0), c(0, 0, 0), c(1, 1, 1))
    position <- c(0, .02, .09)
    result <- estimateLdSize(x, position, distanceBreaks = c(.001, .1), seed = 9L)
    expected <- vapply(seq_len(nrow(result$pairs)), function(j) {
        pair <- result$pairs[j, ]
        stats::cor(x[, pair$first], x[, pair$second])^2
    }, numeric(1))
    expect_equal(result$pairs$rho2, expected)
    r <- (1 - exp(-2 * result$pairs$distanceMorgan)) / 2
    gamma <- ((1 - r)^2 + r^2) / (2 * r * (2 - r))
    expect_equal(result$pairs$gamma, gamma)
    adjusted <- mean(expected) - 1 / nrow(x)
    expect_equal(result$Ne, mean(gamma) / adjusted)
    expect_equal(result$bins$Ne, result$Ne)
    expect_equal(result$bins$pairs, 3)
    # Exact Haldane parity across two adjacent intervals.
    adjacent <- (1 - exp(-2 * diff(position))) / 2
    expect_equal((1 - prod(1 - 2 * adjacent)) / 2,
                 (1 - exp(-2 * diff(range(position)))) / 2)
})

test_that("LD noise and cross-chromosome pairs cannot create a finite estimate", {
    x <- rbind(c(0, 0), c(0, 1), c(1, 0), c(1, 1))
    result <- estimateLdSize(x, c(0, .1), distanceBreaks = c(.01, .2))
    expect_equal(result$pooled["correctedRho2"], c(correctedRho2 = -.25))
    expect_identical(result$Ne, Inf)
    expect_error(estimateLdSize(x, c(0, .1), chromosome = c("a", "b")), "within-chromosome")
    expect_error(estimateLdSize(x, c(0, .1), distanceBreaks = c(.1, .01)), "increasing")
})

test_that("pair subsampling is deterministic and restores the caller's RNG state", {
    set.seed(91)
    x <- matrix(rbinom(40 * 50, 1, .5), 40, 50)
    state <- .Random.seed
    a <- estimateLdSize(x, seq(0, .5, length.out = 50), maxPairs = 50L, seed = 3L)
    expect_identical(.Random.seed, state)
    b <- estimateLdSize(x, seq(0, .5, length.out = 50), maxPairs = 50L, seed = 3L)
    expect_identical(a, b)
    expect_lte(nrow(a$pairs), 50)
    # Exercise the bounded-memory pair sampler independently of LD values.
    set.seed(1)
    pairs <- hpgR:::inferencePairIndices(2000L, 1000L)
    expect_equal(nrow(pairs), 1000)
    expect_true(all(pairs[, 1] < pairs[, 2]))
    expect_false(anyDuplicated(pairs) > 0)
})

test_that("candidate comparison reads final summaries and independent calibration", {
    training <- list(stable = rbind(c(1, 2), c(2, 3), c(3, 4)),
                     bottleneck = rbind(c(5, 6), c(6, 7), c(7, 8)))
    calibration <- list(stable = rbind(c(1, 2), c(2, 3), c(3, 4)),
                        bottleneck = rbind(c(5, 6), c(6, 7), c(7, 8)))
    training <- lapply(training, function(x) {colnames(x) <- c("rare", "common"); x})
    calibration <- lapply(calibration, function(x) {colnames(x) <- c("rare", "common"); x})
    library <- buildHistoryLibrary(training, calibration)
    observed <- c(common = 3, rare = 2)
    fit <- comparePopulationHistories(observed, library)
    expect_identical(fit$selected, "stable")
    expect_equal(fit$profile$discrepancy, c(0, 32))
    expect_equal(fit$profile$compatibility, c(1, .25))
    expect_identical(fit$observedSummary, c(rare = 2, common = 3))
    expect_equal(library$centers["bottleneck", ], c(rare = 6, common = 7))
    expect_error(comparePopulationHistories(c(rare = 1), library), "match")
    expect_error(buildHistoryLibrary(training, calibration[-1]), "same named")
    # A positive variance floor handles summaries that are constant in a bank.
    flat <- lapply(training, function(x) {x[, ] <- 1; x})
    expect_true(all(buildHistoryLibrary(flat, flat)$scale > 0))
})

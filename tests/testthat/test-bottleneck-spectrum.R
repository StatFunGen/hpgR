test_that("profile cells match direct source-kernel likelihoods", {
    spectra <- rbind(c(30, 12, 8), c(28, 14, 7))
    fit <- fitBottleneckSpectrum(spectra, 1e-5, 10000,
        bottleneckSizes = c(10, 20), growthGenerations = c(20, 40))
    for (i in seq_len(nrow(fit$profile))) {
        row <- fit$profile[i, ]
        growth <- round(exp(seq(log(row$bottleneckSize), log(1000),
            length.out = row$growthGenerations + 1L)))[-1L]
        parents <- head(c(row$bottleneckSize, rep(row$bottleneckSize, 20), growth), -1L)
        expected <- expectedHistorySFS(4, c(row$ancestralSize, parents), 1e-5, 20000)
        expect_equal(row$logLikelihood,
            sum(dpois(colSums(spectra), expected, log = TRUE)), tolerance = 1e-9)
    }
    row <- fit$profile[which.max(fit$profile$logLikelihood), ]
    expect_equal(fit$estimate[["bottleneckBegan"]], row$growthGenerations + 20)
    expect_equal(fit$estimate[["bottleneckSize"]], row$bottleneckSize)
    expect_true(all(fit$profile$ancestralSize >= fit$profile$bottleneckSize))
})

test_that("a spectrum generated under a declared history recovers that history", {
    growth <- round(20 * 50^((1:40) / 40))
    parents <- head(c(20, rep(20, 20), growth), -1L)
    expected <- expectedHistorySFS(12, c(100, parents), 1e-5, 1e7)
    spectra <- rbind(round(expected), round(expected))
    fit <- fitBottleneckSpectrum(spectra, 1e-5, 1e7,
        bottleneckSizes = c(10, 20, 40), growthGenerations = c(20, 40, 60))
    expect_equal(fit$estimate[["bottleneckSize"]], 20)
    expect_equal(fit$estimate[["growthGenerations"]], 40)
    expect_equal(fit$estimate[["ancestralSize"]], 100, tolerance = .002)
})

test_that("whole-chromosome refits preserve the caller RNG and return date intervals", {
    spectra <- rbind(c(30, 12, 8), c(28, 14, 7), c(20, 13, 9))
    set.seed(8409)
    before <- .Random.seed
    args <- list(spectra = spectra, mutationRate = 1e-5, sequenceLength = 10000,
        bottleneckSizes = c(10, 20), growthGenerations = c(20, 40),
        bootstrapReplicates = 5L, seed = 321L)
    fit <- do.call(fitBottleneckSpectrum, args)
    expect_identical(.Random.seed, before)
    expect_identical(fit, do.call(fitBottleneckSpectrum, args))
    expect_equal(fit$bootstrap[, "bottleneckBegan"],
                 fit$bootstrap[, "growthGenerations"] + 20)
    expect_equal(fit$interval,
        apply(fit$bootstrap, 2L, quantile, c(.025, .975)))
    expect_error(fitBottleneckSpectrum(spectra, 1e-5, 10000,
        bottleneckSizes = 1000), "Bottleneck sizes")
    expect_error(fitBottleneckSpectrum(spectra / 7, 1e-5, 10000), "integer")
    expect_error(fitBottleneckSpectrum(spectra[1, , drop = FALSE], 1e-5, 10000,
        bootstrapReplicates = 5), "two chromosomes")
})

test_that('classroom interfaces reject invalid values and length recycling', {
    expect_error(hweProbability(NA_real_), 'finite')
    expect_error(hweProbability(1+1i), 'finite')
    expect_error(frequencyKinship(0:2, 0:2, .5), 'length 3')
    expect_error(kingKinship(0:2, 1), 'same length')
    expect_error(standardizeGenotypes(matrix(0, 2, 2)), 'strictly between')
    expect_error(genotypePCA(matrix(0, 5, 4)), 'variation')
    expect_error(fitGaussianMixture(numeric(), list(pi=.5, mu=c(-1,1))), 'finite')
    expect_error(fitGenotypeError(0:2, max_iterations=0), 'integer')
    expect_error(copyingTransition(.8, c(.4,.4)), 'sum to 1')
    expect_error(binaryHMM(c(0,1), c(.1,.2), .1), 'length 1')
    expect_error(selectedFrequency(.5, c(0,0,0)), 'positive fitness')
    expect_error(selectedFrequency(1, c(0,1,1)), 'mean fitness')
    expect_error(simulateFrequencyRecords(N=1.5), 'integer')
    expect_error(copyHistory(c(2,2), 1:3), 'length 4')
})

test_that('single-person and single-marker simulations retain matrix dimensions', {
    set.seed(123)
    expect_equal(dim(simulateGenotypes(1, c(.2,.5))), c(1L,2L))
    expect_equal(dim(simulateGenotypes(3, .2)), c(3L,1L))
    expect_equal(dim(simulateAncestryAssociation(1,2,1)$local), c(1L,2L))
    expect_equal(dim(simulateAncestryAssociation(3,1,1)$local), c(3L,1L))
    fit <- fitAdmixture(matrix(c(0,1,2),3,1),1,matrix(1:3,3,1))
    expect_equal(dim(fit$p), c(1L,1L))
    expect_equal(fit$q, matrix(1,3,1))
})

test_that('log likelihood and probability bounds remain stable at large counts', {
    expected <- (9000 * log(.9) + 1000 * log(.1) - 10000 * log(.5)) / log(10)
    expect_equal(linkageLOD(.1,9000,1000), expected, tolerance=1e-12)
    expect_identical(linkageLOD(0,0,0), 0)
    expect_identical(linkageLOD(0,0,10), -Inf)
    expect_gt(zeroCountUpperBound(1e20), 0)
    expect_equal(zeroCountUpperBound(1e20) * 1e20, -log(.05))
    expect_true(is.na(kingKinship(c(0,2), c(0,2))))
    expect_true(is.na(pairwiseFst(0,0)))
})

test_that('returned marker effects generate the returned genetic values', {
    set.seed(42)
    Z <- matrix(rnorm(100),20,5)
    trait <- simulatePolygenicTrait(Z,.5,.5)
    expect_equal(drop(Z %*% trait$marker_effect), trait$genetic_value)
    expect_equal(trait$marker_effect, trait$drawn_effect * trait$effect_scale)
    expect_equal(var(trait$genetic_value), .5)
    zero <- simulatePolygenicTrait(Z,0,.5)
    expect_equal(zero$genetic_value, rep(0,20))
    expect_equal(zero$marker_effect, rep(0,5))
})

test_that('seeded demonstrations restore the caller RNG and generator kind', {
    restore <- hpgR:::rngRestorer()
    on.exit(restore(), add=TRUE)
    RNGkind('L\'Ecuyer-CMRG')
    set.seed(73)
    before <- .Random.seed
    kind <- RNGkind()
    for (draw in list(samplingExample, populationProcessExample, coalescentExample,
        function() simulateFrequencyRecords(generations=2),
        function() simulateNewMutationPaths(generations=2))) {
        invisible(draw())
        expect_identical(.Random.seed, before)
        expect_identical(RNGkind(), kind)
    }
    rm('.Random.seed', envir=.GlobalEnv)
    invisible(samplingExample())
    expect_false(exists('.Random.seed', envir=.GlobalEnv, inherits=FALSE))
    expect_error(hpgR:::pgWithSeed(1, function() stop('draw failed')), 'draw failed')
    expect_false(exists('.Random.seed', envir=.GlobalEnv, inherits=FALSE))
})

test_that('EM reports an exhausted iteration limit', {
    expect_warning(x <- fitGaussianMixture(c(-2,-1,1,2),
        list(pi=.3,mu=c(-.5,.5)), max_iterations=1), 'max_iterations')
    expect_false(x$converged)
    expect_equal(x$iterations,1L)
    expect_warning(x <- fitGenotypeError(c(0,0,0,1,2), max_iterations=1), 'max_iterations')
    expect_false(x$converged)
})

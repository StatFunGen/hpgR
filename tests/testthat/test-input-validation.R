test_that("named checks accept defaults, bounds and dependent dimensions", {
    checkInputs <- hpgR:::checkInputs
    scalarInteger <- hpgR:::scalarInteger
    probability <- hpgR:::probability
    numericMatrix <- hpgR:::numericMatrix
    validate <- function(n, frequency, genotype) {
        checkInputs(n = scalarInteger, frequency = list(probability, length = n),
                    genotype = list(numericMatrix, columns = n))
    }
    expect_invisible(validate(2, c(0, 1), matrix(0, 1, 2)))
    expect_error(validate(2.5, c(.2, .8), matrix(0, 1, 2)), "n must be an integer")
    expect_error(validate(2, .2, matrix(0, 1, 2)), "frequency.*length 2")
    expect_error(validate(2, c(.2, .8), matrix(0, 1, 1)), "genotype.*2 column")
    expect_error(validate(2, c(.2, Inf), matrix(0, 1, 2)), "frequency.*finite")
})

test_that("checks leave argument values and the random stream unchanged", {
    checkInputs <- hpgR:::checkInputs
    scalarInteger <- hpgR:::scalarInteger
    genotypeMatrix <- hpgR:::genotypeMatrix
    validate <- function(n, genotype) {
        before <- list(n = n, genotype = genotype)
        checkInputs(n = scalarInteger, genotype = genotypeMatrix)
        expect_identical(list(n = n, genotype = genotype), before)
    }
    set.seed(8105)
    before <- .Random.seed
    validate(2, matrix(c(0, 1), 1, 2))
    expect_identical(.Random.seed, before)
    expect_identical(typeof(2), "double")
})

test_that("validation stops in declaration order before later defaults are forced", {
    checkInputs <- hpgR:::checkInputs
    numericMatrix <- hpgR:::numericMatrix
    probability <- hpgR:::probability
    validate <- function(genotype, frequency = stop("later default evaluated")) {
        checkInputs(genotype = numericMatrix,
                    frequency = list(probability, length = ncol(genotype)))
    }
    expect_error(validate(NULL), "genotype.*numeric matrix")
    expect_error(validate(matrix(0, 1, 1)), "later default evaluated")
})

test_that("invalid declarations identify the check or argument", {
    checkInputs <- hpgR:::checkInputs
    scalarInteger <- hpgR:::scalarInteger
    n <- 2
    expect_invisible(checkInputs())
    expect_error(checkInputs(scalarInteger), "unique argument names")
    expect_error(checkInputs(n = scalarInteger, n = scalarInteger), "unique argument names")
    expect_error(checkInputs(n = 1), "invalid input check for n")
    expect_error(checkInputs(n = list()), "invalid input check for n")
    expect_error(checkInputs(n = list(1)), "invalid input check for n")
    expect_error(checkInputs(unknownArgument = scalarInteger), "unknownArgument")
})

test_that("vector rules distinguish scalar, integer and probability constraints", {
    checkInputs <- hpgR:::checkInputs
    scalarInteger <- hpgR:::scalarInteger
    numericVector <- hpgR:::numericVector
    probability <- hpgR:::probability
    validate <- function(n, p, positions, weights) {
        checkInputs(n = list(numericVector, minimum = 1, integer = TRUE),
                    p = list(probability, open = TRUE),
                    positions = list(numericVector, nondecreasing = TRUE),
                    weights = list(probability, length = NULL, sumToOne = TRUE))
    }
    expect_invisible(validate(c(1, 1e20), .5, c(0, 0, 1), c(.2, .8)))
    expect_error(validate(1.5, .5, 0:1, c(.2, .8)), "n.*integers")
    expect_error(validate(1, 0, 0:1, c(.2, .8)), "p.*strictly between")
    expect_error(validate(1, 1, 0:1, c(.2, .8)), "p.*strictly between")
    expect_error(validate(1, .5, 1:0, c(.2, .8)), "positions.*nondecreasing")
    expect_error(validate(1, .5, 0:1, c(.2, .7)), "weights.*sum to 1")
    for (n in list(NULL, numeric(), NA_real_, NaN, Inf, 1 + 1i, "2", TRUE, c(1, 2))) {
        expect_error(checkInputs(n = scalarInteger), "n must be an integer")
    }
})

test_that("binary rules preserve missing alleles and reject dimensions and recoding", {
    checkInputs <- hpgR:::checkInputs
    binaryVector <- hpgR:::binaryVector
    binaryMatrix <- hpgR:::binaryMatrix
    logicalScalar <- hpgR:::logicalScalar
    numericMatrix <- hpgR:::numericMatrix
    validate <- function(reference, observed, frequency, log) {
        checkInputs(reference = binaryMatrix,
                    observed = list(binaryVector, length = nrow(reference), missing = TRUE),
                    frequency = list(numericMatrix, minimum = 0, maximum = 1),
                    log = logicalScalar)
    }
    reference <- matrix(c(0, 1, 1, 0), 2, 2)
    expect_invisible(validate(reference, c(NA, 1), reference, FALSE))
    expect_error(validate(reference, c(TRUE, FALSE), reference, FALSE), "observed.*binary allele")
    expect_error(validate(reference, c(0, 2), reference, FALSE), "observed.*0 or 1")
    expect_error(validate(reference, matrix(0:1, 2), reference, FALSE), "observed.*binary allele")
    expect_error(validate(reference + .1, c(0, 1), reference, FALSE), "reference.*0 or 1")
    expect_error(validate(reference, c(0, 1), reference + .1, FALSE), "frequency.*values in")
    for (log in list(NA, logical(), c(TRUE, FALSE), 1, "TRUE"))
        expect_error(validate(reference, c(0, 1), reference, log), "log must be TRUE or FALSE")
    observed <- c(TRUE, NA)
    expect_invisible(checkInputs(observed = list(binaryVector, missing = TRUE, logical = TRUE)))
    observed <- c(0, NA)
    expect_error(checkInputs(observed = binaryVector), "observed.*finite")
})

test_that("failed public input checks precede simulation draws", {
    set.seed(8149)
    before <- .Random.seed
    expect_error(simulateGenotypes(3, c(.2, NA)), "frequency.*finite")
    expect_identical(.Random.seed, before)
    expect_error(simulateSiblingTraits(3, .5, 0, -1), "V_E.*nonnegative")
    expect_identical(.Random.seed, before)
    expect_error(simulateCopyingPath(matrix(.5, 2, 2), .1, .01), "reference.*0 or 1")
    expect_identical(.Random.seed, before)
})

test_that("HMM input checks retain singleton and varying-transition models", {
    one <- matrix(c(.2, .8), 1, 2)
    constant <- matrix(c(.9, .1, .1, .9), 2, 2)
    empty <- array(numeric(), c(2, 2, 0))
    expect_equal(hmmPosterior(one, empty)$posterior, hmmPosterior(one, constant)$posterior)
    emission <- rbind(one, rev(one), one)
    varying <- array(rep(constant, 2), c(2, 2, 2))
    expect_equal(hmmPosterior(emission, varying), hmmPosterior(emission, constant))
    expect_error(hmmPosterior(emission, array(1, c(2, 2, 1))), "transition dimensions")
    expect_error(hmmPosterior(emission, matrix(.4, 2, 2)), "transition rows")
    expect_error(hmmPosterior(emission, constant, c(.2, .7)), "initial.*sum")
    expect_error(hmmPosterior(emission, matrix(NA_real_, 2, 2)), "transition dimensions")
    expect_error(hmmPosterior(emission, matrix("a", 2, 2)), "transition dimensions")
    expect_error(hmmPosterior(emission * -1, constant), "emission.*values in")
    expect_error(hmmPosterior(matrix(numeric(), 0, 2), constant), "emission.*numeric matrix")
})

test_that("genealogy parameters use the shared scalar and index rules", {
    set.seed(8105L)
    before <- .Random.seed
    expect_error(pgForwardSfs(4, length_bp = NA_real_), "length_bp.*finite")
    expect_error(pgForwardSfs(4, mu = c(.1, .2)), "mu.*finite")
    expect_error(pgRecordAncestry(4, 100, 4, map_length = -1), "map_length.*nonnegative")
    expect_error(pgRecordAncestry(4, 100, 4, sample_indices = c(1, 1, 2, 3)),
                 "sample_indices.*distinct")
    expect_error(pgForwardExposure(4, 4, sample_indices = 6:9), "sample_indices.*values in")
    expect_error(pgKingmanTree(0, 4), "Ne.*positive")
    expect_error(pgDiscreteGenealogy(4, 4, max_extension = c(1, 2)), "max_extension.*length 1")
    expect_error(pgRecordAncestry(4, 100, 4, complete_roots = NA), "complete_roots.*TRUE or FALSE")
    expect_identical(.Random.seed, before)
    expect_identical(sum(pgForwardSfs(4, n = 4, mu = 0)$sfs), 0L)
})

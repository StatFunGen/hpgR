test_that("population simulation checks its units and dimensions", {
    expect_error(simulatePopulation(0, 2, 0, 1), "sizes")
    expect_error(simulatePopulation(10, c(2, 0), 0, 1), "loci")
    expect_error(simulatePopulation(10, 2, 2, 1), "probabilities")
    expect_error(simulatePopulation(10, 2, 0, 11), "sampleSize")
    expect_error(simulatePopulation(10, 2, 0, 1, founderSize = 11), "founderSize")
    expect_error(simulatePopulation(10, 2, 0, 1,
        geneticPositions = list(c(.1, 0))), "geneticPositions")
    expect_error(simulatePopulation(10, 2, 0, 1,
        positions = list(c(1, 1))), "positions")
})

test_that("native population records are separate from final observations", {
    args <- list(sizes = c(12, 4, 8), loci = c(3, 2), mutationRate = .01,
        sampleSize = 4, geneticPositions = list(c(0, .1, .2), c(0, .1)),
        burnin = 5, founderSize = 4, seed = 8149)
    x <- do.call(simulatePopulation, args)
    expect_identical(x, do.call(simulatePopulation, args))
    coarse_args <- args
    coarse_args$recordEvery <- 2L
    expect_identical(x$observed, do.call(simulatePopulation, coarse_args)$observed)
    expect_identical(dim(x$observed$haplotypes), c(8L, 5L))
    expect_identical(dim(x$observed$genotypes), c(4L, 5L))
    expect_identical(dim(x$records$allele_frequency), c(3L, 5L))
    expect_equal(x$records$generation$size, c(4, 4, 8))
    expect_equal(x$records$generation$generation, 0:2)
    expect_equal(x$records$burnin$generation, -5:0)
    expect_length(unique(x$truth$founder_indices), 4L)
    expect_length(unique(x$observed$sample_ids), 4L)
    expect_identical(x$observed$sample_ids, sprintf("S%03d", 1:4))
    expect_length(unique(x$truth$sample_indices), 4L)
    expect_identical(x$truth$ancestralSize, 12L)
    expect_identical(x$truth$reproductiveSizes, c(4L, 4L, 8L))
    expect_equal(x$observed$genotypes,
        x$observed$haplotypes[c(1, 3, 5, 7), , drop = FALSE] +
        x$observed$haplotypes[c(2, 4, 6, 8), , drop = FALSE])
    expect_true(all(x$observed$haplotypes %in% 0:1))
    expect_equal(x$calibration$recombination_fractions,
        c(rep((1 - exp(-.2)) / 2, 2), 0, (1 - exp(-.2)) / 2, 0))
    expect_equal(x$records$generation$heterozygosity,
        rowMeans(2 * x$records$allele_frequency * (1 - x$records$allele_frequency)))
    expect_equal(nrow(x$records$chromosome), 6L)
    expect_false(any(c("sizes", "burnin", "founder_indices", "allele_frequency") %in%
        names(x$observed)))
})

test_that("native mutation acts before reproduction with explicit probabilities", {
    zero <- simulatePopulation(c(10, 10), c(2, 2), 0, 5, seed = 11)
    expect_true(all(zero$observed$haplotypes == 0))
    expect_true(all(zero$records$generation$heterozygosity == 0))
    one <- simulatePopulation(c(10, 10), c(2, 2), 1, 5,
        reverseMutationRate = 0, seed = 11)
    expect_true(all(one$observed$haplotypes == 1))
    expect_equal(one$records$allele_frequency[1, ], rep(0, 4))
    expect_equal(one$records$allele_frequency[2, ], rep(1, 4))
    coarse <- simulatePopulation(rep(10, 6), 1, 0, 2, recordEvery = 2)
    expect_equal(coarse$records$generation$generation, c(0, 2, 4, 5))
})

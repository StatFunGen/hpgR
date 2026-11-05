test_that("R likelihood kernels reproduce original C++ numerical results", {
    fixtures <- readRDS(test_path("fixtures", "tsdate-reference.rds"))
    for (x in fixtures$kernels) {
        actual <- hpgR:::.td_tsdate_precalc_mutation_liks(x$muts, x$spans,
            x$fixed, x$dt, x$mu, TRUE, FALSE)
        expect_equal(actual, x$expected, tolerance = 1e-12)
    }
})

test_that("complete dating matches the original R and C++ calculation", {
    fixtures <- readRDS(test_path("fixtures", "tsdate-reference.rds"))
    for (x in fixtures$pipelines) {
        actual <- dateAncestry(x$ancestry, x$mu, x$Ne)
        expect_equal(actual$dating$mean, x$mean, tolerance = 1e-12)
        expect_equal(actual$dating$variance, x$variance, tolerance = 1e-12)
    }
})

test_that("truth comparisons respect sample labels and genomic length", {
    x <- inferAncestry(rbind(c(1L,1L), c(1L,0L),c(0L,0L)),c(10,20),100)
    result <- compareAncestry(x, x)$lengthWeighted
    expect_equal(unname(result), c(1,1,0,0))
    expect_error(compareAncestry(x, simulateAncestry(2L,100,100)), "same ordered")
})

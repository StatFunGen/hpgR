test_that("shared sequence kernels match the original C fixtures", {
    path <- system.file("reference", "sequence_statistics.tsv", package = "hpgR")
    reference <- read.delim(path)
    actual <- vapply(seq_len(nrow(reference)), function(i) {
        r <- reference[i, ]
        switch(r$kind,
            diversity = hpgR:::tskitDiversity(r$a, r$n),
            D = hpgR:::tskitD(r$a/r$n, (r$a+r$b)/r$n, (r$a+r$c)/r$n),
            r2 = hpgR:::tskitR2(r$a/r$n, (r$a+r$b)/r$n, (r$a+r$c)/r$n),
            polarised = hpgR:::tskitSfs(r$a, r$n)[r$b],
            folded = hpgR:::tskitSfs(r$a, r$n, FALSE)[r$b])
    }, numeric(1))
    expect_identical(is.na(actual), is.na(reference$value))
    expect_equal(actual, reference$value, tolerance = 2e-15)
})

test_that("sequence summaries share their source kernels", {
    for (f in list(summarizePopulationSample, estimateDiversitySize, pgObservedSummaries))
        expect_true("tskitDiversity" %in% all.names(body(f)))
    for (f in list(summarizePopulationSample, pgObservedSummaries, pgMutateTree, sfs))
        expect_true("tskitSfs" %in% all.names(body(f)))
    expect_true("tskitR2" %in% all.names(body(estimateLdSize)))
    expect_true("tskitD" %in% all.names(body(ld)))
    expect_true("simupopParentMaps" %in% all.names(body(simulateForward)))
    expect_true("simupopParentMaps" %in% all.names(body(pgForwardExposure)))
    expect_identical(sfs(matrix(integer(), 4L, 0L)), integer(3L))
    expect_identical(hpgR:::tskitSfs(c(0,4,2,2,1), 4, FALSE), c(1L,2L))
})

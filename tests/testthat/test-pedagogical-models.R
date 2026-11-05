test_that("sampling recipes retain the recorded parental draws", {
    x <- samplingExample()
    expect_identical(x$draws, c(3L,10L,9L,4L,7L,8L,9L,1L,3L,1L))
    expect_identical(x$paths[[1]]$counts[2], 4)
    expect_equal(length(x$paths), 2L)
})

test_that("copying and mutation preserve the saved process invariants", {
    x <- populationProcessExample()
    expect_equal(tail(x$bottleneck$allele_counts, 1), 3)
    expect_identical(tail(x$mutation$heterozygosity, 1), .5)
    expect_identical(x$ancestry$meeting_generation, 2L)
    for (history in x[c("bottleneck", "mutation", "ancestry")]) {
        expect_identical(lengths(history$copies), 2L * history$sizes)
        expect_true(all(history$heterozygosity >= 0 & history$heterozygosity <= 1))
        seen <- unique(history$copies[[1L]])
        for (i in seq_along(history$parents)) {
            expected <- history$copies[[i]][history$parents[[i]] + 1L]
            stable <- setdiff(seq_along(expected), history$mutations[[i]] + 1L)
            expect_identical(expected[stable], history$copies[[i + 1L]][stable])
            fresh <- history$copies[[i + 1L]][history$mutations[[i]] + 1L]
            expect_false(any(fresh %in% seen))
            expect_false(anyDuplicated(fresh) > 0L)
            seen <- union(seen, fresh)
        }
    }
    expect_true(all(diff(x$bottleneck$allele_counts) <= 0))
    expect_true(all(diff(x$ancestry$allele_counts) <= 0))
    expect_true(all(x$bottleneck$allele_counts[8:12] <= 6L))
    for (s in x$recombination) {
        expect_identical(sum(s$counts), s$m)
        expect_true(s$shared[1L] <= .5 && s$shared[2L] >= .5)
    }
    spans <- vapply(x$recombination, function(s) diff(s$shared), 0.0)
    expect_true(all(diff(spans) <= 0))
})

test_that("generation records feed the same conditional variance inverse", {
    x <- simulateFrequencyRecords()
    expect_identical(dim(x$counts), c(20L, 105L))
    expect_true(all(x$counts[1, ] == 16L))
    estimate <- estimate_from_generations(x$records$frequency[1, ],
                                          x$records$frequency[2, ])
    expect_equal(unname(estimate['Ne_V']), 16.07360, tolerance = 1e-5)
    law <- frequencyCountDistribution()
    expect_equal(rowSums(law$transition), rep(1, 33), tolerance = 1e-12)
    expect_equal(rowSums(law$probability), rep(1, 20), tolerance = 1e-12)
})

test_that("small source-model controls retain simple limiting cases", {
    x <- driftMoments(c(16,100), 0:600)
    expect_equal(2 * x$variance + x$diversity,
                 matrix(.5, nrow(x$variance), ncol(x$variance)), tolerance = 1e-12)
    expect_identical(LD(rep(.25,4)), 0)
    expect_identical(LD(c(3,1,1,3)/8), 1/8)
    expect_length(pti_rooted_trees(1:4), 15L)
    G <- rbind(c(0,0,0), c(1,0,0), c(0,1,0), c(0,0,1))
    expect_length(pti_compatible_trees(G), 15L)
    f <- simulateNewMutationPaths()
    expect_identical(dim(f$paths), c(1001L,200L))
    expect_false(anyNA(f$paths[1001, ]))
})

test_that("genealogy animation retains the seeded coalescent realization", {
    x <- coalescentExample()
    expect_identical(x$n, 6L)
    expect_length(x$nodes, 11L)
    expect_length(x$events, 5L)
    expect_length(x$mutations, 9L)
    expect_equal(tail(x$events, 1)[[1]]$time, 2509.8537487382164)
    expect_true(all(vapply(x$mutations, function(m)
        identical(m$carriers, x$nodes[[m$child]]$descendants), TRUE)))
})

test_that("classroom Kingman trees use the shared source-derived event kernel", {
    for (n in c(2L, 3L, 5L, 12L)) {
        source <- hpgR:::.msprimeAncestryCopies(n, 80, 1L, seed = 820001L)
        tree <- pg_kingman_tree(80, n, ancestry_seed = 820001L)
        expect_identical(tree$nodes$id, source$nodes$id + 1L)
        expect_identical(tree$nodes$time, source$nodes$time)
        expect_identical(tree$edges$parent, source$edges$parent + 1L)
        expect_identical(tree$edges$child, source$edges$child + 1L)
    }
    for (samples in c(1L, 2L, 3L))
        expect_identical(simulateAncestry(samples, 80, 10, seed = 820002L),
            hpgR:::.msprimeAncestryCopies(2L * samples, 80, 10, seed = 820002L))
    expect_true(".msprimeAncestryCopies" %in% all.names(body(pg_kingman_tree)))
})

test_that("classroom tree mutations use the shared source-derived mutation kernel", {
    tree <- pg_kingman_tree(80, 5L, ancestry_seed = 820003L)
    for (mu in c(0, .0002)) {
        ancestry <- hpgR:::.pgTreeAncestry(tree, 1000)
        source <- mutateAncestry(ancestry, mu, seed = 820004L)
        sample <- pg_mutate_tree(tree, 1000, mu, mutation_seed = 820004L)
        expect_identical(unname(sample$haplotypes), ancestryData(source)$observed$haplotypes)
        expect_identical(sample$positions_bp, source$sites$position)
        expect_identical(sample$mutations$time, source$mutations$time)
        expect_identical(sample$mutations$child, tree$nodes$id[source$mutations$node + 1L])
        expect_equal(sample$ancestral_exposure, 1000 * sum(tree$edges$length))
    }
    expect_true(".msprimeMutations" %in% all.names(body(pg_mutate_tree)))
})

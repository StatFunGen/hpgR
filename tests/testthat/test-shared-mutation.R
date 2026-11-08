test_that("recorded animation events replay through the source kernels", {
    x <- coalescentExample()
    events <- lapply(x$events, function(e)
        list(type = "common_ancestor", time = e$time, lineages = e$pair - 1L))
    ancestry <- hpgR:::msprimeAncestryCopies(x$n, x$Ne, x$L,
        seed = x$seed, events = events)
    expect_identical(ancestry$nodes$time, vapply(x$nodes, `[[`, 0, "time"))
    expect_identical(ancestry$edges$parent + 1L,
                     vapply(x$edges, `[[`, 0L, "parent"))
    expect_identical(ancestry$edges$child + 1L,
                     vapply(x$edges, `[[`, 0L, "child"))
    mutations <- data.frame(edge = match(vapply(x$mutations, `[[`, 0L, "child") - 1L,
                                         ancestry$edges$child),
        position = vapply(x$mutations, `[[`, 0, "position"),
        time = vapply(x$mutations, `[[`, 0, "time"))
    mutated <- hpgR:::msprimeMutations(ancestry, x$mu, x$seed, events = mutations)
    observed <- ancestryData(mutated)$observed
    G <- matrix(0L, x$n, length(x$mutations))
    for (i in seq_along(x$mutations)) G[x$mutations[[i]]$carriers, i] <- 1L
    expect_identical(observed$haplotypes, G)
    expect_identical(observed$positions, mutations$position)
    expect_identical(mutated$mutations$time, mutations$time)
})

test_that("finite-copy mutations use the common placement kernel", {
    for (map in c(0, .04)) {
        record <- pgRecordAncestry(rep(4L, 6L), 100, 4L, map,
                                     ancestry_seed = 812312L)
        edges <- do.call(rbind, lapply(pgLocalTrees(record), `[[`, "edges"))
        ancestry <- structure(list(nodes = data.frame(id = record$nodes$id - 1L,
            time = record$nodes$time, flags = as.integer(record$nodes$id %in% record$sample_nodes)),
            edges = data.frame(left = edges$left, right = edges$right,
                parent = edges$parent - 1L, child = edges$child - 1L),
            sampleIds = record$sample_nodes - 1L, sequenceLength = record$length_bp,
            timeUnits = "generations"), class = c("hpgR_ancestry", "list"))
        for (mu in c(0, .002)) {
            source <- hpgR:::msprimeMutations(ancestry, mu, 812313L, batchEdges = TRUE)
            x <- pgMutateAncestry(record, mu, 812313L)
            expect_identical(x$positions_bp, source$sites$position)
            expect_identical(x$mutations$time, source$mutations$time)
            expect_identical(x$haplotypes, ancestryData(source)$observed$haplotypes)
            expect_identical(x$haplotypes, pgReplay(record, x$mutations))
            expect_identical(x$ancestral_exposure, sum((edges$right - edges$left) *
                (record$nodes$time[edges$parent] - record$nodes$time[edges$child])))
        }
    }
    expect_true("msprimeMutations" %in% all.names(body(pgMutateAncestry)))
    expect_false(any(c("rpois", "runif") %in% all.names(body(pgMutateAncestry))))
    for (seed in list(-1L, 1.5, NA_real_, Inf, c(1L, 2L)))
        expect_error(pgMutateAncestry(record, .001, seed), "seed")
})


test_that("discrete mutation uses integer sites inside fractional edge spans", {
    ancestry <- structure(list(
        nodes = data.frame(id = 0:2, time = c(0, 0, 10), flags = c(1L, 1L, 0L)),
        edges = data.frame(left = .2, right = 2.2, parent = 2L, child = 0:1),
        sampleIds = 0:1, sequenceLength = 3, timeUnits = "generations"),
        class = c("hpgR_ancestry", "list"))
    set.seed(9181L)
    before <- .Random.seed
    mutated <- mutateAncestry(ancestry, 1, seed = 9182L, model = "binary")
    expect_identical(.Random.seed, before)
    expect_equal(mutated$sites$position, c(1, 2))
    expect_true(all(mutated$sites$position >= .2 & mutated$sites$position < 2.2))

    # [0.2, 0.8) contains no integer position and has zero mutation exposure.
    ancestry$edges$right <- .8
    empty <- mutateAncestry(ancestry, 100, seed = 9182L, model = "binary")
    expect_equal(nrow(empty$sites), 0L)
    expect_equal(nrow(empty$mutations), 0L)
})

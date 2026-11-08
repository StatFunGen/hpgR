test_that("site recording preserves DNA and follows actual transmitted copies", {
    args <- list(sizes = c(5L, 2L, 4L), loci = c(4L, 4L), mutationRate = .2,
        reverseMutationRate = .1, sampleSize = 4L, burnin = 3L, founderSize = 2L,
        geneticPositions = list(c(0, .4, .8, 1.2), c(0, .6, .9, 1.2)), seed = 217L)
    ordinary <- do.call(simulatePopulation, args)
    native <- hpgR:::simupopMeiosis
    steps <- list()
    local_mocked_bindings(simupopMeiosis = function(...) {
        inputs <- list(...)
        result <- native(...)
        steps[[length(steps) + 1L]] <<- list(parental = inputs[[1L]],
            offspring = if (is.list(result)) result$haplotypes else result)
        result
    }, .package = "hpgR")
    sites <- c(1L, 2L, 4L, 5L, 8L)
    set.seed(307L)
    rng <- .Random.seed
    recorded <- do.call(simulatePopulation, c(args, list(recordSites = sites)))
    expect_identical(.Random.seed, rng)
    expect_identical(recorded$observed, ordinary$observed)
    expect_identical(recorded$calibration, ordinary$calibration)
    expect_identical(recorded$truth, ordinary$truth)
    expect_identical(recorded$records[names(ordinary$records)], ordinary$records)
    expect_length(recorded$records$copyParents$parents, length(steps))
    for (g in seq_along(steps)) for (j in seq_along(sites)) {
        parent <- recorded$records$copyParents$parents[[g]][, j]
        # Read the mutated parental matrix directly, independently of crossover parsing.
        expect_identical(steps[[g]]$parental[parent, sites[j]],
                         steps[[g]]$offspring[, sites[j]])
    }
    genealogy <- populationGenealogy(recorded, site = 4L, samples = c(1, 4, 5, 8))
    expect_silent(pgValidateAncestry(genealogy, require_complete = FALSE))
    expect_equal(genealogy$sample_generation, 2)
    expect_equal(min(genealogy$nodes$generation), -3)
    expect_identical(genealogy$sample_ids, c("H1", "H4", "H5", "H8"))
    initial <- genealogy$nodes$id[genealogy$nodes$generation == 0]
    first <- genealogy$edges[genealogy$nodes$generation[genealogy$edges$child] == 1, ]
    founderRows <- as.vector(rbind(2L * recorded$truth$founder_indices - 1L,
                                  2L * recorded$truth$founder_indices))
    localParents <- recorded$records$copyParents$parents[[4L]][, 3L]
    expect_identical(first$parent, initial[founderRows[localParents]])
    expect_error(populationGenealogy(ordinary), "recordSites")
    expect_error(populationGenealogy(recorded, site = 3L), "not recorded")
})

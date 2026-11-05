test_that("complete-example source is registered separately from reusable models", {
    manifest <- utils::read.delim(system.file("provenance", "modules.tsv", package = "hpgR"),
                                 stringsAsFactors = FALSE)
    expect_true(all(c("simupop.R", "pedagogical_effective_size.R",
                      "pedagogical_demographic_history_example.R") %in% manifest$file))
    entry <- manifest[manifest$file == "pedagogical_demographic_history_example.R", ]
    expect_identical(entry$kind, "pedagogical")
    expect_match(entry$validation_scope, "complete classroom")
})

test_that("the example schedule preserves the declared founder and growth history", {
    candidate <- list(founder = 20L, small = 20L, duration = 20L, growth = 40L)
    design <- list(ancestral_size = 100L, present_size = 1000L)
    actual <- populationSizeSchedule(candidate, design)
    expected <- as.integer(c(100L, rep(20L, 20L),
        round(exp(seq(log(20), log(1000), length.out = 41)))[-1L]))
    expect_identical(actual, expected)
    expect_identical(tail(actual, 1L), 1000L)
    expect_length(actual, 61L)
})

test_that("complete-example options are explicit before any simulation", {
    expect_error(populationHistoryExample(NULL, NULL, NULL), "predictionData")
    expect_error(populationHistoryExample(NULL, NULL, list(), simulate = NA), "logical")
    expect_error(populationHistoryExample(NULL, NULL, list(), simulate = 1), "logical")
    expect_error(populationHistoryExample(NULL, NULL, list(),
        regeneratePredictions = c(TRUE, FALSE)), "logical")
    invalid <- list(design = list(), candidates = list(a = list()),
                    labels = "a", banks = list(b = list()))
    expect_error(populationHistoryExample(NULL, NULL, invalid), "same named histories")
})

test_that("the population-history recipe fixes independent R replicate seeds", {
    x <- populationHistoryDesign()
    expect_identical(x$design$observation_seed, 20261103L)
    expect_length(x$candidates, 4L)
    seeds <- unlist(lapply(x$banks, function(b)
        c(b$training_seeds, b$calibration_seeds)))
    expect_length(seeds, 128L)
    expect_false(anyDuplicated(seeds) > 0L)
    expect_false(x$design$observation_seed %in% seeds)
})

test_that("serial and parallel R prediction rebuilds use identical seeded draws", {
    x <- populationHistoryDesign()
    x$design$ancestral_size <- 16L
    x$design$present_size <- 16L
    x$design$chromosomes <- 2L
    x$design$sites_per_chromosome <- 100L
    x$design$mutation_rate <- .02
    x$design$burnin_generations <- 30L
    x$design$sample_diploids <- 8L
    x$design$fitting_chromosomes <- 1L
    x$design$checking_chromosomes <- 2L
    x$candidates <- list(a = list(founder = 6L, small = 6L, duration = 2L, growth = 2L),
        b = list(founder = 10L, small = 10L, duration = 2L, growth = 2L))
    x$labels <- c(a = "Six founders", b = "Ten founders")
    x$banks <- list(a = list(training_seeds = 38001:38002, calibration_seeds = 38101:38102),
        b = list(training_seeds = 38201:38202, calibration_seeds = 38301:38302))
    serial <- populationHistoryPredictions(x)
    expect_identical(serial$engine$language, "R")
    expect_identical(dim(serial$banks$a$train_fit), c(2L, 6L))
    expect_true(all(is.finite(serial$banks$a$train_fit)))
    if (.Platform$OS.type != "windows") {
        parallel <- populationHistoryPredictions(x, workers = 2L)
        expect_identical(parallel, serial)
    }
})

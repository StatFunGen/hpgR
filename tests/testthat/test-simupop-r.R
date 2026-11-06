test_that("forward copying is self-contained and replays its recorded parents", {
    initial <- matrix(seq_len(24L), nrow = 6L)
    x <- simulateForward(initial, c(3, 2, 5, 3), seed = 143L)
    expect_identical(x, simulateForward(initial, c(3, 2, 5, 3), seed = 143L))
    expect_identical(x$states, replayCopying(initial, x$parents))
    expect_identical(lapply(x$states, dim),
        list(c(6L, 4L), c(4L, 4L), c(10L, 4L), c(6L, 4L)))
    expect_true(all(vapply(seq_along(x$parents), function(i)
        all(x$parents[[i]] >= 1L & x$parents[[i]] <= nrow(x$states[[i]])),
        logical(1))))
})

test_that("forward simulation restores the caller's random generator", {
    set.seed(725L)
    before <- .Random.seed
    kind <- RNGkind()
    invisible(simulateForward(matrix(0:1, ncol = 1L), c(1, 3, 2), 91L))
    expect_identical(.Random.seed, before)
    expect_identical(RNGkind(), kind)
    invisible(simulatePopulation(c(4, 3), 3, .05, 2,
        burnin = 4L, seed = 91L))
    expect_identical(.Random.seed, before)
    expect_identical(RNGkind(), kind)
})

test_that("seeded draws restore the caller's generator after an error", {
    restore <- hpgR:::rngRestorer()
    on.exit(restore(), add = TRUE)
    RNGkind("L'Ecuyer-CMRG", "Box-Muller", "Rejection")
    set.seed(725L)
    before <- .Random.seed
    kind <- RNGkind()
    expect_error(hpgR:::simupopWithSeed(91L, function() {
        runif(3L)
        stop("draw failed")
    }), "draw failed")
    expect_identical(.Random.seed, before)
    expect_identical(RNGkind(), kind)
    result <- withVisible(hpgR:::simupopWithSeed(91L, function() invisible(2L)))
    expect_identical(result, list(value = 2L, visible = FALSE))
    expect_identical(.Random.seed, before)
})

test_that("a simulation does not leave an RNG state when the caller had none", {
    # Observe entry and exit before returning to the test reporter.
    inspectAbsentSeed <- function(draw) {
        old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            get(".Random.seed", envir = .GlobalEnv) else NULL
        on.exit(if (is.null(old)) {
            if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
                rm(".Random.seed", envir = .GlobalEnv)
        } else assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
        before <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        invisible(draw())
        after <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        c(before = before, after = after)
    }
    forward <- inspectAbsentSeed(function()
        simulateForward(matrix(0:1, ncol = 1L), c(1, 2), 91L))
    population <- inspectAbsentSeed(function()
        simulatePopulation(c(4, 3), 2, .1, 2, seed = 91L))
    expect_identical(forward, c(before = FALSE, after = FALSE))
    expect_identical(population, c(before = FALSE, after = FALSE))
})

test_that("the public population record has exact dimensions and isolated truth", {
    args <- list(sizes = c(12, 4, 8, 6), loci = c(3, 2), mutationRate = .03,
        sampleSize = 4, geneticPositions = list(c(0, .1, .2), c(0, .1)),
        burnin = 7L, founderSize = 5L, seed = 8149L)
    x <- do.call(simulatePopulation, args)
    expect_identical(x, do.call(simulatePopulation, args))
    args$recordEvery <- 2L
    coarse <- do.call(simulatePopulation, args)
    expect_identical(x$observed, coarse$observed)
    expect_identical(dim(x$observed$haplotypes), c(8L, 5L))
    expect_identical(dim(x$observed$genotypes), c(4L, 5L))
    expect_identical(dim(x$records$allele_frequency), c(4L, 5L))
    expect_equal(x$records$generation$generation, 0:3)
    expect_equal(coarse$records$generation$generation, c(0, 2, 3))
    expect_equal(x$records$generation$size, c(5, 4, 8, 6))
    expect_equal(x$records$burnin$generation, -7:0)
    expect_identical(x$truth$ancestralSize, 12L)
    expect_identical(x$truth$reproductiveSizes, c(5L, 4L, 8L, 6L))
    expect_length(unique(x$truth$founder_indices), 5L)
    expect_true(all(x$truth$founder_indices %in% 1:12))
    expect_length(unique(x$truth$sample_indices), 4L)
    expect_true(all(x$truth$sample_indices %in% 1:6))
    expect_identical(x$observed$sample_ids, sprintf("S%03d", 1:4))
    expect_equal(x$observed$genotypes,
        x$observed$haplotypes[c(1, 3, 5, 7), , drop = FALSE] +
        x$observed$haplotypes[c(2, 4, 6, 8), , drop = FALSE])
    expect_true(all(x$observed$haplotypes %in% 0:1))
    expect_equal(x$records$generation$heterozygosity,
        rowMeans(2 * x$records$allele_frequency * (1 - x$records$allele_frequency)))
    expect_equal(x$calibration$recombination_fractions,
        c(rep(-expm1(-.2) / 2, 2), 0, -expm1(-.2) / 2, 0))
    expect_equal(nrow(x$records$chromosome), 8L)
    expect_false(any(c("sizes", "burnin", "founder_indices", "allele_frequency",
        "sample_indices", "reproductiveSizes") %in% names(x$observed)))
})

test_that("mutation occurs on the parental pool before a shared mating draw", {
    # Starting from all zeros, certain forward mutation must be visible in the
    # first offspring. Certain reverse mutation must not undo that same event.
    x <- simulatePopulation(rep(6, 5), c(3, 2), 1, 4,
        reverseMutationRate = 1, seed = 391L)
    expect_equal(x$records$allele_frequency,
        matrix(rep(c(0, 1, 0, 1, 0), 5L), ncol = 5L))
    expect_true(all(x$observed$haplotypes == 0L))
    y <- simulatePopulation(c(6, 6), 7, 1, 4,
        reverseMutationRate = 0, seed = 391L)
    expect_true(all(y$observed$haplotypes == 1L))
    z <- simulatePopulation(c(6, 6), 7, 0, 4, burnin = 10L, seed = 391L)
    expect_true(all(z$observed$haplotypes == 0L))
    expect_true(all(z$records$generation$heterozygosity == 0))
})

test_that("zero-generation sampling and founder selection do not copy an individual twice", {
    x <- simulatePopulation(12L, 4L, 0, 5L, founderSize = 7L, seed = 392L)
    expect_identical(dim(x$observed$haplotypes), c(10L, 4L))
    expect_equal(x$records$generation$size, 7L)
    expect_equal(x$records$generation$generation, 0L)
    expect_length(unique(x$truth$sample_indices), 5L)
    expect_true(all(x$truth$sample_indices %in% 1:7))
    expect_length(unique(x$truth$founder_indices), 7L)
    expect_true(all(x$truth$founder_indices %in% 1:12))
    expect_identical(x$truth$reproductiveSizes, 7L)
})

test_that("binary mutation replays recorded events simultaneously", {
    mutate <- getFromNamespace("simupopMutate", "hpgR")
    h <- matrix(rep(0:1, 6L), nrow = 4L)
    expected <- h
    changed <- c(1L, 4L, 7L, 12L)
    expected[changed] <- 1L - expected[changed]
    expect_identical(mutate(h, .2, .3, events = changed), expected)
    expect_identical(mutate(h, 1, 1), 1L - h)
    expect_identical(mutate(h, 0, 0), h)
    set.seed(395L)
    recorded <- mutate(h, .2, .3, returnEvents = TRUE)
    expect_identical(mutate(h, .2, .3, events = recorded$events),
        recorded$haplotypes)
})

test_that("specified gamete starts and crossover switches give exact mosaics", {
    meiosis <- getFromNamespace("simupopMeiosis", "hpgR")
    h <- rbind(0:5, 10:15, 20:25, 30:35)
    events <- list(starts = matrix(c(0L, 1L, 1L, 0L, 0L, 0L, 1L, 1L),
            nrow = 4L, byrow = TRUE),
        switches = list(matrix(c(1L, 1L, 2L, 2L), ncol = 2L, byrow = TRUE),
            matrix(c(3L, 1L, 4L, 2L), ncol = 2L, byrow = TRUE)))
    expected <- rbind(c(0L, 11L, 12L, 13L, 14L, 15L),
        c(10L, 11L, 2L, 3L, 4L, 5L),
        c(20L, 21L, 22L, 23L, 34L, 35L),
        c(30L, 31L, 32L, 33L, 34L, 25L))
    expect_identical(meiosis(h, c(1L, 1L, 2L, 2L), c(3L, 3L),
        list(c(.2, .2, 0), c(.2, .2, 0)), events = events), expected)
    set.seed(396L)
    recorded <- meiosis(h, c(1L, 1L, 2L, 2L), c(3L, 3L),
        list(c(.2, .2, 0), c(.2, .2, 0)), returnEvents = TRUE)
    expect_identical(meiosis(h, c(1L, 1L, 2L, 2L), c(3L, 3L),
        list(c(.2, .2, 0), c(.2, .2, 0)), events = recorded$events),
        recorded$haplotypes)
})

test_that("gametes choose parental homologues independently including selfing", {
    meiosis <- getFromNamespace("simupopMeiosis", "hpgR")
    choose <- getFromNamespace("simupopParents", "hpgR")
    set.seed(393L)
    p <- choose(1L, 10000L)
    expect_identical(p, rep(1L, 20000L))
    h <- meiosis(matrix(0:1, ncol = 1L), p, 1L, list(0))
    dosage <- h[seq(1L, nrow(h), 2L), 1L] + h[seq(2L, nrow(h), 2L), 1L]
    expect_lt(max(abs(as.numeric(tabulate(dosage + 1L, nbins = 3L) / 10000) -
        c(.25, .5, .25))), .018)
    set.seed(398L)
    multiple <- choose(7L, 10000L)
    expect_lt(max(abs(as.numeric(tabulate(multiple, nbins = 7L) /
        length(multiple)) - 1 / 7)), .01)
    expect_lt(abs(mean(multiple[seq(1L, length(multiple), 2L)] ==
        multiple[seq(2L, length(multiple), 2L)]) - 1 / 7), .012)
})

test_that("chromosomes assort independently and adjacent switches use supplied rates", {
    meiosis <- getFromNamespace("simupopMeiosis", "hpgR")
    h <- rbind(rep(0L, 4L), rep(1L, 4L))
    parents <- rep(1L, 20000L)
    set.seed(394L)
    y <- meiosis(h, parents, c(2L, 2L), list(c(.2, 0), c(.35, 0)))
    expect_lt(abs(mean(y[, 1L] != y[, 2L]) - .2), .015)
    expect_lt(abs(mean(y[, 3L] != y[, 4L]) - .35), .015)
    expect_lt(abs(mean(y[, 1L] != y[, 3L]) - .5), .015)
    set.seed(394L)
    zero <- meiosis(h, parents, c(2L, 2L), list(c(0, 0), c(0, 0)))
    expect_identical(zero[, 1L], zero[, 2L])
    expect_identical(zero[, 3L], zero[, 4L])
})

test_that("supplied parent and individual-sampling records replay exactly", {
    choose <- getFromNamespace("simupopParents", "hpgR")
    sample_individuals <- getFromNamespace("simupopSample", "hpgR")
    expect_identical(choose(4L, 2L, parents = c(1L, 1L, 3L, 4L)),
        c(1L, 1L, 3L, 4L))
    h <- matrix(seq_len(32L), nrow = 8L)
    x <- sample_individuals(h, 2L, indices = c(4L, 2L))
    expect_identical(x$indices, c(2L, 4L))
    expect_identical(x$haplotypes, unname(h[c(3L, 4L, 7L, 8L), , drop = FALSE]))
    set.seed(397L)
    recorded <- sample_individuals(h, 3L)
    expect_length(unique(recorded$indices), 3L)
    expect_identical(sample_individuals(h, 3L, indices = recorded$indices), recorded)
})

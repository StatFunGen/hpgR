# Population-history recipes compose the shared simulation and inference functions.

#' Population sizes for a founder, bottleneck and expansion history
#' @param candidate Founder size, small-population size, duration and growth.
#' @param design Ancestral and present size and sequence/sampling design.
#' @return Integer diploid sizes, starting with the ancestral population.
#' @export
populationSizeSchedule <- function(candidate, design) {
    growth <- round(exp(seq(log(candidate$small), log(design$present_size),
                           length.out = candidate$growth + 1L)))[-1L]
    as.integer(c(design$ancestral_size,
                 rep(candidate$small, candidate$duration), growth))
}

#' Simulate the declared population history
#' @param candidate Founder size, small-population size, duration and growth.
#' @param design Sequence, rates and sampling specification.
#' @param seed R random seed.
#' @param recordEvery Interval between stored generation records.
#' @return The observed DNA, calibration and separate simulation records.
#' @export
simulatePopulationHistory <- function(candidate, design, seed,
                                         recordEvery = 10000L) {
    simulatePopulation(sizes = populationSizeSchedule(candidate, design),
        loci = rep(design$sites_per_chromosome, design$chromosomes),
        mutationRate = design$mutation_rate,
        reverseMutationRate = design$mutation_rate,
        sampleSize = design$sample_diploids,
        geneticPositions = rep(list(seq(0, design$genetic_length_morgan,
            length.out = design$sites_per_chromosome)), design$chromosomes),
        burnin = design$burnin_generations,
        founderSize = candidate$founder, recordEvery = recordEvery, seed = seed)
}

#' Sequence summaries for population-history comparison
#' @param observed Final-sample haplotypes and marker coordinates.
#' @param design LD distance bins and other declared design choices.
#' @param chromosomes Chromosomes used for this calculation.
#' @return Grouped SFS proportions, distance-binned LD and long IBD sharing.
#' @export
populationSampleFeatures <- function(observed, design, chromosomes) {
    markers <- observed$markers
    take <- markers$chromosome %in% chromosomes
    H <- observed$haplotypes[, take, drop = FALSE]
    groups <- if (!is.null(design$sfs_group_breaks)) {
        edges <- design$sfs_group_breaks
        counts <- seq_len(floor(nrow(H) / 2L))
        bins <- cut(counts, edges, labels = FALSE, right = TRUE)
        stats::setNames(lapply(seq_len(length(edges) - 1L), function(b)
            counts[bins == b]), paste0("bin", seq_len(length(edges) - 1L)))
    } else NULL
    dna <- summarizePopulationSample(H, chromosome = markers$chromosome[take],
                                      sfsGroups = groups)
    ld <- estimateLdSize(H[seq(1L, nrow(H), 2L), , drop = FALSE],
        geneticPositions = markers$genetic_position[take],
        chromosome = markers$chromosome[take],
        distanceBreaks = if (is.null(design$history_ld_distance_breaks))
            design$ld_distance_breaks else design$history_ld_distance_breaks,
        maxPairs = 500000L, seed = 101L)
    ibd <- if (!is.null(design$ibd_length_breaks))
        populationIbdSummary(observed, design, chromosomes)$bins$fraction else numeric()
    c(stats::setNames(dna$groupedSfsPerSite,
                     paste0("sfs_", names(dna$groupedSfsPerSite))),
      stats::setNames(ld$bins$meanRho2, paste0("ld_", seq_len(nrow(ld$bins)))),
      stats::setNames(ibd, paste0("ibd_", seq_along(ibd))))
}

# The classroom DNA has known phase. One sampled haplotype per person is used
# for both LD and long-segment matching; invariant and rare sites are excluded
# from matching, while diversity and the SFS use every analyzed base.
populationIbdSummary <- function(observed, design, chromosomes) {
    markers <- observed$markers
    take <- markers$chromosome %in% chromosomes
    h <- observed$haplotypes[seq.int(1L, nrow(observed$haplotypes), 2L), take, drop = FALSE]
    m <- markers[take, , drop = FALSE]
    lengths <- vapply(split(m$genetic_position, m$chromosome),
        function(p) 100 * diff(range(p)), numeric(1))
    p <- colMeans(h)
    use <- p >= design$ibd_minimum_maf & p <= 1 - design$ibd_minimum_maf
    segments <- if (any(use)) detectIbdSegments(h[, use, drop = FALSE],
        m$genetic_position[use], m$chromosome[use], wordSize = 1L,
        minimumCm = design$ibd_length_breaks[1L]) else
        data.frame(lengthCm = numeric())
    pairs <- choose(nrow(h), 2L)
    list(segments = segments, haplotypePairs = pairs, chromosomeLengthsCm = lengths,
        markersUsed = sum(use), bins = summarizeIbdSharing(segments, pairs, lengths,
                                                          design$ibd_length_breaks))
}

examplePopulationPredictions <- function(candidates, design, banks, workers = 1L) {
    checkInputs(workers = scalarInteger)
    workers <- as.integer(workers)
    jobs <- unlist(lapply(names(candidates), function(id) {
        seeds <- c(banks[[id]]$training_seeds, banks[[id]]$calibration_seeds)
        lapply(seq_along(seeds), function(j) list(id = id, index = j, seed = seeds[j]))
    }), recursive = FALSE)
    generate <- function(job) {
        draw <- simulatePopulationHistory(candidates[[job$id]], design, job$seed)
        list(fit = populationSampleFeatures(draw$observed, design,
                                              design$fitting_chromosomes),
             check = populationSampleFeatures(draw$observed, design,
                                                design$checking_chromosomes))
    }
    draws <- if (workers == 1L) lapply(jobs, generate) else {
        if (.Platform$OS.type == "windows")
            stop("Use workers = 1 on Windows", call. = FALSE)
        parallel::mclapply(jobs, generate, mc.cores = min(workers, length(jobs)),
                          mc.set.seed = FALSE)
    }
    if (any(vapply(draws, inherits, logical(1), "try-error")))
        stop("A population prediction failed; cache was not completed", call. = FALSE)
    result <- lapply(names(candidates), function(id) {
        bank <- banks[[id]]
        seeds <- c(bank$training_seeds, bank$calibration_seeds)
        ntrain <- length(bank$training_seeds)
        take <- which(vapply(jobs, function(job) job$id == id, logical(1)))
        fit <- lapply(draws[take], `[[`, "fit")
        check <- lapply(draws[take], `[[`, "check")
        train <- seq_len(ntrain)
        calibration <- seq.int(ntrain + 1L, length(seeds))
        list(train_fit = do.call(rbind, fit[train]),
             train_check = do.call(rbind, check[train]),
             cal_fit = do.call(rbind, fit[calibration]),
             cal_check = do.call(rbind, check[calibration]),
             training_seeds = bank$training_seeds,
             calibration_seeds = bank$calibration_seeds)
    })
    names(result) <- names(candidates)
    result
}

#' Define the population-history worked example
#' @return Generating design, four candidate histories, labels, and disjoint
#'   training/calibration seed lists. Population sizes count diploid people;
#'   genetic distances are in Morgans and mutation rates are per site per generation.
#' @export
populationHistoryDesign <- function() {
    design <- list(ancestral_size = 100L, founder_size = 20L,
        bottleneck_generations = 20L, expansion_generations = 40L,
        present_size = 1000L, chromosomes = 8L, sites_per_chromosome = 10000L,
        mutation_rate = 2.5e-5, genetic_length_morgan = .2,
        sample_diploids = 64L, burnin_generations = 1000L,
        observation_seed = 20261103L, ld_distance_breaks = c(.03, .07, .12, .2),
        fitting_chromosomes = 1:6, checking_chromosomes = 7:8,
        sfs_group_breaks = c(0, 1, 2, 4, 8, 16, 32, Inf),
        history_ld_distance_breaks = c(.03, .05, .07, .10, .13, .16, .20),
        ibd_length_breaks = c(8, 10, 14, Inf), ibd_minimum_maf = .1)
    candidates <- list(expansion_only = list(founder = NULL, small = 100L,
        duration = 20L, growth = 40L), shallow_founding = list(founder = 40L,
        small = 40L, duration = 20L, growth = 40L),
        recent_founding = list(founder = 20L, small = 20L, duration = 20L, growth = 40L),
        older_founding = list(founder = 20L, small = 20L, duration = 20L, growth = 80L))
    labels <- stats::setNames(c("Expansion only", "40 founders, 60 generations ago",
        "20 founders, 60 generations ago", "20 founders, 100 generations ago"),
        names(candidates))
    banks <- lapply(seq_along(candidates), function(j)
        list(training_seeds = as.integer(100000L + 1000L * j + 1:40),
             calibration_seeds = as.integer(200000L + 1000L * j + 1:40)))
    names(banks) <- names(candidates)
    list(design = design, candidates = candidates, labels = labels, banks = banks)
}

#' Generate population-history predictions in R
#' @param predictionData Design and seeds from populationHistoryDesign.
#' @param workers Number of R worker processes. Each replicate has its own seed.
#' @return Design, candidate labels, simulated feature banks and the complete
#'   seed ledger. Fitting and checking chromosomes share each
#'   replicate's realized pedigree.
#' @export
populationHistoryPredictions <- function(predictionData = populationHistoryDesign(),
                                         workers = 1L) {
    required <- c("design", "candidates", "labels", "banks")
    if (!is.list(predictionData) || !all(required %in% names(predictionData)))
        stop("predictionData must contain design, candidates, labels and banks",
             call. = FALSE)
    result <- predictionData[required]
    result$banks <- examplePopulationPredictions(result$candidates,
        result$design, result$banks, workers)
    result$seed_ledger <- unlist(lapply(result$banks, function(bank)
        c(bank$training_seeds, bank$calibration_seeds)), use.names = TRUE)
    result
}

#' Run the complete population-history worked example
#'
#' Use the same final DNA for diversity, LD, IBD and candidate-history comparison.
#' The example selects one haplotype per individual for LD and reserves the
#' chromosomes declared in the supplied design for a prediction check.
#' Simulation, estimation and comparison call the reusable hpgR functions.
#'
#' @param observed Final DNA object distributed with the worked example.
#' @param calibration Mutation rate and genetic-map calibration for that DNA.
#' @param predictionData Worked-example object containing design, candidates,
#'   labels and independent R prediction/calibration banks.
#' @param simulate Logical. Reproduce the final DNA and calibration using
#'   the design's observation seed and check their exact agreement.
#' @param regeneratePredictions Logical. Regenerate every candidate bank using
#'   its stored seeds in R.
#' @return A list with simulation records (when requested), DNA summaries,
#'   diversity and LD estimates, fitting/checking features, candidate libraries,
#'   fitted history and held-out check.
#' @details This is the complete classroom recipe, with its fixed LD filtering
#'   and chromosome split. General-purpose calculations are documented under
#'   simulatePopulation, summarizePopulationSample, estimateDiversitySize,
#'   estimateLdSize, buildHistoryLibrary and comparePopulationHistories.
#' @export
populationHistoryExample <- function(observed, calibration, predictionData,
                                     simulate = FALSE,
                                     regeneratePredictions = FALSE) {
    if (!is.logical(simulate) || length(simulate) != 1L || is.na(simulate) ||
        !is.logical(regeneratePredictions) || length(regeneratePredictions) != 1L ||
        is.na(regeneratePredictions)) {
        stop("simulate and regeneratePredictions must be single logical values",
             call. = FALSE)
    }
    required <- c("design", "candidates", "labels", "banks")
    if (!is.list(predictionData) || !all(required %in% names(predictionData))) {
        stop("predictionData must contain design, candidates, labels and banks",
             call. = FALSE)
    }
    design <- predictionData$design
    candidates <- predictionData$candidates
    banks <- predictionData$banks
    if (!is.list(candidates) || !length(candidates) ||
        is.null(names(candidates)) || anyDuplicated(names(candidates)) ||
        !setequal(names(candidates), names(banks))) {
        stop("candidates and banks must identify the same named histories",
             call. = FALSE)
    }
    simulation <- NULL
    if (simulate) {
        if (!"recent_founding" %in% names(candidates)) {
            stop("the worked example requires the recent_founding candidate",
                 call. = FALSE)
        }
        simulation <- simulatePopulationHistory(candidates$recent_founding,
            design, design$observation_seed, recordEvery = 1L)
        if (!identical(simulation$observed, observed) ||
            !identical(simulation$calibration, calibration)) {
            stop("supplied DNA or calibration does not match the generating design and seed",
                 call. = FALSE)
        }
    }
    if (regeneratePredictions) {
        banks <- examplePopulationPredictions(candidates, design, banks)
    }
    dna_summary <- summarizePopulationSample(observed$haplotypes,
        chromosome = observed$markers$chromosome)
    diversity_fit <- estimateDiversitySize(observed$haplotypes,
        mutationRate = calibration$mutationRate,
        sequenceLength = calibration$sequenceLength)
    finite_state_fit <- estimateDiversitySize(observed$haplotypes,
        mutationRate = calibration$mutationRate,
        sequenceLength = calibration$sequenceLength,
        mutationModel = "symmetric_parental")
    ld_fit <- estimateLdSize(
        observed$haplotypes[seq(1L, nrow(observed$haplotypes), 2L), , drop = FALSE],
        geneticPositions = observed$markers$genetic_position,
        chromosome = observed$markers$chromosome,
        distanceBreaks = design$ld_distance_breaks,
        maxPairs = 500000L, seed = 101L)
    ibd_summary <- if (!is.null(design$ibd_length_breaks))
        populationIbdSummary(observed, design, seq_len(design$chromosomes)) else NULL
    ibd_fit <- if (!is.null(ibd_summary)) estimateIbdSize(ibd_summary$segments,
        ibd_summary$haplotypePairs, ibd_summary$chromosomeLengthsCm,
        lengthBreaks = design$ibd_length_breaks) else NULL
    fit_features <- populationSampleFeatures(observed, design,
                                              design$fitting_chromosomes)
    check_features <- populationSampleFeatures(observed, design,
                                                design$checking_chromosomes)
    history_library <- buildHistoryLibrary(lapply(banks, `[[`, "train_fit"),
                                          lapply(banks, `[[`, "cal_fit"))
    heldout_library <- buildHistoryLibrary(lapply(banks, `[[`, "train_check"),
                                          lapply(banks, `[[`, "cal_check"))
    list(simulation = simulation, dna_summary = dna_summary,
         diversity_fit = diversity_fit, finite_state_fit = finite_state_fit,
         ld_fit = ld_fit, ibd_summary = ibd_summary, ibd_fit = ibd_fit, fit_features = fit_features, check_features = check_features,
         history_library = history_library, heldout_library = heldout_library,
         history_fit = comparePopulationHistories(fit_features, history_library),
         history_check = comparePopulationHistories(check_features, heldout_library))
}

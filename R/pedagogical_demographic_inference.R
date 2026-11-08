# Demographic-history comparison from final-sample sequence summaries.
# Training predictions fix centers and scales; independent calibration draws
# give a Monte Carlo discrepancy check. Check fractions describe the declared
# candidate's predictive distribution. Linked loci remain dependent within runs.

historySummaryMatrix <- function(x, name) {
    if (!is.matrix(x) || !is.numeric(x) || nrow(x) < 2L || !ncol(x) ||
        any(!is.finite(x)) || is.null(colnames(x)) ||
        anyDuplicated(colnames(x)) || any(!nzchar(colnames(x)))) {
        stop(name, " must be a finite numeric replicate-by-summary matrix with named columns",
             call. = FALSE)
    }
    x
}

#' Calibrate a comparison of candidate population histories
#'
#' Freeze candidate means and common variance scales using training simulations.
#' Keep independent calibration simulations for the discrepancy check. Generate
#' both banks by calling the original forward engine with the same sequence,
#' mutation, recombination, sample and chromosome design as the observation.
#' @param training Named list of numeric matrices. Rows are independent population
#'   simulation replicates; columns are named final-DNA summaries such as grouped
#'   folded SFS counts per site.
#' @param calibration Separate named list of matrices with the same candidate
#'   names and summary columns, generated with disjoint simulation seeds.
#' @param varianceFloor Positive numerical lower bound on common variance scales.
#' @return A prediction library containing frozen means, variance scales and
#'   independent calibration discrepancies and their per-summary contributions.
#'   Candidate labels describe hypotheses;
#'   no generating history or generation record is read by this calculation.
#' @details The discrepancy is the sum of squared residuals divided by common
#'   simulation variances. Its calibration includes correlations among summaries.
#'   It is a comparison within the declared candidates, not a full demographic
#'   likelihood or proof that the smallest-discrepancy model generated the data.
#' @export
buildHistoryLibrary <- function(training, calibration, varianceFloor = 1e-12) {
    if (!is.list(training) || length(training) < 2L || is.null(names(training)) ||
        anyDuplicated(names(training)) || any(!nzchar(names(training))) ||
        !is.list(calibration) || length(calibration) != length(training) ||
        is.null(names(calibration)) || anyDuplicated(names(calibration)) ||
        !setequal(names(training), names(calibration))) {
        stop("training and calibration require the same named candidate histories",
             call. = FALSE)
    }
    floor <- positiveScalar(varianceFloor, "varianceFloor")
    training <- lapply(training, historySummaryMatrix, name = "training")
    calibration <- lapply(calibration[names(training)], historySummaryMatrix, name = "calibration")
    features <- colnames(training[[1L]])
    if (any(!vapply(c(training, calibration), function(x)
        identical(colnames(x), features), logical(1)))) {
        stop("all prediction matrices must have identical named summary columns", call. = FALSE)
    }
    centers <- do.call(rbind, lapply(training, colMeans))
    variances <- do.call(rbind, lapply(training, function(x) apply(x, 2L, stats::var)))
    scale <- pmax(colMeans(variances), floor)
    contributions <- lapply(seq_along(calibration), function(j) {
        residual <- sweep(calibration[[j]], 2L, centers[j, ], "-")
        sweep(residual^2, 2L, scale, "/")
    })
    names(contributions) <- names(training)
    discrepancy <- lapply(contributions, rowSums)
    list(centers = centers, scale = scale, calibrationDiscrepancy = discrepancy,
         calibrationContributions = contributions,
         training = training, calibration = calibration, summaryNames = features,
         method = "simulation_calibrated_squared_distance")
}

#' Compare final DNA with simulated candidate histories
#'
#' Measure the discrepancy between observed summaries and each frozen prediction.
#' Report how often independent simulations from that candidate are at least as
#' discrepant. A held-out chromosome check uses a separately constructed library
#' with that chromosome design.
#' @param observedSummary Finite named numeric vector of final-DNA summaries.
#' @param library Result of buildHistoryLibrary.
#' @return List with the nearest candidate mean (selected), ranked discrepancies,
#'   simulation check
#'   fractions, frozen predictions, per-summary contributions and the observations.
#'   Check fractions
#'   are not posterior probabilities over the candidate histories.
#' @export
comparePopulationHistories <- function(observedSummary, library) {
    if (!is.numeric(observedSummary) || is.null(names(observedSummary)) ||
        any(!is.finite(observedSummary)) || anyDuplicated(names(observedSummary)) ||
        !setequal(names(observedSummary), library$summaryNames)) {
        stop("observedSummary must match the finite named library summaries", call. = FALSE)
    }
    observation <- observedSummary[library$summaryNames]
    residual <- sweep(library$centers, 2L, observation, "-")
    contributions <- sweep(residual^2, 2L, library$scale, "/")
    distance <- rowSums(contributions)
    checks <- vapply(seq_along(distance), function(j) {
        simulated <- library$calibrationDiscrepancy[[j]]
        (1 + sum(simulated >= distance[j])) / (1 + length(simulated))
    }, numeric(1))
    profile <- data.frame(candidate = rownames(library$centers), discrepancy = distance,
                          compatibility = checks, row.names = NULL)
    profile <- profile[order(profile$discrepancy), , drop = FALSE]
    list(selected = profile$candidate[1L], profile = profile,
         observedSummary = observation, predictions = library$centers,
         scale = library$scale, contributions = contributions, method = library$method)
}

# Closed-form Poisson composite-likelihood fit for Pritchard's specified
# terminal-extension approximation. Columns are derived counts 1,...,n-1.
#' @rdname models_genealogy_controls
#' @export
terminalSpectrumFit <- function(spectra, mutationRate, sequenceLength) {
    if (!is.matrix(spectra) || !is.numeric(spectra) || nrow(spectra) < 1L ||
        ncol(spectra) < 2L || any(!is.finite(spectra)) || any(spectra < 0))
        stop("spectra must contain nonnegative finite counts for n >= 3")
    mu <- positiveScalar(mutationRate, "mutationRate")
    L <- positiveScalar(sequenceLength, "sequenceLength")
    n <- ncol(spectra) + 1L
    B <- nrow(spectra)
    total <- colSums(spectra)
    constant <- sum(total) / (B * 4 * mu * L * sum(1 / seq_len(n - 1L)))
    ancestral <- sum(total[-1L]) / (B * 4 * mu * L * sum(1 / seq.int(2L, n - 1L)))
    extension <- (total[1L] - B * 4 * mu * L * ancestral) / (B * mu * L * n)
    estimate <- if (extension >= 0) c(ancestral, extension) else c(constant, 0)
    expected <- function(N, g) B * mu * L *
        (4 * N / seq_len(n - 1L) + c(n * g, rep(0, n - 2L)))
    list(estimate = setNames(estimate, c("ancestral_Ne", "recent_generations")),
        expected = expected(estimate[1L], estimate[2L]), constantNe = constant,
        constantExpected = expected(constant, 0), observed = total)
}
#' Fit ancestral size and growth duration to chromosome spectra
#'
#' A Poisson composite likelihood uses the demographic SFS expectation.
#' Founder size, small-population duration and present size are supplied;
#' ancestral size and the growth duration are fitted. Sizes are diploid.
#' @param spectra Chromosome rows, derived-count columns 1 through n - 1.
#' @param mutationRate Mutation rate per base per generation.
#' @param sequenceLength Analyzed bases per chromosome, including invariant sites.
#' @param founderSize Number of founding diploid individuals.
#' @param bottleneckGenerations Generations at the founder size.
#' @param presentSize Present diploid population size.
#' @param growthGenerations Candidate positive integer growth durations.
#' @param bounds Lower and upper ancestral sizes for optimization.
#' @return Fitted parameters, expected and observed counts, constant-size fit,
#'   and profile composite log likelihood over the growth durations.
#' @export
fitHistorySpectrum <- function(spectra, mutationRate, sequenceLength,
    founderSize = 20L, bottleneckGenerations = 20L, presentSize = 1000L,
    growthGenerations = seq.int(10L, 80L, 5L), bounds = c(5, 1000)) {
    checkInputs(spectra = list(numericMatrix, minimum = 0),
        mutationRate = positiveScalar, sequenceLength = positiveScalar,
        founderSize = scalarInteger, bottleneckGenerations = scalarInteger,
        presentSize = scalarInteger,
        growthGenerations = list(numericVector, minimum = 1, integer = TRUE, distinct = TRUE),
        bounds = list(numericVector, length = 2L, minimum = .Machine$double.xmin))
    if (bounds[1L] >= bounds[2L] || ncol(spectra) < 2L)
        stop("ordered bounds and at least two SFS bins are required", call. = FALSE)
    n <- ncol(spectra) + 1L
    observed <- colSums(spectra)
    exposures <- nrow(spectra) * sequenceLength
    expectations <- function(ancestral, growth) {
        schedule <- c(rep(founderSize, bottleneckGenerations),
            round(exp(seq(log(founderSize), log(presentSize), length.out = growth + 1L)))[-1L])
        parents <- head(c(founderSize, schedule), -1L)
        expectedHistorySFS(n, c(ancestral, parents), mutationRate, exposures)
    }
    fits <- lapply(growthGenerations, function(g) stats::optimize(function(logN)
        -sum(stats::dpois(observed, expectations(exp(logN), g), log = TRUE)),
        log(bounds)))
    profile <- data.frame(growthGenerations = growthGenerations,
        ancestralSize = exp(vapply(fits, `[[`, 0, "minimum")),
        logLikelihood = -vapply(fits, `[[`, 0, "objective"))
    best <- profile[which.max(profile$logLikelihood), ]
    constant <- sum(observed) / (4 * mutationRate * exposures * sum(1 / seq_len(n - 1L)))
    list(estimate = c(ancestralSize = best$ancestralSize,
                     growthGenerations = best$growthGenerations), profile = profile,
        observed = observed, expected = expectations(best$ancestralSize, best$growthGenerations),
        constantNe = constant,
        constantExpected = expectedSFS(n, constant, mutationRate, exposures))
}

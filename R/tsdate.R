# Public R-list interface to the copied tsdate inside--outside routines in
# Daniel Nachun's tskitR/R/tsdate-interface.R (inside_outside).

#' Date inferred ancestral nodes using mutations
#'
#' Run tsdate's discrete inside--outside algorithm with conditional coalescent
#' priors and edgewise Poisson mutation likelihoods. Node times in the input
#' specify only topological order; sample nodes are fixed at zero. The output
#' contains marginal grid probabilities, moments and constrained node ages.
#' @param ancestry Inferred hpgR_ancestry node/edge/site/mutation tables.
#' @param mutationRate Per-base per-generation mutation rate.
#' @param populationSize Constant diploid effective size used for the prior.
#' @param timepoints Integer grid resolution, or increasing timepoints in generations.
#' @param minBranchLength Minimum positive branch length in generations.
#' @return An hpgR_ancestry object with estimated node times and a dating list
#'   containing the posterior grid, means, variances, mutation counts and method.
#'   Grid probabilities condition on inferred topology, rate and size calibration.
#' @export
dateAncestry <- function(ancestry, mutationRate, populationSize,
                         timepoints = 40L, minBranchLength = 1e-8) {
    ancestryCheckTables(ancestry)
    mutationRate <- positiveScalar(mutationRate, "mutationRate")
    populationSize <- positiveScalar(populationSize, "populationSize")
    minBranchLength <- positiveScalar(minBranchLength, "minBranchLength")
    if (!is.numeric(timepoints) || any(!is.finite(timepoints)) ||
        (length(timepoints) == 1L && (timepoints < 5 || timepoints != floor(timepoints))) ||
        (length(timepoints) > 1L && (length(timepoints) < 5L ||
             timepoints[1L] != 0 || any(diff(timepoints) <= 0))))
        stop("use at least five ordered timepoints, beginning at zero", call. = FALSE)
    x <- tsiSimplify(ancestry)
    x$nodes$time[x$sampleIds + 1L] <- 0
    priors <- tdBuildPriorGrid(x, populationSize, timepoints = timepoints)
    lik <- tdDiscreteLikelihoods(x, priors$timepoints,
        mutation_rate = mutationRate, eps = minBranchLength, log_space = TRUE)
    lik <- tdLikPrecalculate(lik)
    fit <- tdBeliefPropagation(priors, lik)
    fit <- tdBpInsidePass(fit, cache_inside = TRUE)
    fit <- tdBpOutsidePass(fit, standardize = TRUE)
    posterior <- tdNtvToProbabilities(tdNtvForceProbabilitySpace(
        tdNtvStandardize(fit$posterior_grid), .td_LIN_GRID))
    times <- priors$timepoints
    means <- numeric(nrow(x$nodes))
    variance <- numeric(nrow(x$nodes))
    grid <- matrix(0, nrow(x$nodes), length(times))
    grid[x$sampleIds + 1L, 1L] <- 1
    for (node in priors$nonfixed_nodes) {
        probability <- tdNtvGet(posterior, node)
        means[node + 1L] <- sum(probability * times)
        variance[node + 1L] <- sum((means[node + 1L] - times)^2 * probability)
        grid[node + 1L, ] <- probability
    }
    if (any(!is.finite(means)) || any(!is.finite(grid)) || any(abs(rowSums(grid) - 1) > 1e-8))
        stop("node-age posterior did not normalize", call. = FALSE)
    order <- order(x$nodes$time[x$edges$parent + 1L], x$edges$parent)
    x$nodes$time <- tdConstrainAges(means, x$nodes$id %in% x$sampleIds,
        x$edges$parent[order], x$edges$child[order], minBranchLength)
    x$timeUnits <- "generations"
    x$dating <- list(method = "tsdate inside-outside", timepoints = times,
        posterior = grid, mean = means, variance = variance,
        mutationCounts = lik$mut_edges, mutationRate = mutationRate,
        populationSize = populationSize, logMarginalLikelihood = fit$marginal_lik)
    ancestryCheckTables(x)
    x
}

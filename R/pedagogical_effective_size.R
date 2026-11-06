# Effective-size inverses using shared sequence statistics and explicit recurrences.
# Sources: Pritchard, An Owner's Guide to the Human Genome, Sec. 2.1
# (pairwise diversity and mutation--drift balance); Hill (1981), Genetical
# Research 38:209--216 (LD and population size); Handbook of Statistical
# Genomics, Eq. 16.4 (the gamma(r) moment used in the course).
# The finite-state diversity correction below follows directly by conditioning
# on the same/different parental copy, with independent symmetric mutation of
# the parental copies before Wright--Fisher reproduction. It matches the
# SNPMutator preOps timing used by simulatePopulation.
# Assumptions: neutral biallelic sequences, known phase/map/rate, random mating.

populationHaplotypes <- function(haplotypes) {
    x <- genotypeMatrix(haplotypes)
    if (nrow(x) < 2L || any(x > 1L)) {
        stop("haplotypes must have at least two rows and contain only 0 and 1",
             call. = FALSE)
    }
    x
}

populationLengths <- function(chromosome, sequenceLength) {
    ids <- unique(chromosome)
    represented <- vapply(ids, function(id) sum(chromosome == id), integer(1))
    if (!is.numeric(sequenceLength) || any(!is.finite(sequenceLength)) ||
        any(sequenceLength <= 0) || any(sequenceLength != floor(sequenceLength))) {
        stop("sequenceLength must contain positive analyzed base counts", call. = FALSE)
    }
    if (length(sequenceLength) == 1L && length(ids) > 1L) {
        if (sequenceLength != sum(represented)) {
            stop("supply a named analyzed length for each chromosome", call. = FALSE)
        }
        sequenceLength <- represented
    } else if (length(sequenceLength) > 1L) {
        if (is.null(names(sequenceLength)) || anyDuplicated(names(sequenceLength)) ||
            !setequal(names(sequenceLength), ids)) {
            stop("sequenceLength names must match the chromosomes", call. = FALSE)
        }
        sequenceLength <- sequenceLength[ids]
    }
    if (length(sequenceLength) != length(ids) || any(sequenceLength < represented)) {
        stop("analyzed length cannot be smaller than the represented site count",
             call. = FALSE)
    }
    stats::setNames(as.numeric(sequenceLength), ids)
}

#' Summarize the DNA in a present-day population sample
#'
#' Count pairwise differences and the folded site frequency spectrum (SFS).
#' Rows are sampled haploid sequences and columns are biallelic sites. Invariant
#' bases contribute to the analyzed length, including bases omitted from a
#' variant-only matrix. Allele 1 need not be the derived allele because the SFS
#' is folded using the minor allele count.
#' @param haplotypes Matrix with haploid sequences in rows and 0/1 sites in columns.
#' @param chromosome Chromosome labels, one per represented site.
#' @param sequenceLength Analyzed bases, including invariants. For several
#'   chromosomes supply named lengths, or keep the default when every base is
#'   represented in the matrix.
#' @param sfsGroups Named list of minor allele counts to pool. The default groups
#'   singletons, counts 2 through floor(n/4), and the remaining polymorphic counts.
#' @return List containing pairwise diversity per site, segregating-site counts,
#'   folded SFS, grouped SFS, chromosome summaries and the explicit sample sizes.
#' @export
summarizePopulationSample <- function(haplotypes, chromosome = NULL,
                                      sequenceLength = ncol(haplotypes),
                                      sfsGroups = NULL) {
    x <- populationHaplotypes(haplotypes)
    n <- nrow(x)
    if (is.null(chromosome)) chromosome <- rep("1", ncol(x))
    if (length(chromosome) != ncol(x) || anyNA(chromosome)) {
        stop("chromosome must label every represented site", call. = FALSE)
    }
    chromosome <- as.character(chromosome)
    lengths <- populationLengths(chromosome, sequenceLength)
    k <- colSums(x)
    minor <- pmin(k, n - k)
    maxMinor <- floor(n / 2)
    if (is.null(sfsGroups)) {
        low <- seq_len(maxMinor)
        sfsGroups <- list(singleton = 1L,
            low_frequency = low[low >= 2L & low <= floor(n / 4)],
            common = low[low > max(1L, floor(n / 4))])
    }
    if (!is.list(sfsGroups) || !length(sfsGroups) || is.null(names(sfsGroups)) ||
        any(!nzchar(names(sfsGroups))) || anyDuplicated(names(sfsGroups))) {
        stop("sfsGroups must be a named list", call. = FALSE)
    }
    groups <- unlist(sfsGroups, use.names = FALSE)
    if (length(groups) && (!is.numeric(groups) || any(!is.finite(groups)) ||
        any(groups != floor(groups)) || any(groups < 1L | groups > maxMinor) ||
        anyDuplicated(groups))) {
        stop("SFS groups must contain nonoverlapping polymorphic minor counts",
             call. = FALSE)
    }
    byChromosome <- lapply(names(lengths), function(id) {
        selected <- chromosome == id
        sfs <- tskitSfs(k[selected], n, polarised = FALSE)
        names(sfs) <- as.character(seq_len(maxMinor))
        differences <- sum(tskitDiversity(k[selected], n))
        list(chromosome = id, sequenceLength = lengths[[id]],
             pairwiseDifferences = differences,
             piPerSite = differences / lengths[[id]],
             segregatingSites = sum(minor[selected] > 0), foldedSfs = sfs,
             groupedSfs = vapply(sfsGroups, function(group) sum(sfs[group]), numeric(1)))
    })
    names(byChromosome) <- names(lengths)
    sfs <- Reduce(`+`, lapply(byChromosome, `[[`, "foldedSfs"))
    grouped <- vapply(sfsGroups, function(group) sum(sfs[group]), numeric(1))
    totalLength <- sum(lengths)
    differences <- sum(vapply(byChromosome, `[[`, numeric(1), "pairwiseDifferences"))
    list(nHaplotypes = n, sequenceLength = totalLength,
         piPerSite = differences / totalLength,
         pairwiseDifferences = differences, segregatingSites = sum(sfs),
         foldedSfs = sfs, groupedSfs = grouped,
         groupedSfsPerSite = grouped / totalLength,
         sfsGroups = sfsGroups, byChromosome = byChromosome)
}

#' Estimate effective size using diversity and a mutation rate
#'
#' Compute the constant-size neutral equilibrium inverse pi/(4 mu). The optional
#' symmetric-parental model accounts for recurrent 0/1 mutation before copying.
#' Under a changing size history these are diversity-based effective summaries,
#' not estimates of the current census or a uniquely identified size trajectory.
#' @param haplotypes Matrix of sampled 0/1 haplotypes.
#' @param mutationRate Mutation probability per site per generation.
#' @param sequenceLength Analyzed bases, including invariant bases.
#' @param mutationModel One of infinite_sites or symmetric_parental.
#' @return List with Ne, observed piPerSite, mutationRate, analyzed length,
#'   sample size and the model used for inversion.
#' @details With symmetric parental mutation let a=2u(1-u), b=(1-2u)^2,
#'   c=1/(2Ne). Expected pair diversity follows H'=(1-c)(a+bH). At equilibrium
#'   Ne=(a+bH)/(2a(1-2H)). This includes the finite parental pool in the
#'   SNPMutator-before-reproduction model. The infinite-sites inverse is the
#'   low-mutation, large-population approximation. Values H at least 1/2 give
#'   no finite estimate under the symmetric-parental equilibrium inverse.
#'   Here H is diversity between distinct sampled copies. A full population's
#'   raw value 2p(1-p) differs by the finite-pool factor (2N-1)/(2N).
#' @export
estimateDiversitySize <- function(haplotypes, mutationRate,
                                  sequenceLength = ncol(haplotypes),
                                  mutationModel = c("infinite_sites", "symmetric_parental")) {
    x <- populationHaplotypes(haplotypes)
    mu <- positiveScalar(mutationRate, "mutationRate")
    if (mu >= .5) stop("mutationRate must be smaller than 1/2", call. = FALSE)
    length <- positiveScalar(sequenceLength, "sequenceLength")
    if (length < ncol(x)) stop("sequenceLength is smaller than the site count", call. = FALSE)
    model <- match.arg(mutationModel)
    k <- colSums(x)
    n <- nrow(x)
    H <- sum(tskitDiversity(k, n)) / length
    if (model == "infinite_sites") {
        Ne <- H / (4 * mu)
    } else {
        a <- 2 * mu * (1 - mu)
        b <- (1 - 2 * mu)^2
        Ne <- if (H < .5) (a + b * H) / (2 * a * (1 - 2 * H)) else Inf
    }
    list(Ne = Ne, piPerSite = H, mutationRate = mu, sequenceLength = length,
         nHaplotypes = n, mutationModel = model)
}

inferencePairIndices <- function(m, maxPairs) {
    count <- m * (m - 1) / 2
    if (count <= maxPairs) return(t(utils::combn(m, 2L)))
    if (count <= 1e6) {
        pairs <- t(utils::combn(m, 2L))
        return(pairs[sample.int(nrow(pairs), maxPairs), , drop = FALSE])
    }
    # Uniform unordered pairs without replacement, with bounded memory.
    key <- numeric()
    while (length(key) < maxPairs) {
        draw <- matrix(sample.int(m, 4L * (maxPairs - length(key)), replace = TRUE), ncol = 2L)
        draw <- draw[draw[, 1L] != draw[, 2L], , drop = FALSE]
        key <- unique(c(key, pmin(draw[, 1L], draw[, 2L]) +
                            m * (pmax(draw[, 1L], draw[, 2L]) - 1)))
    }
    key <- key[seq_len(maxPairs)]
    cbind(as.integer((key - 1) %% m + 1), as.integer((key - 1) %/% m + 1))
}

#' Estimate effective size using linked sites and a recombination map
#'
#' Invert the weak-LD equilibrium moment E(rho^2) approximately gamma(r)/Ne,
#' where gamma(r)=((1-r)^2+r^2)/(2r(2-r)). The leading finite-sample correction
#' subtracts 1/n for n known-phase haplotype draws. Supply one predetermined
#' randomly chosen homolog per sampled individual when both homologs are paired
#' in the simulator output. The approximation is calibrated using independent
#' population simulations with the same observation process.
#' @param haplotypes Matrix of sampled 0/1 haplotypes, one independently selected
#'   homolog per sampled individual for the course sampling design.
#' @param geneticPositions Site positions in Morgans, one per column.
#' @param chromosome Chromosome labels, one per column.
#' @param distanceBreaks Increasing genetic distances defining bins in Morgans.
#' @param minMinorAlleleFrequency Minor allele frequency filter in the sample.
#' @param maxPairs Maximum number of polymorphic site pairs to examine.
#' @param seed Seed for uniform pair subsampling. The caller's RNG state is restored.
#' @return List containing pooled Ne, a data frame of distance-bin estimates,
#'   individual pair statistics and the explicit sample-size correction. Negative
#'   corrected moments give Inf and do not provide a finite size estimate.
#' @details The Haldane fraction r=(1-exp(-2d))/2 converts map distance d into
#'   a two-locus recombination fraction. Within each bin the estimator is
#'   mean(gamma(r))/(mean(rho^2)-1/n). The sample correction, weak-LD formula
#'   and equilibrium approximation can be inaccurate for rare variants, tight
#'   linkage, small samples, structure and changing size. Pairs share markers
#'   and chromosomes; pair counts are not independent replication counts.
#' @export
estimateLdSize <- function(haplotypes, geneticPositions, chromosome = NULL,
                           distanceBreaks = c(.03, .07, .12, .2),
                           minMinorAlleleFrequency = .1, maxPairs = 50000L,
                           seed = 1L) {
    x <- populationHaplotypes(haplotypes)
    n <- nrow(x)
    if (n < 4L) stop("at least four sampled haplotypes are required", call. = FALSE)
    if (!is.numeric(geneticPositions) || length(geneticPositions) != ncol(x) ||
        any(!is.finite(geneticPositions))) {
        stop("geneticPositions must give one finite Morgan position per site", call. = FALSE)
    }
    if (is.null(chromosome)) chromosome <- rep("1", ncol(x))
    if (length(chromosome) != ncol(x) || anyNA(chromosome)) {
        stop("chromosome must label each site", call. = FALSE)
    }
    if (!is.numeric(distanceBreaks) || length(distanceBreaks) < 2L ||
        any(!is.finite(distanceBreaks)) || any(distanceBreaks < 0) ||
        any(diff(distanceBreaks) <= 0)) {
        stop("distanceBreaks must be increasing nonnegative Morgan distances", call. = FALSE)
    }
    cutoff <- positiveScalar(minMinorAlleleFrequency, "minMinorAlleleFrequency")
    if (cutoff >= .5) stop("minor allele frequency cutoff must be below 1/2", call. = FALSE)
    checkInputs(maxPairs = scalarInteger, seed = list(scalarInteger, minimum = 0L))
    maxPairs <- as.integer(maxPairs)
    seed <- as.integer(seed)
    restore <- rngRestorer()
    on.exit(restore(), add = TRUE)
    set.seed(seed)
    p <- colMeans(x)
    selected <- which(pmin(p, 1 - p) >= cutoff)
    if (length(selected) < 2L) stop("fewer than two sites pass the frequency filter", call. = FALSE)
    indices <- inferencePairIndices(length(selected), maxPairs)
    first <- selected[indices[, 1L]]
    second <- selected[indices[, 2L]]
    distance <- abs(geneticPositions[first] - geneticPositions[second])
    bin <- findInterval(distance, distanceBreaks, rightmost.closed = TRUE)
    keep <- chromosome[first] == chromosome[second] & distance > 0 &
        bin >= 1L & bin < length(distanceBreaks)
    first <- first[keep]; second <- second[keep]
    distance <- distance[keep]; bin <- bin[keep]
    if (!length(first)) stop("no within-chromosome pairs fall in the distance bins", call. = FALSE)
    joint <- vapply(seq_along(first), function(j)
        mean(x[, first[j]] * x[, second[j]]), numeric(1))
    D <- tskitD(joint, p[first], p[second])
    rho2 <- tskitR2(joint, p[first], p[second])
    r <- -expm1(-2 * distance) / 2
    gamma <- ((1 - r)^2 + r^2) / (2 * r * (2 - r))
    correction <- 1 / n
    estimate <- function(index) {
        adjusted <- mean(rho2[index]) - correction
        c(meanRho2 = mean(rho2[index]), correctedRho2 = adjusted,
          meanGamma = mean(gamma[index]),
          Ne = if (adjusted > 0) mean(gamma[index]) / adjusted else Inf)
    }
    bins <- lapply(seq_len(length(distanceBreaks) - 1L), function(j) {
        index <- which(bin == j)
        values <- if (length(index)) estimate(index) else
            c(meanRho2 = NA_real_, correctedRho2 = NA_real_, meanGamma = NA_real_, Ne = NA_real_)
        data.frame(lowerMorgan = distanceBreaks[j], upperMorgan = distanceBreaks[j + 1L],
                   pairs = length(index), meanDistance = if (length(index)) mean(distance[index]) else NA_real_,
                   meanRecombinationFraction = if (length(index)) mean(r[index]) else NA_real_,
                   as.list(values), row.names = NULL, check.names = FALSE)
    })
    pooled <- estimate(seq_along(first))
    list(Ne = unname(pooled["Ne"]), bins = do.call(rbind, bins),
         pairs = data.frame(first, second, distanceMorgan = distance,
                            recombinationFraction = r, rho2, gamma, bin),
         nHaplotypes = n, sitesPassingFilter = length(selected),
         samplingCorrection = correction, minMinorAlleleFrequency = cutoff,
         pooled = pooled, model = "weak_LD_equilibrium")
}

#' Two-locus LD and temporal drift inversion
#'
#' LD is f_AB - p_A p_B. Temporal inversion pools replicate frequency changes that share the same initial p, then equates their variance to p(1-p)/(2Ne).
#'
#' @name models_frequency_controls
#' @param f Four haplotype frequencies ordered AB, Ab, aB, ab.
#' @param parent_frequency Equal, interior parental frequencies across at least two independent transitions.
#' @param offspring_frequency Matching offspring frequencies, one per transition.
#' @returns ld returns scalar D. estimateFromGenerations returns parental_p, transitions, variance, c_V and Ne_V; offspring records must have positive variance.
#' @examples
#' ld(c(.3,.2,.2,.3))
#' estimateFromGenerations(rep(.5,4), c(.4,.5,.6,.5))
NULL

# Frequency-class drift inversion and two-locus D from the course's
# existing computations. Sources: Pritchard 2.1/2.3 and canonical Notes 8.8.

#' @rdname models_frequency_controls
#' @export
ld <- function(f) tskitD(f[1], f[1] + f[2], f[1] + f[3])

#' @rdname models_frequency_controls
#' @export
estimateFromGenerations <- function(parent_frequency, offspring_frequency) {
    stopifnot(length(parent_frequency) == length(offspring_frequency),
              length(parent_frequency) > 1L, length(unique(parent_frequency)) == 1L)
    p <- parent_frequency[1]
    stopifnot(p > 0, p < 1)
    changes <- offspring_frequency - parent_frequency
    s2 <- var(changes)
    stopifnot(s2 > 0)
    c_V <- s2 / (p * (1 - p))
    c(parental_p = p, transitions = length(changes), variance = s2,
      c_V = c_V, Ne_V = 1 / (2 * c_V))
}

#' Small genealogy and diversity experiments
#'
#' These recipes call shared ancestry and mutation functions with explicit experimental designs. Simulation truth is stored separately from observations. The terminal-extension model adds a fixed duration to terminal branches and fits its spectrum by a Poisson composite likelihood.
#'
#' @name models_genealogy_controls
#' @param seed Nonnegative integer random seed.
#' @param spectra Region-by-derived-count matrix with n-1 columns of nonnegative spectrum counts.
#' @param mutationRate Positive mutation rate per base per generation.
#' @param sequenceLength Positive analyzed length per independent region in bases.
#' @param N Positive integer diploid population size.
#' @param mu Nonnegative mutation rate per base per generation.
#' @param replicates Positive integer number of independent simulated populations or genealogies.
#' @param base_seed Nonnegative integer starting seed for region or replicate assignments.
#' @param tree Compressed genealogy returned by the tree constructors.
#' @param selected Distinct sample indices retained in the subtree.
#' @param Ne Positive diploid effective size for a constant-size Kingman model.
#' @param n Integer number of sampled haploid copies, at least two.
#' @param extension Nonnegative terminal extension in generations before the sample.
#' @param regions Integer number of independent regions, at least two.
#' @param bootstrap Integer number of region-resampling bootstrap replicates, at least two.
#' @param generations Nonnegative integer number of forward reproductive updates.
#' @param G Binary haplotype matrix: sampled copies in rows, derived variants in columns.
#' @returns matchedDraws returns replicate results and spectra. recordedSubset returns a nested tree and ages. sfs returns counts for derived-copy classes 1 through n-1. founderDiversityExample returns replicate nonIBD and diversityRatio arrays, their summary and truth. terminalSpectrumFit returns estimate, expected, constantNe, constantExpected and observed. terminalExtensionExample returns observed regions, spectra, fit, bootstrap intervals and truth. specifiedLocalTreeExample returns observed, calibration and truth with two specified local genealogies.
#' @examples
#' founderDiversityExample(N=4, generations=2, replicates=2, seed=42)$summary
#' sfs(rbind(c(0,1),c(1,1),c(0,0)))
NULL

# Sources: Pritchard 2.2; canonical Notes 8.10--8.12.

#' @rdname models_genealogy_controls
#' @export
matchedDraws <- function(N, mu, replicates = 600L, base_seed) {
    rows <- list(); spectra <- matrix(0, 2 * replicates, 11)
    for (model in 1:2) for (i in seq_len(replicates)) {
        a <- pgRegionSeed(base_seed + model * 10000L, i)
        m <- pgRegionSeed(base_seed + model * 10000L + 500000L, i, mutation = TRUE)
        genealogy <- if (model == 1L) pgDiscreteGenealogy(rep(N, 61L), ancestry_seed = a)
                     else pgKingmanTree(N, ancestry_seed = a)
        summary <- pgTreeSummaries(genealogy)
        mutations <- pgMutateTree(genealogy, 1000, mu, mutation_seed = m)
        row <- (model - 1L) * replicates + i
        spectra[row, ] <- mutations$sfs
        weights <- 2 * (1:11) * (12 - (1:11)) / (12 * 11)
        internal_ages <- genealogy$nodes$time[is.na(genealogy$nodes$sample_id)]
        first_age <- min(internal_ages)
        first_parents <- genealogy$nodes$id[genealogy$nodes$time == first_age]
        first_loss <- sum(table(genealogy$edges$parent)[as.character(first_parents)] - 1L)
        rows[[row]] <- data.frame(N = N, model = c("Discrete WF", "Kingman")[model],
            replicate = i, pair_time = summary$pair_time, TMRCA = summary$TMRCA,
            branch_length = summary$total_branch_length, pi_region = sum(mutations$sfs * weights),
            S = sum(mutations$sfs),
            pair_differences = sum(mutations$haplotypes[1L, ] != mutations$haplotypes[2L, ]),
            first_loss = first_loss,
            multifurcation = summary$max_children > 2L, ancestry_seed = a, mutation_seed = m)
    }
    list(results = do.call(rbind, rows), sfs = spectra)
}

# Founder-label and allele summaries of copying records.
#' @rdname models_genealogy_controls
#' @export
founderDiversityExample <- function(N = 40L, generations = 60L,
                                    replicates = 200L, seed = 814906L) {
    checkInputs(N = scalarInteger, generations = list(scalarInteger, minimum = 0L),
                replicates = list(scalarInteger, minimum = 2L), seed = scalarInteger)
    N <- as.integer(N)
    generations <- as.integer(generations)
    replicates <- as.integer(replicates)
    seed <- as.integer(seed)
    copies <- 2L * N
    initial <- cbind(seq_len(copies), rep(0:1, each = N))
    H0 <- 2 * N^2 / (copies * (copies - 1))
    nonIBD <- ratio <- matrix(0, replicates, generations + 1L)
    for (i in seq_len(replicates)) {
        record <- simulateForward(initial, rep(N, generations + 1L), seed + i - 1L)
        for (t in seq_along(record$states)) {
            state <- record$states[[t]]
            counts <- tabulate(state[, 1L], copies)
            nonIBD[i, t] <- 1 - sum(counts * (counts - 1)) / (copies * (copies - 1))
            k <- sum(state[, 2L])
            ratio[i, t] <- 2 * k * (copies - k) / (copies * (copies - 1)) / H0
        }
    }
    list(summary = data.frame(generation = 0:generations,
        nonIBD = colMeans(nonIBD), diversityRatio = colMeans(ratio),
        expectation = (1 - 1 / copies)^(0:generations)),
        nonIBD = nonIBD, diversityRatio = ratio,
        truth = list(N = N, H0 = H0, replicates = replicates,
                     seeds = seed + seq_len(replicates) - 1L))
}

#' @rdname models_genealogy_controls
#' @export
terminalExtensionExample <- function(Ne = 80, n = 12L, mutationRate = 2e-5,
                                     sequenceLength = 1000, extension = 40,
                                     regions = 48L, bootstrap = 120L, seed = 814903L) {
    checkInputs(n = list(scalarInteger, minimum = 3L),
                regions = list(scalarInteger, minimum = 2L),
                bootstrap = list(scalarInteger, minimum = 2L), seed = scalarInteger,
                extension = list(positiveScalar, zero = TRUE))
    n <- as.integer(n)
    regions <- as.integer(regions)
    bootstrap <- as.integer(bootstrap)
    seed <- as.integer(seed)
    trees <- lapply(seq_len(regions), function(j) {
        tree <- pgKingmanTree(Ne, n, pgRegionSeed(seed, j))
        nodes <- tree$nodes[c("id", "time")]
        internal <- !nodes$id %in% tree$sample_nodes
        nodes$time[internal] <- nodes$time[internal] + extension
        pgMakeTree(nodes, tree$edges[c("parent", "child")], tree$sample_nodes,
                     tree$sample_ids, "specified terminal extension")
    })
    observed <- lapply(seq_along(trees), function(j) pgMutateTree(trees[[j]],
        sequenceLength, mutationRate, pgRegionSeed(seed + 100000L, j, TRUE)))
    spectra <- t(vapply(observed, function(x) sfs(x$haplotypes), integer(n - 1L)))
    fit <- terminalSpectrumFit(spectra, mutationRate, sequenceLength)
    draws <- pgWithSeed(seed + 1L, function() t(replicate(bootstrap, {
        rows <- sample.int(regions, regions, replace = TRUE)
        terminalSpectrumFit(spectra[rows, , drop = FALSE], mutationRate, sequenceLength)$estimate
    })))
    interval <- apply(draws, 2L, stats::quantile, c(.025, .975))
    list(observed = observed, spectra = spectra, fit = fit, bootstrap = draws,
        interval = interval, truth = list(Ne = Ne, extension = extension,
            mutationRate = mutationRate, sequenceLength = sequenceLength, n = n,
            regions = regions, seed = seed, trees = trees))
}

# Two local genealogies specified before mutation.
#' @rdname models_genealogy_controls
#' @export
specifiedLocalTreeExample <- function(seed = 814905L) {
    checkInputs(seed = scalarInteger)
    seed <- as.integer(seed)
    make <- function(pair, time) {
        other <- setdiff(1:3, pair)
        pgMakeTree(data.frame(id = 1:7, time = c(rep(0, 4), time, 400, 800)),
            data.frame(parent = c(5L, 5L, 6L, 7L, 6L, 7L),
                       child = c(pair, other, 4L, 5L, 6L)),
            1:4, paste0("S", 1:4), "specified local genealogy")
    }
    trees <- list(left = make(1:2, 120), right = make(c(1L, 3L), 180))
    data <- lapply(seq_along(trees), function(i)
        pgMutateTree(trees[[i]], 50000, 1e-6, mutation_seed = seed + i - 1L))
    position <- c(data[[1L]]$positions_bp, data[[2L]]$positions_bp + 50000)
    order <- order(position)
    observed <- list(haplotypes = cbind(data[[1L]]$haplotypes,
        data[[2L]]$haplotypes)[, order, drop = FALSE],
        positions = position[order], sequenceLength = 100000)
    list(observed = observed, calibration = list(mutationRate = 1e-6),
        truth = list(breakpoint = 50000, TMRCA = 800, trees = trees,
                     mutationSeeds = seed + 0:1))
}

#' A recorded Kingman genealogy with branch mutations
#'
#' Event records for the genealogy animation, with diploid Ne = 1000, mutation rate 1e-8 per base per generation and length 100000 bases. The recorded exponential mergers and branch mutations can be replayed through the shared ancestry and mutation kernels.
#'
#' @name models_coalescent_example
#' @returns List containing seed, n, Ne, mu, L, nodes, edges, merger events and mutations with carrier sets.
#' @examples
#' example <- coalescentExample()
#' example$events[[1]]$wait
NULL

# Source model: Pritchard 2.2, pp.91--94 and 102--103, Fig.2.35;
# Laruson and Reed (2021), Sec.9.1, pp.178--180.
#' @rdname models_coalescent_example
#' @export
coalescentExample <- function() {
    readRDS(system.file("examples", "coalescent.rds", package = "hpgR", mustWork = TRUE))
}


#' Extract a site's genealogy from recorded population transmissions
#'
#' Convert the realized parental-copy indices into a finite-copy ancestry record.
#' Founder sampling maps the first reproductive pool to its original ancestral
#' copies without adding a generation. The DNA and transmission draws are kept.
#' @param population Result of simulatePopulation with recordSites supplied.
#' @param site Recorded marker-column index. Defaults to the first recorded site.
#' @param samples Rows of the observed haplotype matrix to trace.
#' @return A copy record for one site, suitable for pgCompressTree. It includes
#'   population nodes and edges, sampled-copy IDs, sample labels, the selected
#'   marker and root-completion status within the recorded burn-in. No older
#'   unrecorded ancestry is added.
#' @export
populationGenealogy <- function(population, site = NULL,
                                samples = seq_len(nrow(population$observed$haplotypes))) {
    transmission <- population$records$copyParents
    if (is.null(transmission)) stop("recordSites is required when simulating the population",
                                   call. = FALSE)
    if (is.null(site)) site <- transmission$sites[1L]
    checkInputs(site = scalarInteger,
        samples = list(numericVector, minimum = 1,
            maximum = nrow(population$observed$haplotypes), integer = TRUE, distinct = TRUE))
    column <- match(site, transmission$sites)
    if (is.na(column)) stop("site was not recorded", call. = FALSE)
    truth <- population$truth
    burnin <- truth$burnin
    sizes <- c(rep(truth$ancestralSize, burnin + 1L), truth$sizes[-1L])
    generations <- seq_along(sizes) - 1L - burnin
    offsets <- cumsum(c(0L, head(2L * sizes, -1L)))
    pools <- lapply(seq_along(sizes), function(i)
        offsets[i] + seq_len(2L * sizes[i]))
    nodes <- data.frame(id = seq_len(2L * sum(sizes)),
        generation = rep(generations, 2L * sizes),
        time = rep(tail(generations, 1L) - generations, 2L * sizes))
    founderRows <- as.vector(rbind(2L * truth$founder_indices - 1L,
                                   2L * truth$founder_indices))
    edges <- lapply(seq_along(transmission$parents), function(i) {
        parent <- transmission$parents[[i]][, column]
        if (i == burnin + 1L && length(founderRows)) parent <- founderRows[parent]
        data.frame(left = 0, right = 1, parent = pools[[i]][parent], child = pools[[i + 1L]])
    })
    rows <- as.vector(rbind(2L * truth$sample_indices - 1L, 2L * truth$sample_indices))
    sampleNodes <- tail(pools, 1L)[[1L]][rows[samples]]
    record <- list(nodes = nodes, edges = do.call(rbind, edges), length_bp = 1,
        sample_nodes = sampleNodes, sample_indices = rows[samples],
        sample_generation = tail(generations, 1L),
        size_history = data.frame(generation = generations, N = sizes),
        ancestral_N = truth$ancestralSize, extension_generations = 0L,
        ancestry_seed = truth$seed, map_length = 0,
        sample_ids = paste0("H", samples), marker = transmission$markers[column, , drop = FALSE])
    record$local_trees <- pgLocalTrees(record)
    record$root_complete <- all(vapply(record$local_trees, `[[`, logical(1), "complete"))
    pgValidateAncestry(record, require_complete = FALSE)
    record
}
#' Trace a fixed sample at recorded genomic positions
#' @param design Population-history design; defaults to populationHistoryDesign.
#' @param sites Marker columns at which parental copying is recorded.
#' @param samples Sampled haplotype rows, in a fixed order.
#' @inheritParams simulatePopulationHistory
#' @return The original simulation, compressed point trees, and four-copy
#'   projections. Recording does not change the simulated DNA.
#' @export
populationAncestryExample <- function(design = populationHistoryDesign()$design,
    sites = unique(c(5000L, as.integer(rep((0:7) * 10000L, each = 2L) +
        rep(c(125L, 9875L), 8L)), 2500L, 7500L)), samples = 1:12,
    recordEvery = 10000L) {
    specification <- populationHistoryDesign()
    simulation <- simulatePopulationHistory(specification$candidates$recent_founding,
        design, design$observation_seed, recordEvery = recordEvery, recordSites = sites)
    trees <- lapply(sites, function(site) {
        populationTree(simulation, site, samples)
    })
    names(trees) <- as.character(sites)
    list(simulation = simulation, sites = sites, trees = trees,
         four = lapply(trees, recordedSubset))
}

#' Compare discrete and coalescent genealogies for one size history
#' @param design Population-history sequence and sampling design.
#' @param replicates Independent draws per ancestry model.
#' @param scale Multiply sizes and durations by this integer; divide mutation rate by it.
#' @param baseSeed First seed, fixed before drawing the trees.
#' @return Repeated tree ages, finite-site sequence summaries, and their spectra.
#' @export
matchedHistoryDraws <- function(design = populationHistoryDesign()$design,
    replicates = 200L, scale = 1L, baseSeed = 7814900L) {
    checkInputs(replicates = list(scalarInteger, minimum = 2L),
        scale = scalarInteger, baseSeed = scalarInteger)
    candidate <- populationHistoryDesign()$candidates$recent_founding
    schedule <- populationSizeSchedule(candidate, design)
    parents <- head(c(candidate$founder, schedule[-1L]), -1L)
    parents <- rep(parents * scale, each = scale)
    n <- 12L
    L <- design$sites_per_chromosome
    mu <- design$mutation_rate / scale
    rows <- vector("list", 2L * replicates)
    spectra <- matrix(0L, 2L * replicates, n - 1L)
    for (i in seq_along(rows)) {
        discrete <- i <= replicates
        tree <- if (discrete) pgDiscreteGenealogy(c(parents, design$present_size * scale),
            n, baseSeed + i, ancestralSize = design$ancestral_size * scale) else
            pgKingmanTree(design$ancestral_size * scale, n, baseSeed + i,
                          populationSizes = rev(parents))
        summary <- pgTreeSummaries(tree)
        ancestry <- structure(list(nodes = data.frame(id = tree$nodes$id - 1L,
            time = tree$nodes$time, flags = as.integer(tree$nodes$id %in% tree$sample_nodes)),
            edges = data.frame(left = 0, right = L,
                parent = tree$edges$parent - 1L, child = tree$edges$child - 1L),
            sampleIds = tree$sample_nodes - 1L, sequenceLength = L,
            timeUnits = "generations"), class = c("hpgR_ancestry", "list"))
        dna <- ancestryData(mutateAncestry(ancestry, mu, baseSeed + 100000L + i,
                                          model = "binary"))$observed$haplotypes
        spectra[i, ] <- sfs(dna)
        # Pair clock uses observed differences, including recurrent hits.
        rows[[i]] <- data.frame(model = if (discrete) "Discrete WF" else "Kingman",
            scale = scale, pair_time = summary$pair_time, TMRCA = summary$TMRCA,
            branch_length = summary$total_branch_length,
            S = sum(spectra[i, ]), pi_region = sum(tskitDiversity(colSums(dna), n)),
            pair_differences = sum(dna[1L, ] != dna[2L, ]),
            multifurcation = summary$max_children > 2L,
            ancestry_seed = baseSeed + i, mutation_seed = baseSeed + 100000L + i)
    }
    list(results = do.call(rbind, rows), sfs = spectra, parentalSizes = parents,
        n = n, length = L, mutationRate = mu,
        expectation = expectedHistorySFS(n, c(design$ancestral_size * scale, parents), mu, L))
}
#' Read a point tree from recorded parental copying
#' @inheritParams populationGenealogy
#' @return A compressed tree with actual integer node ages and sampled-copy labels.
#' @export
populationTree <- function(population, site = NULL,
    samples = seq_len(nrow(population$observed$haplotypes))) {
    transmission <- population$records$copyParents
    if (is.null(transmission)) stop("recordSites is required when simulating the population",
                                   call. = FALSE)
    if (is.null(site)) site <- transmission$sites[1L]
    checkInputs(site = scalarInteger,
        samples = list(numericVector, minimum = 1,
            maximum = nrow(population$observed$haplotypes), integer = TRUE, distinct = TRUE))
    column <- match(site, transmission$sites)
    if (is.na(column)) stop("site was not recorded", call. = FALSE)
    truth <- population$truth
    n <- length(samples)
    rows <- as.vector(rbind(2L * truth$sample_indices - 1L,
                           2L * truth$sample_indices))[samples]
    founderRows <- as.vector(rbind(2L * truth$founder_indices - 1L,
                                   2L * truth$founder_indices))
    active <- seq_len(n)
    nodes <- data.frame(id = active, time = 0)
    edges <- list()
    total <- length(transmission$parents)
    for (age in seq_len(total)) {
        generation <- total - age + 1L
        rows <- transmission$parents[[generation]][rows, column]
        if (generation == truth$burnin + 1L && length(founderRows))
            rows <- founderRows[rows]
        groups <- split(seq_along(rows), rows)
        nextNodes <- integer(length(groups))
        for (j in seq_along(groups)) {
            group <- groups[[j]]
            if (length(group) == 1L) nextNodes[j] <- active[group] else {
                id <- nrow(nodes) + 1L
                nodes <- rbind(nodes, data.frame(id = id, time = age))
                edges[[length(edges) + 1L]] <- data.frame(parent = id, child = active[group])
                nextNodes[j] <- id
            }
        }
        active <- nextNodes
        rows <- as.integer(names(groups))
        if (length(active) == 1L) break
    }
    if (length(active) > 1L)
        stop("recorded burn-in does not contain the sample MRCA", call. = FALSE)
    result <- pgMakeTree(nodes, do.call(rbind, edges), seq_len(n),
        paste0("H", samples), "recorded point tree from the shared forward population")
    result$marker <- transmission$markers[column, , drop = FALSE]
    result
}

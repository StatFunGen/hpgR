# Based on simuPOP's SNPMutator, RandomParentChooser, Recombinator and
# RandomSampler (Bo Peng and contributors, GPL-3.0-or-later).

simupopWithSeed <- function(seed, draw) {
    restore <- rngRestorer()
    on.exit(restore(), add = TRUE)
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
             sample.kind = "Rejection")
    draw()
}

simupopPopulation <- function(sizes, loci, mutationRate, reverseMutationRate,
                                sampleSize, burnin, founderSize, recordEvery,
                                seed, markers, fractions, recordSites = NULL) {
    simupopWithSeed(seed, function() {
        haplotypes <- matrix(0L, nrow = 2L * sizes[1L], ncol = sum(loci))
        copyParents <- list()
        nextGeneration <- function(h, size) {
            h <- simupopMutate(h, mutationRate, reverseMutationRate)
            parents <- simupopParents(nrow(h) %/% 2L, size)
            transmission <- simupopMeiosis(h, parents, loci, fractions,
                binary = TRUE, returnEvents = !is.null(recordSites))
            if (is.null(recordSites)) return(transmission)
            # Read the parental row from the realized Recombinator events.
            # Recording consumes no random draws and leaves the DNA unchanged.
            copies <- vapply(recordSites, function(site) {
                chromosome <- markers$chromosome[site]
                local <- site - sum(loci[seq_len(chromosome - 1L)])
                switches <- transmission$events$switches[[chromosome]]
                parity <- tabulate(switches[switches[, 2L] < local, 1L],
                                   nbins = length(parents)) %% 2L
                2L * parents - 1L + bitwXor(transmission$events$starts[, chromosome], parity)
            }, integer(length(parents)))
            copyParents[[length(copyParents) + 1L]] <<- copies
            transmission$haplotypes
        }
        summarize <- function(generation, h, p = colMeans(h)) {
            c(generation = generation, size = nrow(h) / 2L,
              heterozygosity = mean(2 * p * (1 - p)),
              segregating_sites = sum(p > 0 & p < 1))
        }
        burnRecords <- list(summarize(-burnin, haplotypes))
        burnStep <- max(1L, burnin %/% 20L)
        for (g in seq_len(burnin)) {
            haplotypes <- nextGeneration(haplotypes, sizes[1L])
            if (g %% burnStep == 0L || g == burnin)
                burnRecords[[length(burnRecords) + 1L]] <-
                    summarize(g - burnin, haplotypes)
        }
        ancestralFrequency <- colMeans(haplotypes)
        founderIndices <- integer()
        if (!is.null(founderSize)) {
            founders <- simupopSample(haplotypes, founderSize)
            haplotypes <- founders$haplotypes
            founderIndices <- founders$indices
        }
        p <- colMeans(haplotypes)
        generationRecords <- list(summarize(0L, haplotypes, p))
        frequencies <- list(p)
        for (g in seq_len(length(sizes) - 1L)) {
            haplotypes <- nextGeneration(haplotypes, sizes[g + 1L])
            if (g %% recordEvery == 0L || g == length(sizes) - 1L) {
                p <- colMeans(haplotypes)
                generationRecords[[length(generationRecords) + 1L]] <-
                    summarize(g, haplotypes, p)
                frequencies[[length(frequencies) + 1L]] <- p
            }
        }
        final <- simupopSample(haplotypes, sampleSize)
        generation <- as.data.frame(do.call(rbind, generationRecords))
        frequencyMatrix <- do.call(rbind, frequencies)
        chromosomeRecords <- do.call(rbind, lapply(seq_along(loci), function(j) {
            p <- frequencyMatrix[, markers$chromosome == j, drop = FALSE]
            data.frame(generation = generation$generation, chromosome = j,
                heterozygosity = rowMeans(2 * p * (1 - p)),
                segregating_sites = rowSums(p > 0 & p < 1))
        }))
        h <- final$haplotypes
        result <- list(observed = list(haplotypes = h,
                 genotypes = h[seq(1L, nrow(h), 2L), , drop = FALSE] +
                     h[seq(2L, nrow(h), 2L), , drop = FALSE],
                 sample_ids = sprintf("S%03d", seq_len(sampleSize)),
                 markers = markers, sequenceLength = sum(loci)),
             calibration = list(mutationRate = mutationRate,
                 reverseMutationRate = reverseMutationRate, markers = markers,
                 recombination_fractions = unlist(fractions, use.names = FALSE),
                 sequenceLength = sum(loci)),
             records = list(generation = generation,
                 allele_frequency = frequencyMatrix,
                 burnin = as.data.frame(do.call(rbind, burnRecords)),
                 ancestral_frequency = ancestralFrequency,
                 chromosome = chromosomeRecords),
             truth = list(sample_indices = final$indices,
                 founder_indices = founderIndices, burnin = burnin, seed = seed,
                 sizes = sizes, ancestralSize = sizes[1L],
                 reproductiveSizes = c(if (is.null(founderSize)) sizes[1L] else
                     founderSize, sizes[-1L])))
        if (!is.null(recordSites)) result$records$copyParents <-
            list(sites = recordSites, markers = markers[recordSites, , drop = FALSE],
                 parents = copyParents)
        result
    })
}

#' Finite-copy histories and reproducible sampling examples
#'
#' Each new copy samples one parent uniformly with replacement. copyHistory uses integer allele labels and gives every mutation a fresh label. Fixed-seed example functions restore the caller random-number state; ordinary unseeded copying consumes the caller stream.
#'
#' @name models_process_examples
#' @param N Positive integer diploid size; driftMoments accepts a vector of sizes.
#' @param generations Integer number of reproductive updates; positive for new-mutation paths, otherwise nonnegative.
#' @param cultures Positive integer number of independent populations.
#' @param p0 Scalar initial frequency in [0,1]; generation-record counts round 2N p0 to an integer.
#' @param seed Nonnegative integer random seed.
#' @param sizes Nonempty vector of positive integer diploid sizes, starting at generation zero.
#' @param initial Nonnegative integer allele labels for exactly 2 times sizes[1] initial copies.
#' @param mu Scalar probability of a fresh allele label per copied allele.
#' @param copies Positive integer haploid copy-pool size.
#' @param replicates Positive integer number of paths.
#' @param time Nonnegative integer vector of generation indices.
#' @returns copyHistory returns sizes, copies, zero-based parent and mutation positions, heterozygosity and allele_counts. samplingExample returns sampled indices and finite-population paths. populationProcessExample returns copying, mutation, ancestry and recombination histories. simulateFrequencyRecords returns counts and generation records. simulateNewMutationPaths returns a time-by-replicate matrix and zero-based introduction times. frequencyCountDistribution returns a count transition matrix and generation-by-count probabilities. driftMoments returns variance and diversity arrays indexed by time and N.
#' @examples
#' x <- simulateFrequencyRecords(N=8, generations=3, cultures=4)
#' x$records$frequency
#' driftMoments(c(10,100), 0:3)
#' set.seed(42)
#' copyHistory(c(2,2), 1:4)$copies
NULL

# Sources: Pritchard 2.1/2.3; Gorroochurn 8.2/8.8; Coop Ch.4 and public R examples.
# Label mutation and focal-interval tracking retain their explicit small models;
# uniform parental draws use the shared simuPOP kernel.

#' @rdname models_process_examples
#' @export
copyHistory <- function(sizes, initial, mu=0) {
    numericVector(sizes, "sizes", minimum = 1)
    if (any(sizes != floor(sizes))) stop("sizes must contain integer population sizes", call. = FALSE)
    numericVector(initial, "initial", 2 * sizes[1L], minimum = 0)
    if (any(initial != floor(initial))) stop("initial must contain integer allele labels", call. = FALSE)
    probability(mu, "mu")
    copies <- list(as.integer(initial))
    parents <- list()
    mutations <- list()
    next.id <- max(initial)+1L
    for (t in seq_len(length(sizes)-1L)) {
        old <- copies[[t]]
        draw <- simupopDrawParents(length(old), 2L*sizes[t+1L])
        new <- old[draw]
        changed <- which(runif(length(new)) < mu)
        if (length(changed)) {
            new[changed] <- seq.int(next.id, length.out=length(changed))
            next.id <- next.id+length(changed)
        }
        parents[[t]] <- as.integer(draw-1L) # zero-based visual positions
        mutations[[t]] <- as.integer(changed-1L)
        copies[[t+1L]] <- as.integer(new)
    }
    he <- vapply(copies, function(x) 1-sum((table(x)/length(x))^2), 0.0)
    list(sizes=as.integer(sizes), copies=copies, parents=parents,
         mutations=mutations, mu=mu, heterozygosity=he,
         allele_counts=vapply(copies, function(x) length(unique(x)), 0L))
}

#' @rdname models_process_examples
#' @export
samplingExample <- function() {
    restore <- rngRestorer()
    on.exit(restore(), add = TRUE)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(8149)
    draws <- sample.int(10L, 10L, replace=TRUE)
    first.count <- sum(draws <= 3L)
    paths <- lapply(c(5L, 50L), function(N) {
        k <- integer(21L)
        k[1L] <- round(2*N*0.3)
        for (t in 2:21) {
            k[t] <- if (N == 5L && t == 2L) first.count else
                rbinom(1L, size=2*N, prob=k[t-1L]/(2*N))
        }
        list(N=N, counts=k)
    })
    list(seed=8149L, p0=0.3, draws=draws, paths=paths, R_version=R.version.string)
}

#' @rdname models_process_examples
#' @export
populationProcessExample <- function() {
    restore <- rngRestorer()
    on.exit(restore(), add = TRUE)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    # Coop's single.crash: five large, two small, five large generations.
    set.seed(814901)
    bottleneck <- copyHistory(c(rep(10L,5L),rep(3L,2L),rep(10L,5L)), 1:20)
    
    # Coop Fig.4.7: N=5, locus mutation probability 0.1, initially one allele.
    set.seed(814902)
    mutation <- copyHistory(rep(5L,30L), rep(0L,10L), mu=0.1)
    
    # Coop Fig.4.2: distinct founder labels track descent without mutation.
    set.seed(814903)
    ancestry <- copyHistory(rep(5L,15L), 1:10)
    pair <- c(0L,1L) # fixed before looking at the history, not a selected pair
    trace <- list(pair)
    current <- pair
    meeting <- NULL
    for (g in seq.int(14L,1L)) {
        current <- ancestry$parents[[g]][current+1L]
        trace[[length(trace)+1L]] <- as.integer(current)
        if (current[1L]==current[2L]) {
            meeting <- g-1L
            break
        }
    }
    ancestry$sampled_pair <- pair
    ancestry$backward_trace <- trace
    ancestry$meeting_generation <- meeting
    
    # Follow the unbroken focal ancestral interval on two descendant paths.
    # Each gametic transmission has a Poisson(1) process on a 1-Morgan region.
    # Condition on inheritance of the ancestral copy at the central site.
    # Other homologs are unrelated to it; no mutation or gene conversion.
    set.seed(814904)
    focal <- 0.5
    intervals <- list(c(0,1),c(0,1))
    meioses <- c(0L,0L)
    recombination <- list(list(m=0L, counts=meioses, intervals=intervals,
                               shared=c(0,1), cuts=numeric(), path=NULL))
    for (m in 1:12) {
        path <- if (m %% 2L) 1L else 2L
        cuts <- sort(runif(rpois(1L,1)))
        old <- intervals[[path]]
        intervals[[path]] <- c(max(c(old[1L],cuts[cuts<focal])),
                              min(c(old[2L],cuts[cuts>focal])))
        meioses[path] <- meioses[path]+1L
        shared <- c(max(vapply(intervals, `[`, 0.0, 1L)),
                    min(vapply(intervals, `[`, 0.0, 2L)))
        recombination[[m+1L]] <- list(m=as.integer(m), counts=meioses,
            intervals=intervals, shared=shared, cuts=cuts, path=path-1L)
    }
    
    result <- list(seed_by_process=c(bottleneck=814901L, mutation=814902L,
        ancestry=814903L,recombination=814904L), R_version=R.version.string,
        bottleneck=bottleneck, mutation=mutation, ancestry=ancestry,
        recombination=recombination,
        sources=list(Pritchard="2.1 pp.76,80-83; 2.3 pp.106-107",
            Coop="Ch.4 pp.4.2,4.5,4.8; Figs.4.2,4.7,4.9",
            Coop_code=paste0("https://github.com/cooplab/popgen-notes/",
                             "blob/master/Rcode/Loss_of_heterozyg_varying_pop.R"),
            Notes="8.2,8.5,8.6,8.9",
            Handbook="16.2.3 pp.465-466"))
    result
}

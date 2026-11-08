# Upstream: Bo Peng and simuPOP contributors,
# GPL-3.0-or-later, https://github.com/BoPeng/simuPOP.
# Source functions:
# __init__.py::SNPMutator; mutator.cpp::BaseMutator::apply and
# MatrixMutator::{MatrixMutator,mutate}; mating.cpp::RandomParentChooser::
# chooseParents and CombinedParentsChooser::chooseParents; transmitter.cpp::
# Recombinator::{initialize,transmitGenotype,applyDuringMating};
# sampling.py::RandomSampler::drawSample and Population::extractIndividuals.

simupopCheckHaplotypes <- function(haplotypes) {
    haplotypes <- genotypeMatrix(haplotypes)
    if (nrow(haplotypes) %% 2L || any(haplotypes > 1L)) {
        stop("haplotypes must have paired diploid rows and binary 0/1 alleles",
             call. = FALSE)
    }
    unname(haplotypes)
}

# Draw precisely the success set of independent Bernoulli(probability) trials.
simupopBernoulli <- function(size, probability) {
    if (!size || probability == 0) return(integer())
    if (probability == 1) return(seq_len(size))
    count <- as.integer(stats::rbinom(1L, size, probability))
    if (!count) return(integer())
    if (count == size) return(seq_len(size))
    sort.int(sample.int(size, count, replace = FALSE,
                        useHash = count <= size / 2L), method = "quick")
}

# Equivalent to the two-state SNPMutator transition matrix. Mutation acts on
# parental chromosomes before they are used by any offspring, so two gametes
# inheriting the same parental entry share that entry's mutation outcome.
simupopMutate <- function(haplotypes, forward, reverse = forward,
                            events = NULL, returnEvents = FALSE) {
    if (!is.matrix(haplotypes) || !is.numeric(haplotypes) ||
        !nrow(haplotypes) || !ncol(haplotypes)) {
        stop("haplotypes must be a nonempty numeric matrix", call. = FALSE)
    }
    checkInputs(forward = list(positiveScalar, zero = TRUE),
                reverse = list(positiveScalar, zero = TRUE))
    if (forward > 1 || reverse > 1) {
        stop("mutation probabilities must lie between zero and one", call. = FALSE)
    }
    if (is.null(events)) {
        proposalRate <- max(forward, reverse)
        events <- simupopBernoulli(length(haplotypes), proposalRate)
        if (length(events) && forward != reverse) {
            old <- haplotypes[events]
            if (anyNA(old) || any(old != 0L & old != 1L)) {
                stop("SNPMutator supports only binary 0/1 alleles", call. = FALSE)
            }
            accept <- ifelse(old == 0L, forward, reverse) / proposalRate
            events <- events[stats::runif(length(events)) < accept]
        }
    } else if (!is.numeric(events) || any(!is.finite(events)) ||
               any(events != floor(events) | events < 1L |
                   events > length(haplotypes)) || anyDuplicated(events)) {
        stop("mutation events must be distinct one-based matrix indices",
             call. = FALSE)
    }
    if (length(events)) {
        old <- haplotypes[events]
        if (anyNA(old) || any(old != 0L & old != 1L)) {
            stop("SNPMutator supports only binary 0/1 alleles", call. = FALSE)
        }
        haplotypes[events] <- 1L - old
    }
    if (returnEvents) {
        return(list(haplotypes = haplotypes, events = as.integer(events)))
    }
    haplotypes
}

# Two independent ANY_SEX RandomParentChooser calls, with replacement and
# allowSelfing=TRUE. Each entry supplies one gamete of an offspring.
simupopParents <- function(parentSize, offspringSize, parents = NULL) {
    checkInputs(parentSize = scalarInteger, offspringSize = scalarInteger)
    parentSize <- as.integer(parentSize)
    offspringSize <- as.integer(offspringSize)
    if (offspringSize > .Machine$integer.max / 2L) {
        stop("offspringSize is too large", call. = FALSE)
    }
    if (is.null(parents)) {
        return(simupopDrawParents(parentSize, 2L * offspringSize))
    }
    if (!is.numeric(parents) || length(parents) != 2L * offspringSize ||
        any(!is.finite(parents)) || any(parents != floor(parents) |
            parents < 1L | parents > parentSize)) {
        stop("parents must index one parental diploid for every gamete",
             call. = FALSE)
    }
    as.integer(parents)
}

simupopSwitches <- function(gametes, probabilities) {
    empty <- matrix(integer(), nrow = 0L, ncol = 2L,
                    dimnames = list(NULL, c("gamete", "after")))
    if (!length(probabilities) || max(probabilities) == 0) return(empty)
    proposalRate <- max(probabilities)
    selected <- simupopBernoulli(gametes * length(probabilities), proposalRate)
    if (!length(selected)) return(empty)
    after <- (selected - 1L) %/% gametes + 1L
    if (any(probabilities != proposalRate)) {
        keep <- stats::runif(length(selected)) < probabilities[after] / proposalRate
        selected <- selected[keep]
        after <- after[keep]
    }
    result <- cbind(gamete = (selected - 1L) %% gametes + 1L, after = after)
    storage.mode(result) <- "integer"
    result
}

# Recombinator copies a homologue until the next successful interval switch.
# At each chromosome boundary its probability-one-half switch makes the next
# chromosome's starting homologue independent. Drawing each chromosome start
# independently therefore gives exactly the same joint distribution.
#
# Replay starts are 0/1, with rows=gametes and columns=chromosomes. Each switches
# matrix has columns gamete and after, where after is a chromosome-local site
# index: switching after site k changes the source beginning at site k+1.
simupopMeiosis <- function(haplotypes, parents, loci, fractions,
                            events = NULL, returnEvents = FALSE, binary = FALSE) {
    if (!is.matrix(haplotypes) || !is.numeric(haplotypes) ||
        !nrow(haplotypes) || nrow(haplotypes) %% 2L || !ncol(haplotypes)) {
        stop("haplotypes must have two rows for each parental diploid",
             call. = FALSE)
    }
    if (!is.numeric(loci) || !length(loci) || any(!is.finite(loci)) ||
        any(loci < 1L | loci != floor(loci)) || sum(loci) != ncol(haplotypes)) {
        stop("loci must partition the haplotype columns into chromosomes",
             call. = FALSE)
    }
    loci <- as.integer(loci)
    gametes <- length(parents)
    if (!is.numeric(parents) || !gametes || gametes %% 2L ||
        any(!is.finite(parents)) || any(parents != floor(parents) |
            parents < 1L | parents > nrow(haplotypes) / 2L)) {
        stop("parents must index a parental diploid for each paired offspring gamete",
             call. = FALSE)
    }
    parents <- as.integer(parents)
    if (!is.list(fractions) || length(fractions) != length(loci)) {
        stop("fractions must give interval-switch probabilities for each chromosome",
             call. = FALSE)
    }
    rates <- lapply(seq_along(loci), function(j) {
        rate <- fractions[[j]]
        if (!is.numeric(rate) || any(!is.finite(rate)) ||
            any(rate < 0 | rate > 0.5) ||
            !length(rate) %in% c(loci[j] - 1L, loci[j]) ||
            (length(rate) == loci[j] && tail(rate, 1L) != 0)) {
            stop("fractions must lie in [0,0.5], with an optional final zero",
                 call. = FALSE)
        }
        head(rate, loci[j] - 1L)
    })
    if (is.null(events)) {
        starts <- matrix(sample.int(2L, gametes * length(loci), replace = TRUE) - 1L,
                         nrow = gametes, ncol = length(loci))
        switches <- lapply(rates, function(rate) simupopSwitches(gametes, rate))
    } else {
        if (!is.list(events) || !is.matrix(events$starts) ||
            !identical(dim(events$starts), c(gametes, length(loci))) ||
            anyNA(events$starts) || any(!events$starts %in% 0:1) ||
            !is.list(events$switches) || length(events$switches) != length(loci)) {
            stop("meiosis events require a binary starts matrix and chromosome switch lists",
                 call. = FALSE)
        }
        starts <- events$starts
        switches <- events$switches
    }
    offspring <- matrix(0L, nrow = gametes, ncol = ncol(haplotypes))
    # All transmission events above use the full map. Constant binary columns
    # have the same offspring state for every parent choice and switch path.
    variable <- rep(TRUE, ncol(haplotypes))
    if (binary) {
        counts <- colSums(haplotypes)
        fixedOne <- which(counts == nrow(haplotypes))
        if (length(fixedOne)) offspring[, fixedOne] <- 1L
        variable <- counts > 0 & counts < nrow(haplotypes)
    }
    offset <- 0L
    for (j in seq_along(loci)) {
        columns <- seq.int(offset + 1L, offset + loci[j])
        active <- columns[variable[columns]]
        source <- 2L * parents - 1L + starts[, j]
        if (length(active)) offspring[, active] <- haplotypes[source, active, drop = FALSE]
        changes <- switches[[j]]
        if (!is.matrix(changes) || !is.numeric(changes) || ncol(changes) != 2L ||
            any(!is.finite(changes)) || any(changes != floor(changes)) ||
            (nrow(changes) && (any(changes[, 1L] < 1L | changes[, 1L] > gametes) ||
                any(changes[, 2L] < 1L | changes[, 2L] >= loci[j]) ||
                anyDuplicated(data.frame(changes))))) {
            stop("switch events require distinct gamete/after-site integer pairs",
                 call. = FALSE)
        }
        if (nrow(changes)) {
            groups <- split(changes[, 2L], changes[, 1L])
            for (id in names(groups)) {
                gamete <- as.integer(id)
                cuts <- sort.int(groups[[id]])
                other <- 2L * parents[gamete] - starts[gamete, j]
                ends <- c(cuts[-1L], loci[j])
                for (k in seq.int(1L, length(cuts), 2L)) {
                    replace <- active[active > offset + cuts[k] & active <= offset + ends[k]]
                    offspring[gamete, replace] <- haplotypes[other, replace]
                }
            }
        }
        offset <- offset + loci[j]
    }
    if (returnEvents) {
        return(list(haplotypes = offspring,
                    events = list(starts = starts, switches = switches)))
    }
    offspring
}

# RandomSampler draws a uniform subset of diploids without replacement.
# extractIndividuals retains their original population order, not shuffle order.
simupopSample <- function(haplotypes, size, indices = NULL) {
    if (!is.matrix(haplotypes) || !is.numeric(haplotypes) ||
        !nrow(haplotypes) || nrow(haplotypes) %% 2L || !ncol(haplotypes)) {
        stop("haplotypes must have paired diploid rows", call. = FALSE)
    }
    checkInputs(size = scalarInteger)
    size <- as.integer(size)
    populationSize <- nrow(haplotypes) %/% 2L
    if (size > populationSize) {
        stop("sample size cannot exceed the diploid population", call. = FALSE)
    }
    if (is.null(indices)) {
        indices <- sample.int(populationSize, size, replace = FALSE)
    } else if (!is.numeric(indices) || length(indices) != size ||
               any(!is.finite(indices)) || any(indices != floor(indices) |
                   indices < 1L | indices > populationSize) || anyDuplicated(indices)) {
        stop("indices must identify distinct diploids in the requested sample",
             call. = FALSE)
    }
    indices <- sort.int(as.integer(indices))
    rows <- as.vector(rbind(2L * indices - 1L, 2L * indices))
    list(haplotypes = unname(haplotypes[rows, , drop = FALSE]), indices = indices)
}

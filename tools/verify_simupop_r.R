#!/usr/bin/env Rscript
# Independent model checks for the self-contained simuPOP operator translation.
# R calculations verify the supported sampling and mutation models.
library(hpgR)

checks <- 0L
assert <- function(value, message) {
    if (!isTRUE(value)) stop(message, call. = FALSE)
    checks <<- checks + 1L
}

initial <- matrix(seq_len(60L), nrow = 10L)
for (seed in c(1L, 8149L, 1956L)) {
    x <- simulateForward(initial, c(5L, 2L, 9L, 3L), seed)
    assert(identical(x$states, replayCopying(initial, x$parents)),
        "Recorded parental indices did not replay exact clone states")
    assert(identical(x, simulateForward(initial, c(5L, 2L, 9L, 3L), seed)),
        "A repeated pure-R seed changed its result")
}

# Each gamete chooses one of the two homologues of an independently selected
# diploid parent. N=1 is therefore sufficient for a direct selfing test.
meiosis <- getFromNamespace(".simupopMeiosis", "hpgR")
parents <- getFromNamespace(".simupopParents", "hpgR")
mutate <- getFromNamespace(".simupopMutate", "hpgR")
set.seed(18301L)
p <- parents(1L, 20000L)
h <- meiosis(rbind(rep(0L, 4L), rep(1L, 4L)), p,
    c(2L, 2L), list(c(.2, 0), c(.35, 0)))
assert(all(p == 1L), "The N=1 parent pool failed to allow selfing")
for (j in seq_len(ncol(h)))
    assert(abs(mean(h[, j]) - .5) < .012, "Homologue starts are not uniform")
assert(abs(mean(h[, 1L] != h[, 2L]) - .2) < .012,
    "Adjacent-marker switching differs from its recombination fraction")
assert(abs(mean(h[, 3L] != h[, 4L]) - .35) < .012,
    "The second chromosome did not use its supplied recombination fraction")
assert(abs(mean(h[, 1L] != h[, 3L]) - .5) < .012,
    "Chromosome starts are not independent")
g <- h[seq(1L, nrow(h), 2L), 1L] + h[seq(2L, nrow(h), 2L), 1L]
assert(max(abs(tabulate(g + 1L, nbins = 3L) / length(g) - c(.25, .5, .25))) < .015,
    "The two gametes do not give the expected independent genotype proportions")

# Direct directional SNP mutation laws; no sequential reversal of a just-
# mutated allele is allowed within the same mutator application.
set.seed(18302L)
zeros <- matrix(0L, nrow = 400L, ncol = 200L)
assert(abs(mean(mutate(zeros, .15, .27)) - .15) < .007,
    "Forward mutation does not have the specified probability")
ones <- zeros + 1L
assert(abs(mean(mutate(ones, .15, .27)) - .73) < .007,
    "Reverse mutation does not have the specified probability")
assert(identical(mutate(ones, 1, 1), zeros),
    "A certain reverse mutation did not flip every ancestral one")
assert(identical(mutate(zeros, 1, 1), ones),
    "A certain forward mutation was reversed during the same step")

# Mutation takes place once on the finite parental copy pool, then offspring
# sample that pool. With symmetric mutation u, m=2N, a=1-1/m, s=1-2u:
# E[H_(t+1)] = a * {2*u*(1-u)*a + s^2*E[H_t]}.
# This tests timing as well as long-run mutation--drift balance. It is not the
# infinite-sites approximation H=4*N*u and does not require exact equilibrium.
N <- 2L
u <- .15
steps <- 120L
a <- 1 - 1 / (2 * N)
s <- 1 - 2 * u
fixed <- 2 * u * (1 - u) * a^2 / (1 - a * s^2)
expected <- fixed * (1 - (a * s^2)^steps)
means <- numeric(32L)
diversity <- numeric(32L)
for (r in seq_along(means)) {
    x <- simulatePopulation(N, rep(1L, 80L), u, N,
        burnin = steps, seed = 18400L + r)
    freq <- x$records$ancestral_frequency
    means[r] <- mean(freq)
    diversity[r] <- mean(2 * freq * (1 - freq))
}
expected_frequency <- (1 - s^steps) / 2
assert(abs(mean(means) - expected_frequency) < .04,
    "Symmetric mutation's expected allele frequency was not recovered")
assert(abs(mean(diversity) - expected) < .03,
    "Finite-parental-mutation diversity does not match its exact recurrence")

# Asymmetric mutation has expectation p*=u/(u+v), not necessarily 1/2.
asymmetric <- vapply(1:8, function(r) {
    x <- simulatePopulation(12L, rep(1L, 80L), .03, 12L,
        reverseMutationRate = .09, burnin = 120L, seed = 18500L + r)
    mean(x$records$ancestral_frequency)
}, numeric(1))
assert(abs(mean(asymmetric) - .25 * (1 - .88^120)) < .04,
    "Asymmetric mutation's expected frequency was not recovered")

cat("PASS:", checks, "self-contained operator and probability-model checks.\n")
cat(sprintf("Symmetric mutation: mean frequency %.5f (expected %.5f); mean H %.5f (expected %.5f).\n",
    mean(means), expected_frequency, mean(diversity), expected))
cat(sprintf("Asymmetric mutation: mean frequency %.5f (expected %.5f).\n",
    mean(asymmetric), .25 * (1 - .88^120)))
cat("Original simuPOP streams and independent ensembles are checked by tools/verify_simupop_reference.R.\n")


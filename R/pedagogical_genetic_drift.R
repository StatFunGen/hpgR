# Finite-population allele counts, mutation paths and analytical moments.

#' @rdname models_process_examples
#' @export
simulateFrequencyRecords <- function(N=16L, generations=19L, cultures=105L, p0=.5, seed=1956L) {
    scalarInteger(N, "N")
    scalarInteger(generations, "generations", 0L)
    scalarInteger(cultures, "cultures")
    scalarInteger(seed, "seed", 0L)
    probability(p0, "p0")
    restore <- rngRestorer()
    on.exit(restore(), add = TRUE)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(seed)
    record_N <- as.integer(N)
    record_generations <- as.integer(generations)
    record_cultures <- as.integer(cultures)
    generation_counts <- matrix(0L, record_generations + 1L, record_cultures)
    generation_counts[1, ] <- as.integer(round(2L * record_N * p0))
    for (t in seq_len(record_generations)) {
        generation_counts[t + 1L, ] <- rbinom(record_cultures, 2L * record_N,
                                            generation_counts[t, ] / (2L * record_N))
    }
    generation_records <- list(
        generation = 0:record_generations,
        frequency = generation_counts / (2L * record_N),
        generating_N = record_N
    )
    list(counts=generation_counts, records=generation_records)
}

#' @rdname models_process_examples
#' @export
simulateNewMutationPaths <- function(copies=100L, replicates=200L, generations=1000L, seed=8202L) {
    scalarInteger(copies, "copies")
    scalarInteger(replicates, "replicates")
    scalarInteger(generations, "generations")
    scalarInteger(seed, "seed", 0L)
    restore <- rngRestorer()
    on.exit(restore(), add = TRUE)
    set.seed(seed)
    mutation_copies <- as.integer(copies)
    mutation_replicates <- as.integer(replicates)
    introduction <- sample(0:(generations - 1L), mutation_replicates, replace = TRUE)
    mutation_paths <- matrix(NA_real_, (generations + 1L), mutation_replicates)
    mutation_paths[cbind(introduction + 1L, seq_len(mutation_replicates))] <- 1 / mutation_copies
    for (t in seq_len(generations)) {
        active <- which(!is.na(mutation_paths[t, ]) & mutation_paths[t, ] > 0 &
                        mutation_paths[t, ] < 1)
        mutation_paths[t + 1L, active] <- rbinom(length(active), mutation_copies,
                                               mutation_paths[t, active]) / mutation_copies
        absorbed <- which(mutation_paths[t, ] %in% c(0, 1))
        mutation_paths[t + 1L, absorbed] <- mutation_paths[t, absorbed]
    }
    list(paths=mutation_paths, introduction=introduction)
}

#' @rdname models_process_examples
#' @export
frequencyCountDistribution <- function(N=16L, generations=19L) {
    scalarInteger(N, "N")
    scalarInteger(generations, "generations", 0L)
    record_N <- as.integer(N)
    record_generations <- as.integer(generations)
    M <- 2L * record_N
    transition <- t(vapply(0:M, function(k) dbinom(0:M, M, k / M), numeric(M + 1L)))
    model_probability <- matrix(0, record_generations + 1L, M + 1L)
    model_probability[1, record_N + 1L] <- 1
    for (t in seq_len(record_generations)) {
        model_probability[t + 1L, ] <- model_probability[t, ] %*% transition
    }
    list(transition=transition, probability=model_probability)
}

#' @rdname models_process_examples
#' @export
driftMoments <- function(N, time, p0=.5) {
    numericVector(N, "N", minimum = 1)
    if (any(N != floor(N))) stop("N must contain integer population sizes", call. = FALSE)
    numericVector(time, "time", minimum = 0)
    if (any(time != floor(time))) stop("time must contain integer generations", call. = FALSE)
    probability(p0, "p0")
    list(variance=sapply(N, function(n) p0*(1-p0)*(1-(1-1/(2*n))^time)),
         diversity=sapply(N, function(n) 2*p0*(1-p0)*(1-1/(2*n))^time))
}

#' @rdname models_demographic_processes
#' @export
zeroCountUpperBound <- function(n, confidence = .95) {
    numericVector(n, "n", minimum = 1)
    if (any(n != floor(n))) stop("n must contain positive integers", call. = FALSE)
    probability(confidence, "confidence")
    if (confidence == 0 || confidence == 1)
        stop("confidence must be strictly between 0 and 1", call. = FALSE)
    -expm1(log1p(-confidence) / n)
}

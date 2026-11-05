# Developer reference checks only. These functions do not set the package RNG.
# Source translations from the GSL code bundled with simuPOP (GPL >= 3).
# gsl/rng/mt.c mt_set, copyright Makoto Matsumoto, Takuji Nishimura and Brian Gough.
# gsl/gsl_rng.h gsl_rng_uniform_int; gsl/randist/geometric.c gsl_ran_geometric,
# copyright James Theiler and Brian Gough. Exact hashes are in the fixture.

simupopTestSeed <- function(seed) {
    stopifnot(seed > 0, seed <= 2^32 - 1, seed == floor(seed))
    key <- numeric(624L)
    key[1L] <- seed
    multiplier <- 1812433253
    for (i in 2:624) {
        old <- key[i - 1L]
        shifted <- bitwXor(as.integer(old %% 65536),
            as.integer(floor(old / 2^30))) + 65536 * floor(old / 65536)
        # Split the product so every integer operation is exact in a double.
        key[i] <- (multiplier * (shifted %% 65536) +
            (multiplier %% 65536) * floor(shifted / 65536) * 65536 + i - 1) %% 2^32
    }
    signed <- suppressWarnings(as.integer(ifelse(key > 2^31 - 1, key - 2^32, key)))
    assign(".Random.seed", c(10403L, 624L, signed), envir = .GlobalEnv)
}

simupopTestUniform <- function(n) {
    # Recover the raw 32-bit word, including R's special handling of word zero.
    floor(stats::runif(n) * 2^32) / 2^32
}

simupopTestInteger <- function(bound, n = 1L) {
    scale <- floor((2^32 - 1) / bound)
    vapply(seq_len(n), function(i) {
        repeat {
            value <- floor(floor(stats::runif(1L) * 2^32) / scale)
            if (value < bound) return(value)
        }
    }, 0)
}

simupopTestGeometric <- function(n, p) {
    vapply(seq_len(n), function(i) {
        repeat {
            u <- simupopTestUniform(1L)
            if (u > 0) break
        }
        if (p == 1) 1 else floor(log(u) / log(1 - p) + 1)
    }, 0)
}

simupopTestParents <- function(seed, sizes) {
    simupopTestSeed(seed)
    bitIndex <- 0L
    parents <- lapply(seq_len(length(sizes) - 1L), function(t) {
        out <- integer(2L * sizes[t + 1L])
        for (i in seq_along(out)) {
            out[i] <- as.integer(simupopTestInteger(2L * sizes[t])) + 1L
            # Native RandomSelection consumes a sex draw before cloning.
            # RNG::randBit refills a 16-bit cache, even though cloning replaces sex.
            if (bitIndex == 16L) bitIndex <<- 0L
            if (bitIndex == 0L) invisible(simupopTestInteger(65535))
            bitIndex <<- bitIndex + 1L
        }
        out
    })
    parents
}

simupopTestRng <- function(fixture) {
    getFromNamespace(".simupopWithSeed", "hpgR")(1L, function() {
        rows <- lapply(fixture$results, function(d) {
            simupopTestSeed(d$seed)
            uniforms <- simupopTestUniform(length(d$uniform))
            integerEqual <- vapply(names(d$integers), function(bound) {
                simupopTestSeed(d$seed)
                identical(simupopTestInteger(as.numeric(bound), length(d$integers[[bound]])),
                    as.numeric(unlist(d$integers[[bound]])))
            }, logical(1))
            simupopTestSeed(d$seed)
            geometric <- simupopTestGeometric(length(d$geometric), .15)
            simupopTestSeed(d$seed)
            binomial <- stats::rbinom(length(d$binomial), 20, .3)
            simupopTestSeed(d$seed)
            rGeometric <- stats::rgeom(length(d$geometric), .15) + 1
            simupopTestSeed(d$seed)
            rPoisson <- stats::rpois(length(d$poisson), 2.3)
            parents <- simupopTestParents(d$seed, c(5L, 2L, 7L, 3L))
            referenceParents <- lapply(d$clone_parents, function(x) as.integer(unlist(x)) + 1L)
            states <- hpgR::replayCopying(cbind(0:9, (0:9) %% 2L), parents)[-1L]
            stateEqual <- vapply(seq_along(states), function(j)
                identical(as.integer(t(states[[j]])),
                    as.integer(unlist(d$clone_states[[j]]))), logical(1))
            set.seed(d$seed)
            data.frame(seed = d$seed,
                uniform_equal = identical(uniforms, as.numeric(unlist(d$uniform))),
                integer_equal = all(integerEqual),
                geometric_equal = identical(geometric, as.numeric(unlist(d$geometric))),
                copy_equal = identical(parents, referenceParents) && all(stateEqual),
                default_uniform_equal = identical(stats::runif(length(d$uniform)),
                    as.numeric(unlist(d$uniform))),
                r_binomial_equal = identical(as.numeric(binomial), as.numeric(unlist(d$binomial))),
                r_poisson_equal = identical(rPoisson, as.numeric(unlist(d$poisson))),
                r_geometric_equal = identical(rGeometric, as.numeric(unlist(d$geometric))))
        })
        do.call(rbind, rows)
    })
}

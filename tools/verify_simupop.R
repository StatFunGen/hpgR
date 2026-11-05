#!/usr/bin/env Rscript
library(hpgR)
checks <- 0L
for (N in c(2L, 5L, 16L)) for (seed in c(1L, 8149L, 1956L)) {
    initial <- matrix(seq_len(2L * N) %% 2L, ncol = 1L)
    sizes <- c(N, N, 2L, N)
    original <- simulateForward(initial, sizes, seed)
    replay <- replayCopying(initial, original$parents)
    stopifnot(identical(replay, original$states))
    repeated <- simulateForward(initial, sizes, seed)
    stopifnot(identical(original, repeated))
    checks <- checks + 1L
}
cat("PASS:", checks, "R histories replayed exactly from recorded parents;",
    "same-seed R repetition passed.\n")

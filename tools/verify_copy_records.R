# Run in a clean R process. Optional --report PATH writes the machine-readable record.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
course <- Sys.getenv("HPGR_COURSE_DIR")
library(hpgR)
checks <- list()
check <- function(name, body) {
    timing <- system.time(answer <- tryCatch({body(); NULL}, error = identity))
    checks[[name]] <<- list(passed = is.null(answer), seconds = unname(timing["elapsed"]),
                           error = if (is.null(answer)) NULL else conditionMessage(answer))
    cat(if (is.null(answer)) "PASS" else "FAIL", name, "\n")
}
fails <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity),
                                         "error"))
core <- pgRecordAncestry(rep(5L, 5L), 1000, 6L, ancestry_seed = 31L)

check("record/replay/reproducibility", function() {
    pgValidateAncestry(core)
    again <- pgRecordAncestry(rep(5L, 5L), 1000, 6L, ancestry_seed = 31L)
    stopifnot(identical(core, again), length(core$local_trees) == 1L)
    mutation <- pgMutateAncestry(core, 0.0001, 71L)
    stopifnot(identical(mutation$haplotypes, pgReplay(core, mutation$mutations)))
    changed <- core
    changed$local_trees <- NULL
    stopifnot(identical(mutation$haplotypes, pgReplay(changed, mutation$mutations)))
})
check("seed is local; mutation leaves ancestry fixed", function() {
    set.seed(135)
    before <- .Random.seed
    a <- pgSimulate(rep(3L, 4L), 100, 4L, 0, ancestry_seed = 9L)
    stopifnot(identical(before, .Random.seed))
    b <- pgSimulate(rep(3L, 4L), 100, 4L, .01, ancestry_seed = 9L)
    stopifnot(identical(a$truth$regions[[1]]$ancestry, b$truth$regions[[1]]$ancestry))
    seeds <- unlist(lapply(1:20, function(base)
        vapply(1:3, function(region) pgRegionSeed(base, region), integer(1))))
    stopifnot(!anyDuplicated(seeds))
})
check("empty mutation matrices; smallest population; invalid inputs", function() {
    d <- pgSimulate(c(1L, 1L), 100, 2L, 0, ancestry_seed = 41L)
    s <- pgObservedSummaries(d$observed, d$calibration)
    stopifnot(identical(dim(d$observed$regions[[1]]$haplotypes), c(2L, 0L)),
              s$regions$S == 0, s$regions$pi_site == 0, sum(s$sfs[[1]]) == 0)
    fails(pgRecordAncestry(c(0, 3), 100, 2))
    fails(pgRecordAncestry(c(2.5, 3), 100, 2))
    fails(pgRecordAncestry(c(2, 3), 100, 7))
    fails(pgRecordAncestry(c(2, 3), 100, 2, map_length = -1))
    fails(pgRecordAncestry(c(2, 3), 100, 2, sample_indices = c(1, 1)))
    fails(pgMutateAncestry(core, -1))
})
check("root completion preserves sample; truncated forests fail loudly", function() {
    forest <- pgRecordAncestry(5L, 100, 4L, ancestry_seed = 51L,
                                 complete_roots = FALSE)
    pgValidateAncestry(forest, require_complete = FALSE)
    fails(pgValidateAncestry(forest))
    fails(pgMutateAncestry(forest, .001))
    fails(pgRecordAncestry(5L, 100, 4L, ancestry_seed = 51L, max_extension = 0L))
    whole <- pgRecordAncestry(5L, 100, 4L, ancestry_seed = 51L)
    stopifnot(whole$root_complete, whole$extension_generations > 0,
              identical(whole$sample_nodes, forest$sample_nodes),
              whole$sample_generation == 0, max(whole$nodes$generation) == 0)
    fixed <- pgRecordAncestry(rep(5L, 5L), 1000, 6L, ancestry_seed = 31L,
                                complete_roots = FALSE)
    forward_edges <- core$edges[core$nodes$generation[core$edges$child] > 0, ]
    rownames(forward_edges) <- NULL
    stopifnot(identical(forward_edges, fixed$edges),
              identical(core$sample_nodes, fixed$sample_nodes))
})
check("linked intervals and replay", function() {
    linked <- pgRecordAncestry(rep(5L, 8L), 1000, 6L, map_length = .08,
                                 ancestry_seed = 63L)
    pgValidateAncestry(linked)
    mutation <- pgMutateAncestry(linked, .0002, 74L)
    stopifnot(identical(mutation$haplotypes, pgReplay(linked, mutation$mutations)),
              all(vapply(linked$local_trees, `[[`, logical(1), "complete")))
    bounds <- sort(unique(c(0, linked$edges$left, linked$edges$right, 1000)))
    # Independently check one root at every raw-edge-defined interval, even if
    # no mutation and no retained-tree boundary occurs there.
    for (x in (head(bounds, -1L) + tail(bounds, -1L)) / 2) {
        keep <- linked$edges$left <= x & x < linked$edges$right
        parent <- integer(nrow(linked$nodes))
        parent[linked$edges$child[keep]] <- linked$edges$parent[keep]
        lineage <- linked$sample_nodes
        while (length(unique(lineage)) > 1L && all(parent[lineage] > 0L))
            lineage <- parent[lineage]
        stopifnot(length(unique(lineage)) == 1L)
    }
})
check("malformed edges/times/mutations rejected", function() {
    bad <- core
    bad$edges$parent[1] <- nrow(core$nodes) + 1L
    fails(pgValidateAncestry(bad))
    bad <- core
    bad$edges$parent[1] <- bad$edges$child[1]
    fails(pgValidateAncestry(bad))
    bad <- core
    bad$edges$left[1] <- .01
    fails(pgValidateAncestry(bad))
    bad <- core
    bad$edges <- rbind(bad$edges, bad$edges[1, ])
    fails(pgValidateAncestry(bad))
    mutation <- pgMutateAncestry(core, .001, 80L)$mutations
    stopifnot(nrow(mutation) > 0)
    mutation$time[1] <- -1
    fails(pgReplay(core, mutation))
})
check("observed-only denominator and SFS", function() {
    G <- cbind(c(1L, 0L, 0L, 0L), c(1L, 1L, 0L, 0L))
    rownames(G) <- paste0("S", 1:4)
    observed <- list(sample_ids = rownames(G), sampling_time = 0,
        regions = list(list(region_id = "R1", analyzed_length_bp = 100,
                            positions_bp = c(2, 7), haplotypes = G)))
    calibration <- list(mu_per_base_per_generation = 1e-5)
    s <- pgObservedSummaries(observed, calibration)
    stopifnot(s$regions$S == 2, abs(s$regions$pi_region - 7/6) < 1e-12,
              abs(s$regions$pi_site - 7/600) < 1e-12,
              identical(s$sfs[[1]], c(1L, 1L, 0L)))
    # Round-trip just observations/calibration into a clean environment. No truth file.
    out <- tempfile(fileext = ".rds")
    saveRDS(list(observed = observed, calibration = calibration), out)
    env <- new.env(parent = baseenv())
    env$pgObservedSummaries <- hpgR::pgObservedSummaries
    env$public <- readRDS(out)
    result <- evalq(pgObservedSummaries(public$observed, public$calibration), env)
    stopifnot(identical(result, s))
    observed$regions <- list()
    fails(pgObservedSummaries(observed, calibration))
    unlink(out)
})
check("one-generation mean variance and parental collision", function() {
    pgWithSeed(95L, function() {
        J <- 30000L
        parent <- pgParentDraws(1:40, 20L * J)
        p_next <- colMeans(matrix(parent <= 12L, 20L, J))
        stopifnot(abs(mean(p_next) - .3) < 6 * sqrt(.3 * .7 / 20 / J),
                  abs(var(p_next) / (.3 * .7 / 20) - 1) < .06)
        first <- pgParentDraws(1:40, J)
        second <- pgParentDraws(1:40, J)
        stopifnot(abs(mean(first == second) - 1/40) < 6 * sqrt(.025 * .975 / J))
        parent4 <- matrix(pgParentDraws(1:10, 4L * J), 4L, J)
        all_distinct <- colSums(parent4[c(1, 1, 1, 2, 2, 3), ] ==
                                parent4[c(2, 3, 4, 3, 4, 4), ]) == 0
        expected <- prod((10 - 0:3) / 10)
        stopifnot(abs(mean(all_distinct) - expected) < 6 * sqrt(expected / J))
    })
})
check("conditional mutation count mean and variance", function() {
    means <- vapply(1:160, function(seed)
        ncol(pgMutateAncestry(core, .0001, seed)$haplotypes), integer(1))
    edges <- core$local_trees[[1]]$edges
    lambda <- sum((edges$right - edges$left) *
                   (core$nodes$time[edges$parent] - core$nodes$time[edges$child])) * .0001
    stopifnot(abs(mean(means) - lambda) < 6 * sqrt(lambda / length(means)),
              abs(var(means) - lambda) < 6 * lambda * sqrt(2 / (length(means) - 1)))
})
check("geometric waiting and discrete multifurcations", function() {
    times <- vapply(1:300, function(seed) {
        record <- pgRecordAncestry(4L, 10, 2L, ancestry_seed = seed)
        record$local_trees[[1]]$root_age
    }, numeric(1))
    stopifnot(abs(mean(times) - 8) < 6 * sqrt(8 * 7 / length(times)),
              abs(mean(times > 8) - (7/8)^8) < 6 / sqrt(length(times)))
    fourway <- pgRecordAncestry(2L, 10, 4L, ancestry_seed = 38L)
    tree <- fourway$local_trees[[1]]
    stopifnot(tree$root_age == 1,
              sum(tree$edges$parent == tree$root) == 4L)
    pgWithSeed(5L, function() {
        only_donor <- pgInheritance(77L, 101:110, 100, .5)
        stopifnot(all(only_donor$parent == 77L))
    })
})

performance <- list()
for (seed in c(101L, 202L, 303L)) {
    check(paste0("multi-seed benchmark ", seed), function() {
        gc(reset = TRUE)
        elapsed <- system.time(result <- pgSimulate(
            c(40L, rep(10L, 10L), rep(40L, 10L)), 1000, n = 12L,
            ancestry_seed = seed, mutation_seed = seed + 1000L))["elapsed"]
        region <- result$truth$regions[[1]]
        pgCheck(identical(unname(result$observed$regions[[1]]$haplotypes),
                           pgReplay(region$ancestry, region$mutations)), "benchmark replay")
        performance[[as.character(seed)]] <<- list(
            elapsed_seconds = unname(elapsed), object_bytes = as.numeric(object.size(result)),
            R_gc_heap_high_water_bound_MiB = sum(gc()[, ncol(gc())]),
            S = ncol(result$observed$regions[[1]]$haplotypes),
            extension = region$ancestry$extension_generations)
    })
}
args <- commandArgs(trailingOnly = TRUE)
if ("--report" %in% args) {
    at <- match("--report", args)
    stopifnot(at < length(args), requireNamespace("jsonlite", quietly = TRUE))
    jsonlite::write_json(list(R = R.version.string, checks = checks,
                             performance = performance), args[at + 1L],
                         pretty = TRUE, auto_unbox = TRUE)
}
if (!all(vapply(checks, `[[`, logical(1), "passed"))) stop("P1 checks failed")
cat(length(checks), "checks passed. No inference recovery or site deployment tested.\n")

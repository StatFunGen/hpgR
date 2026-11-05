msprimeAncestry <- utils::getFromNamespace("msprimeAncestry", "hpgR")

test_that("the restricted Hudson translation gives a complete zero-recombination tree", {
    x <- msprimeAncestry(2L, 100, 10L, 0, 81L, FALSE)
    expect_s3_class(x, "hpgR_ancestry")
    expect_identical(x$sampleIds, 0:3)
    expect_equal(nrow(x$nodes), 7L)
    expect_equal(nrow(x$edges), 6L)
    expect_true(all(x$edges$left == 0 & x$edges$right == 10))
    expect_equal(unname(x$eventCounts["recombination"]), 0)
    expect_equal(unname(x$eventCounts["common_ancestor"]), 3)
    expect_true(all(x$nodes$time[x$edges$parent + 1L] >
                    x$nodes$time[x$edges$child + 1L]))
})

test_that("Hudson event replay reproduces exact tables without sharing the RNG", {
    x <- msprimeAncestry(2L, 100, 100L, 1e-4, 8149L, TRUE)
    y <- msprimeAncestry(2L, 100, 100L, 1e-4, 1L, TRUE, events = x$events)
    expect_identical(x$nodes, y$nodes)
    expect_identical(x$edges, y$edges)
    expect_identical(x$eventCounts, y$eventCounts)
    expect_equal(sum(x$nodes$flags == 2^17),
                 2 * unname(x$eventCounts["recombination"]))
    expect_equal(sum(x$nodes$flags == 2^18),
                 unname(x$eventCounts["nonoverlap_common_ancestor"]))
    expect_gt(unname(x$eventCounts["nonoverlap_common_ancestor"]), 0)
    expect_gt(unname(x$eventCounts["gap_recombination"]), 0)
    expect_true(all(x$edges$left == floor(x$edges$left)))
    expect_true(all(x$edges$right == floor(x$edges$right)))
})

test_that("nonoverlapping common ancestors and full ARG unary edges are retained", {
    events <- list(
        list(type = "recombination", time = 1, lineage = 0L, breakpoint = 4),
        list(type = "common_ancestor", time = 2, lineages = c(0L, 2L)),
        list(type = "recombination", time = 3, lineage = 3L, breakpoint = 4),
        list(type = "common_ancestor", time = 4, lineages = c(1L, 3L)),
        list(type = "common_ancestor", time = 5, lineages = c(4L, 5L)))
    x <- msprimeAncestry(1L, 100, 10L, .001, 1L, TRUE, events = events)
    expect_identical(x$nodes$id, 0:8)
    expect_equal(x$nodes$time, c(0, 0, 1, 1, 2, 3, 3, 4, 5))
    expect_equal(x$nodes$flags, c(1, 1, 2^17, 2^17, 2^18, 2^17, 2^17, 0, 0))
    expected <- data.frame(left = c(0, 4, 0, 4, 0, 4, 0, 0, 4, 4),
        right = c(4, 10, 4, 10, 4, 10, 10, 4, 10, 10),
        parent = as.integer(c(2, 3, 4, 4, 5, 6, 7, 7, 8, 8)),
        child = as.integer(c(0, 0, 2, 3, 4, 4, 1, 5, 6, 7)))
    expect_identical(x$edges, expected)
    expect_equal(unname(x$eventCounts["nonoverlap_common_ancestor"]), 1)
})

test_that("an original msprime full-ARG fixture replays to exact native tables", {
    # Native reference: msprime 1.4.1, sim_ancestry(samples=1,ploidy=2,
    # population_size=5,sequence_length=6,recombination_rate=.02,
    # model="hudson",record_full_arg=TRUE,random_seed=360).
    # The full native tables expose each event's selected ancestral nodes.
    # Its gap recombination has no edge covering (3,5); cutting at 3 gives
    # exactly the same two segment chains as any native integer cut in that gap.
    t <- c(3.4435670429738514, 3.8122050292099168, 4.792741709693215,
           4.906233402598424, 5.296546714101565, 5.515541372702208,
           10.584786549549293)
    events <- list(
        list(type = "recombination", time = t[1L], lineage = 0L, breakpoint = 5),
        list(type = "recombination", time = t[2L], lineage = 1L, breakpoint = 3),
        list(type = "common_ancestor", time = t[3L], lineages = c(0L, 3L)),
        list(type = "recombination", time = t[4L], lineage = 4L, breakpoint = 3),
        list(type = "common_ancestor", time = t[5L], lineages = c(1L, 5L)),
        list(type = "common_ancestor", time = t[6L], lineages = c(4L, 6L)),
        list(type = "common_ancestor", time = t[7L], lineages = c(2L, 7L)))
    x <- msprimeAncestry(1L, 5, 6L, .02, 1L, TRUE, events = events)
    nodes <- data.frame(id = 0:11,
        time = c(0, 0, t[1L], t[1L], t[2L], t[2L], t[3L],
                 t[4L], t[4L], t[5L], t[6L], t[7L]),
        flags = as.integer(c(1, 1, 2^17, 2^17, 2^17, 2^17, 0,
                             2^17, 2^17, 2^18, 0, 0)))
    edges <- data.frame(
        left = c(0, 5, 0, 3, 0, 3, 0, 5, 0, 5, 0, 0, 5, 5, 5),
        right = c(5, 6, 3, 6, 5, 6, 3, 6, 3, 6, 3, 3, 6, 6, 6),
        parent = as.integer(c(2, 3, 4, 5, 6, 6, 7, 8, 9, 9, 10, 10, 10, 11, 11)),
        child = as.integer(c(0, 0, 1, 1, 2, 5, 6, 6, 4, 8, 7, 9, 9, 3, 10)))
    expect_identical(x$nodes, nodes)
    expect_identical(x$edges, edges)
    expect_equal(unname(x$eventCounts["gap_recombination"]), 1)
    expect_equal(unname(x$eventCounts["nonoverlap_common_ancestor"]), 1)
})

test_that("full ARG recording leaves marginal root ages unchanged", {
    full <- msprimeAncestry(2L, 100, 100L, 1e-4, 8149L, TRUE)
    small <- msprimeAncestry(2L, 100, 100L, 1e-4, 1L, FALSE,
                             events = full$events)
    rootTime <- function(x, position) {
        e <- x$edges[x$edges$left <= position & x$edges$right > position, ]
        parent <- rep(NA_integer_, nrow(x$nodes))
        parent[e$child + 1L] <- e$parent
        roots <- vapply(x$sampleIds, function(u) {
            while (!is.na(parent[u + 1L])) u <- parent[u + 1L]
            u
        }, integer(1L))
        expect_equal(length(unique(roots)), 1L)
        x$nodes$time[roots[1L] + 1L]
    }
    breaks <- sort(unique(c(full$edges$left, full$edges$right,
                            small$edges$left, small$edges$right)))
    positions <- head(breaks, -1L) + diff(breaks) / 2
    expect_equal(vapply(positions, function(p) rootTime(full, p), numeric(1L)),
                 vapply(positions, function(p) rootTime(small, p), numeric(1L)))
})

test_that("the two-copy Hudson clock has mean 2Ne", {
    ages <- vapply(seq_len(1000L), function(seed)
        max(msprimeAncestry(1L, 100, 1L, 0, seed, FALSE)$nodes$time), numeric(1L))
    expect_lt(abs(mean(ages) - 200), 200 * 5 / sqrt(length(ages)))
    expect_lt(abs(mean(ages > 200) - exp(-1)), .06)
})

test_that("Hudson simulation preserves the caller RNG and rejects unsupported inputs", {
    set.seed(49L)
    before <- .Random.seed
    msprimeAncestry(1L, 100, 1L, 0, 2L, FALSE)
    expect_identical(.Random.seed, before)
    expect_error(msprimeAncestry(1L, 100, 1.5, 0, 2L, FALSE), "sequenceLength")
    expect_error(msprimeAncestry(1L, 100, 10L, -1, 2L, FALSE), "recombinationRate")
    expect_error(msprimeAncestry(1L, 100, 10L, 0, 2L, FALSE,
        events = list()), "ended before")
    expect_error(msprimeAncestry(1L, 100, 10L, 0, 2L, FALSE,
        events = list(list(type = "recombination", time = 1,
                           lineage = 0L, breakpoint = 4))), "positive rate")
})

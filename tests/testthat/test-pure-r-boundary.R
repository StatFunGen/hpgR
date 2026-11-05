boundaryCallViolations <- function(code) {
    found <- character()
    native <- c(".Call", ".External", ".External2", ".C", ".Fortran",
        "dyn.load", "library.dynam", "system", "system2", "shell", "pipe")
    external <- c("reticulate", "Rcpp", "rPython", "basilisk")
    visit <- function(x) {
        if (!is.call(x)) {
            if (is.pairlist(x) || is.expression(x) || is.list(x))
                invisible(lapply(as.list(x), visit))
            return(invisible(NULL))
        }
        h <- x[[1L]]
        if (is.symbol(h) && as.character(h) %in% native)
            found <<- c(found, paste(deparse(x, width.cutoff = 120L), collapse = " "))
        if (is.symbol(h) && as.character(h) %in% c("::", ":::") &&
            length(x) >= 3L && (as.character(x[[2L]]) %in% external ||
                as.character(x[[3L]]) %in% native))
            found <<- c(found, paste(deparse(x, width.cutoff = 120L), collapse = " "))
        if (is.symbol(h) && as.character(h) %in% c("library", "require",
            "requireNamespace", "getNamespace", "asNamespace", "getExportedValue") &&
            length(x) >= 2L && as.character(x[[2L]])[1L] %in% external)
            found <<- c(found, paste(deparse(x, width.cutoff = 120L), collapse = " "))
        if (is.symbol(h) && as.character(h) == "do.call" && length(x) >= 2L &&
            is.character(x[[2L]]) && x[[2L]] %in% native)
            found <<- c(found, paste(deparse(x, width.cutoff = 120L), collapse = " "))
        invisible(lapply(as.list(x), visit))
    }
    visit(code)
    unique(found)
}

test_that("the executable namespace contains only self-contained R calculations", {
    ns <- asNamespace("hpgR")
    violations <- character()
    embedded_python <- character()
    for (name in ls(ns, all.names = TRUE)) {
        object <- get(name, envir = ns, inherits = FALSE)
        if (is.function(object)) {
            calls <- unique(c(boundaryCallViolations(body(object)),
                boundaryCallViolations(formals(object))))
            if (length(calls)) violations <- c(violations, paste(name, calls, sep = ": "))
        } else if (is.character(object) && length(object) &&
            any(grepl("(^|\n)[[:space:]]*(import[[:space:]]+[[:alnum:]_]+|from[[:space:]]+[[:alnum:]_.]+[[:space:]]+import|def[[:space:]]+[[:alnum:]_]+[[:space:]]*\\()",
                object))) {
            embedded_python <- c(embedded_python, name)
        }
    }
    expect_identical(violations, character())
    expect_identical(embedded_python, character())
})

test_that("installation does not declare a Python or compiled runtime", {
    desc <- utils::packageDescription("hpgR")
    entries <- function(field) {
        value <- desc[[field]]
        if (is.null(value) || is.na(value)) return(character())
        trimws(sub("[[:space:]]*\\(.*$", "", strsplit(value, ",", fixed = TRUE)[[1L]]))
    }
    runtime <- unique(c(entries("Depends"), entries("Imports"), entries("LinkingTo")))
    expect_false(any(c("reticulate", "rPython", "basilisk", "Rcpp", "cpp11",
        "RcppGSL") %in% runtime))
    expect_false("reticulate" %in% entries("Suggests"))
    expect_false(identical(tolower(desc[["NeedsCompilation"]]), "yes"))
})

test_that("the installed package payload has no Python files or native libraries", {
    root <- system.file(package = "hpgR")
    files <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    ignore_file <- file.path(root, ".Rbuildignore")
    if (file.exists(ignore_file)) {
        patterns <- readLines(ignore_file, warn = FALSE)
        patterns <- patterns[nzchar(patterns)]
        if (length(patterns)) files <- files[!vapply(files, function(path) {
            parts <- strsplit(path, "/", fixed = TRUE)[[1L]]
            ancestors <- vapply(seq_along(parts), function(i)
                paste(parts[seq_len(i)], collapse = "/"), character(1))
            any(vapply(patterns, function(pattern)
                any(grepl(pattern, ancestors, ignore.case = TRUE)), logical(1)))
        }, logical(1))]
    }
    # inst/ is flattened by installation; tests cover both the source and the
    # actual installed package path returned by pkgload/system.file.
    expect_identical(files[grepl("[.]py[cow]?$", files, ignore.case = TRUE)], character())
    expect_identical(files[grepl("[.](so|dll|dylib|pyd)$", files, ignore.case = TRUE)],
        character())
})

test_that("sequence inference and dating have implemented source records", {
    available <- c("inferAncestry", "dateAncestry")
    expect_true(all(available %in% getNamespaceExports("hpgR")))
    manifest <- sourceProvenance()
    for (name in available) {
        rows <- which(vapply(strsplit(manifest$functions, ";", fixed = TRUE),
            function(functions) name %in% functions, logical(1)))
        expect_true(length(rows) > 0L,
            info = paste(name, "must retain an honest source-status record"))
        if (length(rows)) {
            text <- paste(unlist(manifest[rows,
                intersect(c("status", "adaptation", "validation_scope"), names(manifest)),
                drop = FALSE]), collapse = " ")
            expect_match(tolower(text), "r (translation|implementation|adaptation)")
        }
    }
})

test_that("ancestry simulation and mutation return ordinary zero-based R tables", {
    ancestry <- simulateAncestry(2L, 50, 100, 0, seed = 19401L,
        recordFullArg = TRUE)
    expect_true(is.list(ancestry))
    expect_true(is.data.frame(ancestry$nodes))
    expect_true(is.data.frame(ancestry$edges))
    expect_identical(as.integer(ancestry$sampleIds), 0:3)
    expect_true(all(c("id", "time", "flags") %in% names(ancestry$nodes)))
    expect_true(all(c("left", "right", "parent", "child") %in% names(ancestry$edges)))
    expect_true(all(ancestry$nodes$id >= 0L))
    expect_true(all(ancestry$edges$parent %in% ancestry$nodes$id))
    expect_true(all(ancestry$edges$child %in% ancestry$nodes$id))
    age <- setNames(ancestry$nodes$time, ancestry$nodes$id)
    expect_true(all(age[as.character(ancestry$edges$parent)] >
        age[as.character(ancestry$edges$child)]))
    zero <- ancestryData(mutateAncestry(ancestry, 0, seed = 19402L))
    expect_identical(dim(zero$observed$haplotypes), c(4L, 0L))
    expect_length(zero$observed$positions, 0L)
    expect_equal(zero$observed$sequenceLength, 100)
    expect_equal(zero$nTrees, 1)
    expect_false(any(c("nodes", "edges", "events", "populationSize", "Ne") %in%
        names(zero$observed)))
    hasExternal <- function(x) {
        if (typeof(x) %in% c("externalptr", "environment")) return(TRUE)
        if (is.list(x)) return(any(vapply(x, hasExternal, logical(1))))
        FALSE
    }
    expect_false(hasExternal(ancestry))
    expect_false(hasExternal(zero))
})

test_that("pure-R ancestry calls are reproducible without changing caller state", {
    set.seed(19499L)
    before <- .Random.seed
    a <- simulateAncestry(2L, 20, 50, .001, seed = 19801L,
        recordFullArg = TRUE)
    expect_identical(.Random.seed, before)
    expect_identical(a, simulateAncestry(2L, 20, 50, .001,
        seed = 19801L, recordFullArg = TRUE))
    frozen <- serialize(a, NULL)
    m <- mutateAncestry(a, .002, seed = 19498L)
    expect_identical(.Random.seed, before)
    expect_identical(m, mutateAncestry(a, .002, seed = 19498L))
    expect_identical(serialize(a, NULL), frozen)
    replay <- hpgR:::msprimeAncestry(2L, 20, 50, .001,
        seed = 19497L, recordFullArg = TRUE, events = a$events)
    expect_identical(replay$nodes, a$nodes)
    expect_identical(replay$edges, a$edges)
    expect_identical(replay$eventCounts, a$eventCounts)
})

test_that("two-copy waiting times follow the neutral coalescent probability law", {
    Ne <- 37
    waits <- vapply(seq_len(400L), function(i) {
        a <- simulateAncestry(1L, Ne, 10, 0, seed = 19500L + i)
        max(a$nodes$time)
    }, numeric(1))
    expect_true(all(waits > 0))
    # Four standard errors of a 400-draw exponential mean; seeds and margins
    # are declared before inspecting results, not selected for a target tree.
    expect_lt(abs(mean(waits) / (2 * Ne) - 1), .20)
    expect_lt(abs(mean(waits <= 2 * Ne) - (1 - exp(-1))), .12)
})

test_that("mutations are Poisson events along a fixed sampled tree", {
    a <- simulateAncestry(1L, 50, 100, 0, seed = 19601L)
    ages <- setNames(a$nodes$time, a$nodes$id)
    exposure <- sum((a$edges$right - a$edges$left) *
        (ages[as.character(a$edges$parent)] - ages[as.character(a$edges$child)]))
    mu <- 20 / exposure
    mutants <- lapply(seq_len(300L), function(i)
        ancestryData(mutateAncestry(a, mu, seed = 19700L + i)))
    counts <- vapply(mutants, function(out) as.numeric(out$nSites), numeric(1))
    expect_true(all(vapply(mutants, function(out)
        nrow(out$observed$haplotypes) == 2L, logical(1))))
    expect_true(all(vapply(mutants, function(out)
        all(colSums(out$observed$haplotypes) == 1L), logical(1))))
    expect_true(all(vapply(mutants, function(out)
        all(out$observed$positions >= 0 & out$observed$positions < 100), logical(1))))
    expect_true(all(vapply(mutants, function(out)
        !anyDuplicated(out$observed$positions), logical(1))))
    expect_lt(abs(mean(counts) - 20), 1.25)
    expect_lt(abs(stats::var(counts) - 20), 6)
})

test_that("linked ancestral intervals cover the sequence without gaps", {
    out <- ancestryData(simulateAncestry(2L, 20, 50, .001,
        seed = 19801L, recordFullArg = TRUE))
    intervals <- out$intervals
    expect_gt(nrow(intervals), 0L)
    expect_equal(intervals$left[1L], 0)
    expect_equal(tail(intervals$right, 1L), 50)
    expect_true(all(intervals$right > intervals$left))
    if (nrow(intervals) > 1L)
        expect_equal(tail(intervals$left, -1L), head(intervals$right, -1L))
    expect_true(all(is.finite(intervals$tmrca)))
    expect_true(all(intervals$tmrca >= intervals$pair_time))
})

test_that("hand-checked interval trees give the correct mutation carriers and ages", {
    # [0,5): ((S0,S1),S2), with pair/root ages 2/5.
    # [5,10): ((S0,S2),S1), with pair/root ages 7/7.
    a <- structure(list(nodes = data.frame(id = 0:6,
        time = c(0, 0, 0, 2, 5, 3, 7), flags = c(1L, 1L, 1L, rep(0L, 4L))),
        edges = data.frame(left = rep(c(0, 5), each = 4L),
            right = rep(c(5, 10), each = 4L),
            parent = c(3L, 3L, 4L, 4L, 5L, 5L, 6L, 6L),
            child = c(0L, 1L, 3L, 2L, 0L, 2L, 5L, 1L)),
        sampleIds = 0:2, sequenceLength = 10, timeUnits = "generations"),
        class = c("hpgR_ancestry", "list"))
    events <- data.frame(edge = c(3L, 7L, 8L),
        position = c(1, 5, 9), time = c(3, 5, 2))
    m <- hpgR:::msprimeMutations(a, 0, seed = 19901L, events = events)
    out <- ancestryData(m)
    expect_identical(out$observed$haplotypes,
        rbind(c(1L, 1L, 0L), c(1L, 0L, 1L), c(0L, 1L, 0L)))
    expect_identical(out$observed$positions, c(1, 5, 9))
    expect_equal(out$intervals$left, c(0, 5))
    expect_equal(out$intervals$right, c(5, 10))
    expect_equal(out$intervals$pair_time, c(2, 7))
    expect_equal(out$intervals$tmrca, c(5, 7))
    expect_equal(out$nTrees, 2)
    expect_equal(out$nRecordedTrees, 2)
})

test_that("suppressing a unary display node preserves mutations on its branch", {
    a <- structure(list(nodes = data.frame(id = 0:3,
        time = c(0, 0, 1, 4), flags = c(1L, 1L, 0L, 0L)),
        edges = data.frame(left = 0, right = 10,
            parent = c(2L, 3L, 3L), child = c(0L, 2L, 1L)),
        sampleIds = 0:1, sequenceLength = 10, timeUnits = "generations"),
        class = c("hpgR_ancestry", "list"))
    m <- hpgR:::msprimeMutations(a, 0, seed = 19902L,
        events = data.frame(edge = 2L, position = 2, time = 2))
    out <- ancestryData(m)
    expect_identical(out$observed$haplotypes, matrix(c(1L, 0L), ncol = 1L))
    expect_false(2L %in% out$localTrees[[1L]]$node)
    expect_true(2L %in% out$nodes$id)
    expect_equal(out$intervals$pair_time, 4)
    expect_equal(out$intervals$tmrca, 4)
    bad <- a
    bad$edges <- rbind(bad$edges,
        data.frame(left = 0, right = 10, parent = 3L, child = 0L))
    expect_error(ancestryData(bad), "multiple parents")
})

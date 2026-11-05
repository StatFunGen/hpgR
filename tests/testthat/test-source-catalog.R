test_that("the supplied source catalog can be searched by package and symbol", {
    packages <- sourceCatalog(table = "packages")
    expect_setequal(packages$package,
        c("ancestryTools", "ARGtools", "demographyTools", "haplotypeTools", "tskitR"))
    symbols <- sourceCatalog(".truncated_sfs", package = "demographyTools")
    expect_true(any(symbols$symbol == ".truncated_sfs"))
    expect_true(all(symbols$package == "demographyTools"))
    expect_equal(sourceCatalog(".TRUNCATED_SFS", package = "demographyTools"), symbols)
    files <- sourceCatalog("R/demography_sfs.R", package = "demographyTools", table = "files")
    expect_identical(nrow(files), 1L)
    expect_identical(nchar(files$sha256), 64L)
    expect_identical(nrow(sourceCatalog("not_a_source_symbol_8149")), 0L)
    expect_error(sourceCatalog(NA_character_), "query must be one string")
})

test_that("implementation types are recorded separately from validation scope", {
    path <- system.file("provenance", "modules.tsv", package = "hpgR")
    modules <- utils::read.delim(path, stringsAsFactors = FALSE)
    expect_false(anyDuplicated(modules$file) > 0L)
    expect_true(all(nzchar(modules$validation_scope)))
    expect_identical(modules$kind[modules$file == "simupop.R"], "upstream")
    expect_match(modules$source_ids[modules$file == "simupop.R"], "transmission", fixed = TRUE)
    expect_match(modules$validation_scope[modules$file == "weighted_ld.R"],
                 "strict identity remains open", fixed = TRUE)
})

test_that("original course code and upstream implementations have distinct files", {
    modules <- utils::read.delim(system.file("provenance", "modules.tsv", package = "hpgR"),
                               stringsAsFactors = FALSE)
    sources <- sourceProvenance()
    expect_true(all(modules$kind %in% c("upstream", "pedagogical", "documentation")))
    executable <- modules$kind != "documentation"
    expect_identical(startsWith(modules$file[executable], "pedagogical_"),
                     modules$kind[executable] == "pedagogical")
    expect_identical(modules$file[modules$kind == "documentation"], "hpgR-package.R")
    expect_true(all(sources$upstream %in% c("haplotypeTools", "ancestryTools",
        "demographyTools", "ARGtools", "tskitR", "simuPOP")))
    expect_true(all(grepl("^[[:xdigit:]]{64}$", sources$sha256)))
    for (i in which(modules$kind == "upstream")) {
        ids <- strsplit(modules$source_ids[i], ";", fixed = TRUE)[[1L]]
        expect_true(length(ids) > 0L && all(nzchar(ids)) && all(ids %in% sources$id),
                    info = modules$file[i])
    }
})

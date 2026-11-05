# Rebuild the R sequence-control datasets.
# Rscript --vanilla tools/generate_pedigree_control.R --output-dir DIR --truth-dir PRIVATE_DIR
# Both directories must be explicit. Existing artifact files are not overwritten.

controlChecksum <- function(path) unname(tools::md5sum(path))

pgBuildP1 <- function(core_path, output_dir, truth_dir) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite needed for manifest")
    builder_path <- normalizePath(file.path(dirname(core_path), "../tools/generate_pedigree_control.R"))
    core_sha <- controlChecksum(core_path)
    builder_sha <- controlChecksum(builder_path)
    library(hpgR)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(truth_dir, recursive = TRUE, showWarnings = FALSE)
    public <- normalizePath(output_dir)
    private <- normalizePath(truth_dir)
    if (public == private || startsWith(private, paste0(public, "/")))
        stop("truth directory must be outside the observed-data directory")
    files <- file.path(public, c("observed.rds", "calibration.rds", "manifest.json"))
    private_files <- file.path(private, c("truth.rds", "run_report.json"))
    if (any(file.exists(c(files, private_files)))) stop("artifact exists; use fresh directories")
    histories <- list(D1 = rep(40L, 61L), D2 = rep(80L, 61L),
                      D3 = c(80L, rep(20L, 40L), rep(80L, 20L)),
                      D4 = rep(80L, 61L))
    observed <- calibration <- truth <- list()
    report <- list()
    for (i in seq_along(histories)) {
        id <- names(histories)[i]
        linked <- i == 4L
        lengths <- if (linked) 1000 else rep(1000, 3L)
        # Scaled linked extension: 3 cM per 1 kb makes ancestry changes visible.
        # Every dataset and region has a fixed ancestry and mutation seed.
        maps <- rep(if (linked) 0.03 else 0, length(lengths))
        gc(reset = TRUE)
        timing <- system.time({
            result <- pgSimulate(histories[[i]], lengths, n = 12L, mu = 2e-5,
                                  map_lengths = maps, ancestry_seed = 814900L + 1000L * i,
                                  mutation_seed = 914900L + 1000L * i)
        })
        memory <- gc()
        replay_time <- system.time({
            for (j in seq_along(result$truth$regions)) {
                region <- result$truth$regions[[j]]
                stopifnot(identical(unname(result$observed$regions[[j]]$haplotypes),
                                    pgReplay(region$ancestry, region$mutations)))
            }
        })
        observed[[id]] <- result$observed
        calibration[[id]] <- result$calibration
        truth[[id]] <- result$truth
        report[[id]] <- list(
            elapsed_seconds = unname(timing["elapsed"]),
            replay_seconds = unname(replay_time["elapsed"]),
            retained_object_bytes = as.numeric(object.size(result)),
            R_gc_heap_high_water_bound_MiB = sum(memory[, ncol(memory)]),
            summaries = pgObservedSummaries(result$observed, result$calibration)$regions,
            ancestral_extensions = vapply(result$truth$regions,
                function(x) x$ancestry$extension_generations, integer(1)),
            local_tree_intervals = vapply(result$truth$regions,
                function(x) length(x$ancestry$local_trees), integer(1)))
        message(id, " generated in ", round(timing["elapsed"], 2),
                " s; independent replay passed")
    }
    seeds <- unlist(lapply(truth, function(d) lapply(d$regions,
        function(r) c(ancestry = r$ancestry$ancestry_seed, mutation = r$mutation_seed))))
    stopifnot(!anyDuplicated(seeds))
    stopifnot(identical(core_sha, controlChecksum(core_path)),
              identical(builder_sha, controlChecksum(builder_path)))
    saveRDS(observed, files[1], version = 2)
    saveRDS(calibration, files[2], version = 2)
    saveRDS(truth, private_files[1], version = 2)
    manifest <- list(schema_version = 1L,
        role = "R-generated sequence controls",
        files = setNames(lapply(files[1:2], function(f)
            list(md5 = controlChecksum(f), bytes = unname(file.info(f)$size))),
            basename(files[1:2])),
        source = list(core = basename(core_path), md5 = core_sha,
                      builder = basename(builder_path), builder_md5 = builder_sha),
        package = "hpgR", package_version = as.character(packageVersion("hpgR")),
        language = "R", R = R.version.string, RNG = c("Mersenne-Twister", "Inversion", "Rejection"),
        dataset_ids = names(observed), truth_stored_separately = TRUE,
        generations = "sampling at t=0; true g and node mappings in separate validation record",
        design = "three independent-region datasets plus a separate linked-region extension")
    jsonlite::write_json(manifest, files[3], pretty = TRUE, auto_unbox = TRUE)
    private_report <- list(
        source_md5 = controlChecksum(core_path), R = R.version.string,
        history_labels = list(D1 = "constant40", D2 = "constant80",
                              D3 = "80 then20for40 then80for20", D4 = "constant80linked"),
        truth_md5 = controlChecksum(private_files[1]), report = report, seed_ledger = seeds,
        root_rule = "prepend older discrete-WF maps; no seed rejection or later sample",
        status = "R generation and exact event replay passed")
    jsonlite::write_json(private_report, private_files[2], pretty = TRUE, auto_unbox = TRUE)
    invisible(list(public = files, private = private_files, report = private_report))
}

if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    option <- function(name) {
        at <- match(name, args)
        if (is.na(at) || at == length(args)) stop(paste("missing", name))
        args[at + 1L]
    }
    script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
    core <- file.path(dirname(normalizePath(script)), "../R/pedagogical_wright_fisher_genealogies.R")
    pgBuildP1(normalizePath(core), option("--output-dir"), option("--truth-dir"))
}

#!/usr/bin/env Rscript
# Rebuild the population-history observations and prediction cache in R.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || length(args) > 2L)
    stop("Usage: Rscript tools/rebuild_population_cache.R OUTPUT_DIRECTORY [WORKERS]")
library(hpgR)
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite for the manifest")
directory <- normalizePath(args[1L], mustWork = FALSE)
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
workers <- if (length(args) == 2L) as.integer(args[2L]) else 1L
specification <- populationHistoryDesign()
cat("Generating observed DNA with seed", specification$design$observation_seed, "\n")
simulation <- simulatePopulationHistory(
    specification$candidates$recent_founding, specification$design,
    specification$design$observation_seed, recordEvery = 1L)
saveRDS(simulation$observed, file.path(directory, "observed.rds"), version = 3)
saveRDS(simulation$calibration, file.path(directory, "calibration.rds"), version = 3)
replicateCount <- sum(vapply(specification$banks, function(b)
    length(b$training_seeds) + length(b$calibration_seeds), integer(1)))
cat("Generating", replicateCount, "independent R replicates using", workers, "workers\n")
started <- proc.time()
predictions <- populationHistoryPredictions(specification, workers)
saveRDS(predictions, file.path(directory, "history_predictions.rds"), version = 3)
analysis <- populationHistoryExample(simulation$observed, simulation$calibration, predictions)
saveRDS(analysis, file.path(directory, "analysis.rds"), version = 3)
saveRDS(simulation, file.path(directory, "generation_records.rds"), version = 3)
design <- specification$design
manifest <- list(schema_version = 2L, engine = "hpgR", language = "R",
    package_version = as.character(packageVersion("hpgR")),
    canonical_notebook = "P8149_Lecture_9_Effective_Size_and_History.ipynb",
    files = c("observed.rds", "calibration.rds", "history_predictions.rds"),
    generation_design = design, candidates = specification$candidates,
    prediction_replicates_per_candidate = length(specification$banks[[1L]]$training_seeds),
    calibration_replicates_per_candidate = length(specification$banks[[1L]]$calibration_seeds),
    seed_ledger = predictions$seed_ledger,
    observation = list(diploid_individuals = design$sample_diploids,
        phased_haplotypes = 2L * design$sample_diploids,
        chromosomes = design$chromosomes,
        modeled_sites = design$chromosomes * design$sites_per_chromosome),
    mutation_model = "symmetric biallelic mutation before reproduction",
    reproduction = "uniform independent diploid parents with replacement",
    recombination = "Haldane interval switching; independent chromosome starts",
    founder_sampling = "diploid individuals without replacement",
    dependence = "chromosomes share individuals and the realized pedigree",
    rebuild = "Rscript tools/rebuild_population_cache.R OUTPUT_DIRECTORY WORKERS")
jsonlite::write_json(manifest, file.path(directory, "manifest.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("Elapsed seconds:", unname((proc.time() - started)["elapsed"]), "\n")
print(analysis$history_fit$profile)
print(analysis$history_check$profile)
cat("Completed R cache at", directory, "\n")

#!/usr/bin/env Rscript
# Check the maintained public inputs, code records and figure ownership.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Supply the P8149_2026_vignettes directory")
site <- normalizePath(args[1L], mustWork = TRUE)
library(hpgR)
popgen_files <- c("P8149_Lecture_8_Sampling_Drift_Mutation.ipynb",
    "P8149_Lecture_8_Effective_Size_Generation_Records.ipynb",
    "P8149_Lecture_9_Effective_Size_and_History.ipynb",
    "P8149_Lecture_9_Two_DNA_Copies.ipynb",
    "P8149_Lecture_9_Mutations_on_a_Tree.ipynb",
    "P8149_Lecture_9_Ancestral_Graphs.ipynb")
files <- sort(list.files(site, '^P8149_Lecture_.*[.]ipynb$'))
stopifnot(all(popgen_files %in% files))
books <- setNames(lapply(file.path(site, files), jsonlite::fromJSON,
    simplifyVector = FALSE), files)
source_text <- function(cell) paste(cell$source, collapse = "")
code_cells <- function(book) Filter(function(c) c$cell_type == "code", book$cells)
sha_program <- Sys.which("shasum")
sha_args <- c("-a", "256")
if (!nzchar(sha_program)) {
    sha_program <- Sys.which("sha256sum")
    sha_args <- character()
}
if (!nzchar(sha_program)) stop("A SHA-256 command is required for source records")
sha <- function(path = NULL, text = NULL) {
    out <- if (is.null(path)) system2(sha_program, sha_args, input = text, stdout = TRUE) else
        system2(sha_program, c(sha_args, shQuote(path)), stdout = TRUE)
    strsplit(out[1L], "[[:space:]]+")[[1L]][1L]
}
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
paths <- list.files(file.path(root, "R"), "[.]R$", full.names = TRUE)
source_hashes <- setNames(lapply(paths, sha), paste0("R/", basename(paths)))
rows <- list()
for (name in names(books)) {
    book <- books[[name]]
    cells <- code_cells(book)
    code <- vapply(cells, source_text, "")
    stopifnot(identical(book$metadata$kernelspec$name, "ir"),
        identical(book$metadata$kernelspec$language, "R"),
        identical(book$metadata$hpgR_version, as.character(packageVersion("hpgR"))))
    invisible(lapply(code, function(text) parse(text = text)))
    imports <- sum(vapply(gregexpr("library(hpgR)", code, fixed = TRUE),
        function(x) sum(x > 0L), 0L))
    stopifnot(imports == 1L,
        !any(trimws(code) == "library(hpgR)"),
        !any(grepl(".p7_eval(", code, fixed = TRUE)),
        all(vapply(cells[nzchar(trimws(code))], function(c) !is.null(c$execution_count), TRUE)))
    plotting_helpers <- c('draw_figure', 'draw_tree', 'draw_carriers',
        'clade_nodes', 'tree_positions', 'draw_recorded_tree')
    for (part in code) for (expr in parse(text = part)) {
        if (is.call(expr) && identical(expr[[1L]], as.name('<-')) &&
            is.symbol(expr[[2L]]) && is.call(expr[[3L]]) &&
            identical(expr[[3L]][[1L]], as.name('function')) &&
            !as.character(expr[[2L]]) %in% plotting_helpers)
            stop('A reusable implementation remains in ', name, ': ', expr[[2L]])
    }
    packed <- as.character(jsonlite::toJSON(lapply(cells,
        function(c) list(c$id, source_text(c))), auto_unbox = TRUE))
    stopifnot(identical(book$metadata$notebook_code_sha256, sha(text = packed)),
        setequal(names(book$metadata$hpgR_source_sha256), names(source_hashes)))
    for (path in names(source_hashes))
        stopifnot(identical(book$metadata$hpgR_source_sha256[[path]], source_hashes[[path]]))
    inputs <- book$metadata$notebook_input_sha256
    for (path in names(inputs)) stopifnot(identical(inputs[[path]], sha(file.path(site, path))))
    outputs <- unlist(lapply(cells, function(c) c$outputs), recursive = FALSE)
    errors <- Filter(function(o) identical(o$output_type, "error"), outputs)
    warnings <- Filter(function(o) identical(o$output_type, "stream") &&
        any(grepl("(^|\\n)(Warning|Error)(:| in| message)",
            paste(o$text, collapse = ""))), outputs)
    figures <- sum(vapply(outputs, function(o) "image/png" %in% names(o$data), TRUE))
    stopifnot(!length(errors), !length(warnings))
    rows[[name]] <- list(notebook = name, code_cells = length(cells),
        figures = figures, errors = list(), warnings = list(),
        package_version = as.character(packageVersion("hpgR")))
}

public_data <- function(name) {
    cell <- Filter(function(c) c$id == "public-data", books[[name]]$cells)
    stopifnot(length(cell) == 1L)
    env <- new.env(parent = .GlobalEnv)
    eval(parse(text = source_text(cell[[1L]])), envir = env)
    env
}
population <- public_data(popgen_files[3L])
pair <- public_data(popgen_files[4L])
tree <- public_data(popgen_files[5L])
stopifnot(identical(pair$observed, tree$observed),
    identical(pair$observed, readRDS(file.path(site, "data/popgen_forward_p1/observed.rds"))),
    identical(tree$calibration, readRDS(file.path(site, "data/popgen_forward_p1/calibration.rds"))))
for (name in c("observed", "calibration"))
    stopifnot(identical(get(name, population),
        readRDS(file.path(site, "data/population_history", paste0(name, ".rds")))))
stopifnot(identical(population$prediction_data,
    readRDS(file.path(site, "data/population_history/history_predictions.rds"))))
owners <- vapply(books, function(book) any(vapply(code_cells(book), function(c)
    grepl('draw_figure("p7_common_sample_validation"', source_text(c), fixed = TRUE), TRUE)), TRUE)
stopifnot(identical(names(owners)[owners], popgen_files[3L]),
    !file.exists(file.path(site, "P8149_2026_Lecture_9_Simulation_Inference.ipynb")))

# The browser only renders frames calculated by the package.
widget <- readLines(file.path(site, 'widgets/hmm-basics.mjs'), warn = FALSE)
data_line <- widget[startsWith(widget, 'const demonstrations = ')]
stopifnot(length(data_line) == 1L,
    !any(grepl('forwardBackward|function normalize|Math[.]random', widget)))
frames <- jsonlite::fromJSON(sub(';$', '', sub('^const demonstrations = ', '', data_line)))
for (i in seq_along(frames$switch$values))
    stopifnot(max(abs(frames$switch$posterior[i, ] -
        binaryHMM(frames$observed, frames$switch$values[i], .15))) < 1e-12)
for (i in seq_along(frames$error$values))
    stopifnot(max(abs(frames$error$posterior[i, ] -
        binaryHMM(frames$observed, .08, frames$error$values[i]))) < 1e-12)

if (length(args) > 2L) {
    course <- normalizePath(args[3L], mustWork = TRUE)
    visited <- character()
    figures <- character()
    matches <- function(text, pattern) regmatches(text, gregexpr(pattern, text, perl = TRUE))[[1L]]
    visit <- function(path) {
        if (path %in% visited || !file.exists(file.path(course, path))) return(invisible(NULL))
        visited <<- c(visited, path)
        lines <- readLines(file.path(course, path), warn = FALSE)
        text <- paste(lines[!grepl("^[[:space:]]*%", lines)], collapse = "\n")
        text <- gsub("(?s)\\\\begin\\{comment\\}.*?\\\\end\\{comment\\}", "", text, perl = TRUE)
        includes <- matches(text, "\\\\(?:input|include)\\{[^}]+\\}")
        for (item in includes) {
            child <- sub("^.*\\{([^}]+)\\}$", "\\1", item)
            if (!grepl("[.]tex$", child)) child <- paste0(child, ".tex")
            visit(child)
        }
        image <- matches(text, "\\\\includegraphics(?:\\[[^]]*\\])?\\{[^}]+\\}")
        for (item in image) {
            value <- sub("^.*\\{([^}]+)\\}$", "\\1", item)
            if (startsWith(value, "figures/coalescent/simulation_inference/"))
                figures <<- unique(c(figures, sub("[.]pdf$", "", basename(value))))
        }
    }
    visit("lectures/columbia_2026_P8149/P8149_2026_Lecture_9.tex")
    textbook <- c(p7_pritchard_f2_57_complete = "Pritchard, Fig. 2.57")
    for (figure in setdiff(figures, names(textbook))) {
        owner <- vapply(books, function(book) sum(vapply(code_cells(book), function(c)
            grepl(paste0('draw_figure("', figure, '"'), source_text(c), fixed = TRUE), TRUE)), 0L)
        if (sum(owner) != 1L) stop("Figure needs one public generating cell: ", figure)
    }
    cat("Active lecture figure sources:", length(setdiff(figures, names(textbook))),
        "notebook outputs and", length(intersect(figures, names(textbook))),
        "credited textbook crop checked across", length(visited), "TeX inputs.\n")
}

result <- unname(rows)
if (length(args) > 1L) jsonlite::write_json(result, args[2L], auto_unbox = TRUE, pretty = TRUE)
cat('PASS:', length(files), 'independently executed R notebooks, matching inputs and source hashes.\n')
cat("Code cells", sum(vapply(rows, `[[`, 0L, "code_cells")),
    "; figures", sum(vapply(rows, `[[`, 0L, "figures")), "\n")
cat("Shared genealogy observations and the single closing-figure source match.\n")

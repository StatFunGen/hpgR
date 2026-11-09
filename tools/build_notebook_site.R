#!/usr/bin/env Rscript
# Build only the notebooks selected for the course website.
main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) && !identical(args, '--stage-only'))
        stop('Usage: Rscript tools/build_notebook_site.R [--stage-only]')
    if (!requireNamespace('jsonlite', quietly = TRUE)) stop('Install jsonlite')
    script <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())])
    root <- normalizePath(file.path(dirname(script), '..'), mustWork = TRUE)
    source <- file.path(root, 'notebooks')
    baseUrl <- sub('/+$', '', Sys.getenv('BASE_URL', ''))
    if (nzchar(baseUrl) && !grepl('^(/[A-Za-z0-9._-]+)+$', baseUrl))
        stop('BASE_URL must be an absolute URL path, such as /hpgR')
    sha <- function(path = NULL, text = NULL) {
        program <- Sys.which('shasum')
        arguments <- c('-a', '256')
        if (!nzchar(program)) {
            program <- Sys.which('sha256sum')
            arguments <- character()
        }
        if (!nzchar(program)) stop('Install a SHA-256 command')
        output <- if (is.null(path)) system2(program, arguments, input = text, stdout = TRUE) else
            system2(program, c(arguments, shQuote(path)), stdout = TRUE)
        strsplit(output[1L], '[[:space:]]+')[[1L]][1L]
    }
    packageFiles <- list.files(file.path(root, 'R'), '[.]R$', full.names = TRUE)
    if (!length(packageFiles)) stop('Missing hpgR source files')
    packageHashes <- setNames(lapply(packageFiles, sha), paste0('R/', basename(packageFiles)))
    catalog <- jsonlite::fromJSON(file.path(source, 'catalog.json'), simplifyVector = FALSE)
    allEntries <- catalog$notebooks
    files <- vapply(allEntries, `[[`, '', 'file')
    if (anyDuplicated(files) || any(basename(files) != files) ||
        !setequal(files, list.files(source, '[.]ipynb$')))
        stop('The notebook catalog must describe each maintained notebook once')
    selected <- Filter(function(entry) isTRUE(entry$website), allEntries)
    if (!length(selected)) stop('The website publication list is empty')
    if (any(vapply(selected, function(entry)
        !entry$lecture %in% c(2L, 3L, 8L, 9L) ||
        !identical(entry$role, 'teaching_demonstration') ||
        grepl('homework|answer|solution|midterm|instructor', entry$file, ignore.case = TRUE), TRUE)))
        stop('A notebook is outside the approved publication scope')
    selectedFiles <- vapply(selected, `[[`, '', 'file')
    withheld <- setdiff(files, selectedFiles)
    build <- file.path(source, '_build')
    dir.create(build, showWarnings = FALSE)
    stage <- tempfile('site-', tmpdir = build)
    dir.create(stage)
    on.exit(if (dir.exists(stage)) message('Build files: ', stage), add = TRUE)
    copy <- function(relative, destination = relative) {
        path <- normalizePath(file.path(source, relative), mustWork = TRUE)
        if (!startsWith(path, paste0(source, '/'))) stop('Input outside notebooks: ', relative)
        output <- file.path(stage, destination)
        dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
        if (!file.copy(path, output)) stop('Could not copy ', relative)
    }
    quoteYaml <- function(text) as.character(jsonlite::toJSON(text, auto_unbox = TRUE))
    title <- function(book) {
        markdown <- Filter(function(cell) cell$cell_type == 'markdown', book$cells)
        lines <- unlist(strsplit(paste(vapply(markdown, function(cell)
            paste(cell$source, collapse = ''), ''), collapse = '\n'), '\n'))
        headings <- lines[grepl('^# ', lines)]
        if (!length(headings)) stop('Notebook needs a title')
        sub('^# ', '', headings[1L])
    }
    datasets <- catalog$datasets
    dataFiles <- vapply(datasets, `[[`, '', 'file')
    if (anyDuplicated(dataFiles) ||
        any(!grepl('^data/[A-Za-z0-9_/-]+[.](rds|csv)$', dataFiles)))
        stop('Invalid data download list')
    dependencies <- dataFiles
    titles <- character()
    for (entry in selected) {
        book <- jsonlite::fromJSON(file.path(source, entry$file), simplifyVector = FALSE)
        code <- Filter(function(cell) cell$cell_type == 'code', book$cells)
        packed <- jsonlite::toJSON(lapply(code, function(cell)
            list(cell$id, paste(cell$source, collapse = ''))), auto_unbox = TRUE)
        recorded <- book$metadata$hpgR_source_sha256
        inputs <- book$metadata$notebook_input_sha256
        if (!setequal(names(recorded), names(packageHashes)) ||
            !all(vapply(names(packageHashes), function(path)
                identical(recorded[[path]], packageHashes[[path]]), TRUE)) ||
            !identical(book$metadata$notebook_code_sha256, sha(text = packed)) ||
            !all(vapply(names(inputs), function(path)
                identical(inputs[[path]], sha(file.path(source, path))), TRUE)))
            stop('Refresh changed notebook code or inputs with tools/sync_vignettes.R: ', entry$file)
        titles[entry$file] <- title(book)
        executable <- Filter(function(cell) cell$cell_type == 'code' &&
            nzchar(trimws(paste(cell$source, collapse = ''))), book$cells)
        if (!identical(book$metadata$kernelspec$name, 'ir') ||
            any(vapply(executable, function(cell) is.null(cell$execution_count) ||
                any(vapply(cell$outputs, function(output)
                    identical(output$output_type, 'error'), TRUE)), TRUE)))
            stop('Notebook has unverified supplied code: ', entry$file)
        markdown <- paste(vapply(Filter(function(cell) cell$cell_type == 'markdown',
            book$cells), function(cell) paste(cell$source, collapse = ''), ''), collapse = '\n')
        local <- regmatches(markdown, gregexpr(
            '(?:data|assets|widgets)/[A-Za-z0-9_./-]+', markdown, perl = TRUE))[[1L]]
        dependencies <- unique(c(dependencies, local,
            names(book$metadata$notebook_input_sha256)))
        # Only publication metadata is added to the staged reading copy.
        # Download links point to the byte-identical maintained notebook.
        copy(entry$file, file.path('downloads', entry$file))
        if (!identical(book$cells[[1L]]$cell_type, 'markdown'))
            stop('The first notebook cell must contain its title')
        first <- paste(book$cells[[1L]]$source, collapse = '')
        if (startsWith(first, '---')) stop('Merge existing frontmatter before publication')
        url <- paste0('https://github.com/StatFunGen/hpgR/')
        header <- c('---', paste0('title: ', quoteYaml(titles[entry$file])),
            paste0('source_url: ', url, 'blob/main/notebooks/', entry$file),
            paste0('edit_url: ', url, 'edit/main/notebooks/', entry$file),
            'downloads:', paste0('  - file: downloads/', entry$file),
            '    title: Download R notebook')
        dataFiles <- names(book$metadata$notebook_input_sha256)
        for (path in dataFiles[startsWith(dataFiles, 'data/')]) {
            dataTitle <- switch(basename(path),
                observed.rds = 'Observed DNA',
                calibration.rds = 'Rate and map calibration',
                history_predictions.rds = 'Simulated history predictions',
                basename(path))
            header <- c(header, paste0('  - file: ', path),
                paste0('    title: ', quoteYaml(dataTitle)))
        }
        header <- c(header, '---', '')
        book$cells[[1L]]$source <- paste(c(header, first), collapse = '\n')
        jsonlite::write_json(book, file.path(stage, entry$file), auto_unbox = TRUE,
            pretty = TRUE, null = 'null', digits = NA)
    }
    if (any(!grepl('^(data|assets|widgets)/', dependencies)) ||
        any(grepl('(^|/)[.][.](/|$)', dependencies))) stop('Unsafe asset path')
    for (path in dependencies) copy(path)
    copy('styles.css')
    home <- readLines(file.path(source, 'index.md'), warn = FALSE)
    toc <- c('  toc:', '    - file: index.md')
    for (lecture in sort(unique(vapply(selected, `[[`, 0L, 'lecture')))) {
        entries <- Filter(function(entry) entry$lecture == lecture, selected)
        toc <- c(toc, paste0('    - title: Lecture ', lecture), '      children:')
        home <- c(home, '', paste0('## Lecture ', lecture), '')
        if (lecture == 9L)
            home <- c(home, 'Read 1 and 2 first. In notebook 3, start at **Fit a population size history from the same sample** for the classroom example. Notebook 4 is for self-study.', '')
        home <- c(home, '::::{grid} 1 1 2 2', ':class: p8149-notebook-grid', '')
        for (i in seq_along(entries)) {
            entry <- entries[[i]]
            label <- titles[entry$file]
            if (lecture == 9L) label <- paste0(i, '. ', label)
            toc <- c(toc, paste0('        - file: ', entry$file),
                paste0('          title: ', quoteYaml(label)))
            home <- c(home, paste0(':::{card} ', label),
                paste0(':link: ', entry$file), '', 'R notebook', ':::', '')
        }
        home <- c(home, '::::')
    }
    home <- c(home, '', '## Data', '')
    for (entry in datasets)
        home <- c(home, paste0('- [', entry$title, '](', entry$file, ')'))
    writeLines(home, file.path(stage, 'index.md'), useBytes = TRUE)
    config <- readLines(file.path(source, 'myst.yml'), warn = FALSE)
    position <- which(config == 'site:')
    if (length(position) != 1L) stop('Expected one site configuration')
    config <- append(config, c(toc, ''), after = position - 1L)
    writeLines(config, file.path(stage, 'myst.yml'), useBytes = TRUE)
    staged <- list.files(stage, recursive = TRUE)
    expected <- c(selectedFiles, file.path('downloads', selectedFiles),
        dependencies, 'styles.css', 'index.md', 'myst.yml')
    stopifnot(setequal(staged, expected))
    if ('--stage-only' %in% args) {
        cat('PASS: isolated publication source, ', length(selectedFiles), ' notebooks.\n',
            stage, '\n', sep = '')
        return(invisible(stage))
    }
    program <- file.path(source, '.pixi/envs/default/bin/jupyter-book')
    if (!file.exists(program)) program <- Sys.which('jupyter-book')
    if (!nzchar(program)) stop('Run pixi install in notebooks/ first')
    old <- getwd()
    on.exit(setwd(old), add = TRUE)
    setwd(stage)
    log <- file.path(stage, 'build.log')
    status <- system2(program, c('build', '--html', '--strict'), stdout = log, stderr = log)
    if (!identical(status, 0L)) {
        cat(tail(readLines(log, warn = FALSE), 60L), sep = '\n')
        stop('Website build failed; see ', log)
    }
    html <- file.path(stage, '_build/html')
    if (!file.exists(file.path(html, 'index.html'))) stop('Missing website index')
    outputFiles <- list.files(html, recursive = TRUE, full.names = TRUE)
    forbidden <- unique(c(withheld, tolower(gsub('_', '-', sub('[.]ipynb$', '', withheld))),
        'P8149_2026_vignettes', 'Midterm_Part_B', 'instructor_truth.rds'))
    textFiles <- outputFiles[grepl('[.](html|json|js|mjs|md|txt|ipynb)$', outputFiles)]
    for (path in textFiles) {
        contents <- paste(readLines(path, warn = FALSE), collapse = '\n')
        if (grepl('(?i)(?:P8149[^"<>[:space:]]*Homework[^"<>[:space:]]*[.]ipynb|p8149-(?:lecture-5-)?homework-)',
            contents, perl = TRUE)) stop('Assignment notebook in ', path)
        found <- forbidden[vapply(forbidden, function(term)
            grepl(term, contents, fixed = TRUE), TRUE)]
        if (length(found)) stop('Unpublished material in ', path, ': ', paste(found, collapse = ', '))
    }
    status <- system2(file.path(R.home('bin'), 'Rscript'),
        c(shQuote(file.path(root, 'tools/verify_notebook_site.R')), shQuote(stage)))
    if (!identical(status, 0L)) stop('Website download and link verification failed')
    jsonlite::write_json(list(notebooks = selectedFiles, assets = dependencies,
        base_url = baseUrl,
        rendered_pages = length(list.files(html, '[.]html$', recursive = TRUE)),
        withheld_notebooks = length(withheld),
        forbidden_matches = 0L), file.path(stage, 'publication-check.json'),
        auto_unbox = TRUE, pretty = TRUE)
    final <- file.path(build, 'site-source')
    previous <- file.path(build, 'previous-site-source')
    if (dir.exists(previous)) unlink(previous, recursive = TRUE)
    if (dir.exists(final) && !file.rename(final, previous)) stop('Could not retain previous preview')
    setwd(old)
    if (!file.rename(stage, final)) stop('Could not save the completed site')
    cat('PASS: ', length(selectedFiles), ' notebooks; ', length(withheld),
        ' withheld; publication scan clear.\n',
        'Preview: cd notebooks && pixi run preview\n',
        'HTML: ', file.path(final, '_build/html'), '\n', sep = '')
}
main()

#!/usr/bin/env Rscript
# Check published pages, original downloads and required local assets.
main <- function() {
    if (!requireNamespace('jsonlite', quietly = TRUE)) stop('Install jsonlite')
    args <- commandArgs(trailingOnly = TRUE)
    script <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())])
    root <- normalizePath(file.path(dirname(script), '..'), mustWork = TRUE)
    source <- file.path(root, 'notebooks')
    stage <- normalizePath(if (length(args)) args[1L] else
        file.path(source, '_build/site-source'), mustWork = TRUE)
    html <- file.path(stage, '_build/html')
    baseUrl <- sub('/+$', '', Sys.getenv('BASE_URL', ''))
    readJson <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
    catalog <- readJson(file.path(source, 'catalog.json'))
    entries <- catalog$notebooks
    selected <- Filter(function(entry) isTRUE(entry$website), entries)
    files <- vapply(selected, `[[`, '', 'file')
    slugs <- tolower(gsub('_', '-', sub('[.]ipynb$', '', files)))
    config <- readJson(file.path(stage, '_build/site/config.json'))
    pages <- Filter(function(page) !is.null(page$slug), config$projects[[1L]]$pages)
    stopifnot(setequal(slugs, vapply(pages, `[[`, '', 'slug')),
        length(pages) == length(files))
    hash <- function(path) {
        program <- Sys.which('shasum')
        arguments <- c('-a', '256', shQuote(path))
        if (!nzchar(program)) {
            program <- Sys.which('sha256sum')
            arguments <- shQuote(path)
        }
        if (!nzchar(program)) stop('Install a SHA-256 command')
        strsplit(system2(program, arguments, stdout = TRUE)[1L], '[[:space:]]+')[[1L]][1L]
    }
    resolve <- function(url) {
        path <- sub('[?#].*$', '', url)
        if (nzchar(baseUrl) && (identical(path, baseUrl) ||
            startsWith(path, paste0(baseUrl, '/')))) {
            path <- substring(path, nchar(baseUrl) + 1L)
        }
        path <- URLdecode(sub('^/', '', path))
        if (grepl('(^|/)[.][.](/|$)', path)) stop('Unsafe output URL: ', url)
        output <- file.path(html, path)
        if (dir.exists(output)) output <- file.path(output, 'index.html')
        output
    }
    localUrls <- function(node) {
        if (!is.list(node)) return(character())
        own <- as.character(unlist(node[intersect(names(node), c('url', 'src', 'esm'))],
            use.names = FALSE))
        own <- own[vapply(own, is.character, TRUE)]
        unique(c(own[startsWith(own, '/') & !startsWith(own, '//')],
            unlist(lapply(node, localUrls), use.names = FALSE)))
    }
    downloads <- 0L
    inputs <- 0L
    localLinks <- character()
    for (i in seq_along(files)) {
        name <- files[i]
        page <- readJson(file.path(html, paste0(slugs[i], '.json')))
        book <- readJson(file.path(source, name))
        expected <- c(name, names(book$metadata$notebook_input_sha256))
        expected <- expected[expected == name | startsWith(expected, 'data/')]
        for (path in expected) {
            matches <- Filter(function(download) identical(download$filename,
                basename(path)), page$frontmatter$downloads)
            stopifnot(length(matches) == 1L,
                startsWith(matches[[1L]]$url, paste0(baseUrl, '/')),
                identical(hash(resolve(matches[[1L]]$url)), hash(file.path(source, path))))
            downloads <- downloads + as.integer(path == name)
            inputs <- inputs + as.integer(path != name)
        }
        stopifnot(file.exists(file.path(html, slugs[i], 'index.html')),
            identical(page$frontmatter$source_url,
                paste0('https://github.com/StatFunGen/hpgR/blob/main/notebooks/', name)),
            identical(page$frontmatter$edit_url,
                paste0('https://github.com/StatFunGen/hpgR/edit/main/notebooks/', name)))
        localLinks <- unique(c(localLinks, localUrls(page)))
    }
    homeLinks <- localUrls(readJson(file.path(html, 'index.json')))
    for (entry in catalog$datasets) {
        expectedHash <- hash(file.path(source, entry$file))
        candidates <- homeLinks[grepl('[.](rds|csv)$', homeLinks)]
        matches <- candidates[vapply(candidates, function(url)
            file.exists(resolve(url)) && identical(hash(resolve(url)), expectedHash), TRUE)]
        if (length(matches) != 1L) stop('Missing or duplicate data download: ', entry$file)
        inputs <- inputs + 1L
    }
    localLinks <- unique(c(localLinks, homeLinks))
    missing <- localLinks[!vapply(localLinks, function(url) file.exists(resolve(url)), TRUE)]
    if (length(missing)) stop('Missing local assets: ', paste(missing, collapse = ', '))
    htmlPages <- list.files(html, '[.]html$', recursive = TRUE)
    stopifnot(length(htmlPages) == length(files) + 1L)
    renderedUrls <- character()
    for (file in file.path(html, htmlPages)) {
        text <- paste(readLines(file, warn = FALSE), collapse = '\n')
        attributes <- regmatches(text, gregexpr(
            '[[:space:]](?:href|src)="/[^" ]*"', text, perl = TRUE))[[1L]]
        urls <- sub('^[^=]+="([^"]*)"$', '\\1', attributes)
        urls <- urls[!startsWith(urls, '//')]
        if (nzchar(baseUrl) && any(!startsWith(urls, paste0(baseUrl, '/')) & urls != baseUrl))
            stop('Rendered HTML contains a link outside BASE_URL: ', file)
        renderedUrls <- unique(c(renderedUrls, urls))
    }
    missing <- renderedUrls[!vapply(renderedUrls, function(url) file.exists(resolve(url)), TRUE)]
    if (length(missing)) stop('Missing rendered HTML resources: ', paste(missing, collapse = ', '))
    report <- list(notebooks = length(files), html_pages = length(htmlPages),
        base_url = baseUrl,
        original_notebook_downloads = downloads, data_downloads = inputs,
        rendered_html_links_checked = length(renderedUrls),
        local_links_checked = length(localLinks), missing_local_files = missing)
    jsonlite::write_json(report, file.path(stage, 'download-check.json'),
        pretty = TRUE, auto_unbox = TRUE)
    cat('PASS: ', length(files), ' pages, ', downloads,
        ' original notebook downloads, ', inputs, ' data downloads; ',
        length(localLinks), ' local links checked.\n', sep = '')
}
main()

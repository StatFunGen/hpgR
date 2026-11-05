#!/usr/bin/env Rscript
# Run maintained notebooks with this package, verify them, then refresh the site.
# Rscript tools/sync_vignettes.R /path/to/P8149_2026_vignettes [--build-site]
main <- function() {
args <- commandArgs(trailingOnly=TRUE)
if (!length(args) || any(!args[-1L] %in% '--build-site'))
    stop('Supply the vignette directory and optionally --build-site')
site <- normalizePath(args[1L], mustWork=TRUE)
build_site <- '--build-site' %in% args
script <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())])
root <- normalizePath(file.path(dirname(script), '..'), mustWork=TRUE)
if (!requireNamespace('jsonlite', quietly=TRUE)) stop('Install jsonlite for notebook maintenance')
rscript <- file.path(R.home('bin'),'Rscript')
r <- file.path(R.home('bin'),'R')
jupyter <- Sys.which('jupyter')
if (!nzchar(jupyter)) stop('jupyter with nbconvert and an ir kernel is required')
sha_program <- Sys.which('shasum')
sha_args <- c('-a','256')
if (!nzchar(sha_program)) {sha_program <- Sys.which('sha256sum'); sha_args <- character()}
if (!nzchar(sha_program)) stop('A SHA-256 command is required')
sha <- function(path=NULL,text=NULL) {
    out <- if (is.null(path)) system2(sha_program,sha_args,input=text,stdout=TRUE) else
        system2(sha_program,c(sha_args,shQuote(path)),stdout=TRUE)
    strsplit(out[1L],'[[:space:]]+')[[1L]][1L]
}
run <- function(command, arguments, log) {
    status <- system2(command, arguments, stdout=log, stderr=log)
    if (!identical(status,0L)) {
        cat(tail(readLines(log,warn=FALSE),40),sep='\n')
        stop('Command failed; log: ',log,call.=FALSE)
    }
}
stage <- tempfile('hpgR-vignettes-', tmpdir=dirname(tempdir()))
dir.create(stage)
cat('Staged validation: ',stage,'\n',sep='')
source_paths <- list.files(file.path(root,'R'),'[.]R$',full.names=TRUE)
source_hashes <- setNames(lapply(source_paths,sha),paste0('R/',basename(source_paths)))
read_book <- function(path) jsonlite::fromJSON(path,simplifyVector=FALSE)
source_cells <- function(book) lapply(book$cells,function(c)
    list(id=c$id,type=c$cell_type,source=c$source))
files <- sort(list.files(site,'^P8149_Lecture_.*[.]ipynb$'))
if (!length(files)) stop('No maintained notebooks found')
original <- setNames(lapply(file.path(site,files),read_book),files)
inputs <- unique(unlist(lapply(original,function(b)names(b$metadata$notebook_input_sha256))))
input_hashes <- setNames(lapply(inputs,function(path) {
    resolved <- normalizePath(file.path(site,path),mustWork=TRUE)
    if (!startsWith(resolved,paste0(site,.Platform$file.sep))) stop('Input outside site: ',path)
    sha(resolved)
}),inputs)
for (path in c(files,inputs)) {
    dest <- file.path(stage,path)
    dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    if (!file.copy(file.path(site,path),dest)) stop('Could not stage ',path)
}
run(r,c('CMD','INSTALL',shQuote(root)),file.path(stage,'install.log'))
version <- read.dcf(file.path(root,'DESCRIPTION'))[1L,'Version']

# Figure files are generated only when explicitly requested by the lecturer.
old_figure_dir <- Sys.getenv('P8149_FIGURE_DIR',unset=NA_character_)
Sys.unsetenv(c('P8149_FIGURE_DIR','R_HOME'))
oldwd <- setwd(stage)
for (file in files) {
    cat('Execute ',file,'\n',sep=''); flush.console()
    run(jupyter,c('nbconvert','--to','notebook','--execute','--inplace',
        '--ExecutePreprocessor.kernel_name=ir','--ExecutePreprocessor.timeout=600',
        shQuote(file)),file.path(stage,paste0(file,'.log')))
}
run(rscript,c(shQuote(file.path(root,'tools/rebuild_hmm_widget.R')),shQuote(stage)),
    file.path(stage,'widgets.log'))
for (file in files) {
    book <- read_book(file.path(stage,file))
    if (!identical(source_cells(original[[file]]),source_cells(book)))
        stop('Execution changed notebook source: ',file)
    for (i in seq_along(book$cells)) if (!length(book$cells[[i]]$metadata))
        book$cells[[i]]$metadata <- structure(list(),names=character())
    cells <- Filter(function(c)c$cell_type=='code',book$cells)
    packed <- jsonlite::toJSON(lapply(cells,function(c)
        list(c$id,paste(c$source,collapse=''))),auto_unbox=TRUE)
    book$metadata$notebook_code_sha256 <- sha(text=packed)
    book$metadata$notebook_code_hash_format <- 'UTF-8 JSON cell-id/source pairs followed by newline'
    book$metadata$hpgR_source_sha256 <- source_hashes
    book$metadata$hpgR_version <- version
    book$metadata$hpgR_validation <- list(status='fresh_R_session_pass',version=version)
    book$metadata$notebook_source <- 'R calculations provided by hpgR'
    paths <- names(book$metadata$notebook_input_sha256)
    book$metadata$notebook_input_sha256 <- if(length(paths))
        setNames(lapply(paths,function(p)sha(file.path(stage,p))),paths) else
        structure(list(),names=character())
    jsonlite::write_json(book,file.path(stage,file),pretty=TRUE,auto_unbox=TRUE,
        null='null',digits=NA)
}
report <- file.path(stage,'notebook-verification.json')
run(rscript,c(shQuote(file.path(root,'tools/verify_notebook_contracts.R')),
    shQuote(stage),shQuote(report)),file.path(stage,'verification.log'))

# Detect concurrent edits before replacing any canonical notebook output.
for (i in seq_along(source_paths)) if (!identical(sha(source_paths[i]),source_hashes[[i]]))
    stop('Package source changed during notebook execution; rerun synchronization')
for (file in files) if (!identical(source_cells(read_book(file.path(site,file))),
                                  source_cells(original[[file]])))
    stop('Notebook source changed during execution: ',file)
for (path in inputs) if (!identical(sha(file.path(site,path)),input_hashes[[path]]))
    stop('Notebook input changed during execution: ',path)

# Render only the standalone HTML exports already maintained by the course.
html_files <- files[file.exists(file.path(site,sub('[.]ipynb$','.html',files)))]
for (file in html_files) run(jupyter,c('nbconvert','--to','html',
    '--HTMLExporter.embed_images=True',shQuote(file)),file.path(stage,paste0(file,'.html.log')))
for (path in c(files,sub('[.]ipynb$','.html',html_files),'widgets/hmm-basics.mjs'))
    if (!file.copy(file.path(stage,path),file.path(site,path),overwrite=TRUE))
        stop('Could not refresh ',path)
if (build_site) {
    setwd(site)
    book_jupyter <- file.path(site,'.pixi/envs/default/bin/jupyter')
    if (!file.exists(book_jupyter)) book_jupyter <- jupyter
    run(book_jupyter,c('book','build','--html'),file.path(stage,'site-build.log'))
}
setwd(oldwd)
if (is.na(old_figure_dir)) Sys.unsetenv('P8149_FIGURE_DIR') else
    Sys.setenv(P8149_FIGURE_DIR=old_figure_dir)
cat('PASS: ',length(files),' notebooks; ',length(html_files),
    ' standalone HTML exports; hpgR ',version,
    if(build_site)'; local website rebuilt' else '',
    '. Logs: ',stage,'\n',sep='')
}
main()

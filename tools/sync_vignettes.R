#!/usr/bin/env Rscript
# Run the course notebooks with this package and verify their saved outputs.
# Rscript tools/sync_vignettes.R [NOTEBOOK_DIRECTORY] [--build-site] [--figure-dir=PATH]
main <- function() {
args <- commandArgs(trailingOnly=TRUE)
script <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())])
root <- normalizePath(file.path(dirname(script), '..'), mustWork=TRUE)
directories <- args[!startsWith(args, '--')]
figure_flags <- args[startsWith(args, '--figure-dir=')]
if (length(directories) > 1L || length(figure_flags) > 1L ||
    any(startsWith(args, '--') & args != '--build-site' &
        !startsWith(args, '--figure-dir=')))
    stop('Supply one notebook directory, optionally --build-site and --figure-dir=PATH')
figure_dir <- if (length(figure_flags)) sub('^--figure-dir=', '', figure_flags) else NULL
if (!is.null(figure_dir) && !nzchar(figure_dir)) stop('Supply a figure output directory')
if (!is.null(figure_dir)) figure_dir <- normalizePath(figure_dir, mustWork=FALSE)
site <- normalizePath(if (length(directories)) directories[1L] else
    file.path(root, 'notebooks'), mustWork=TRUE)
build_site <- '--build-site' %in% args
if (build_site && !file.exists(file.path(site, 'myst.yml')))
    stop('A MyST site configuration is required for --build-site')
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
readBook <- function(path) jsonlite::fromJSON(path,simplifyVector=FALSE)
sourceCells <- function(book) lapply(book$cells,function(c)
    list(id=c$id,type=c$cell_type,source=paste(c$source,collapse='')))
catalog_path <- file.path(site, 'catalog.json')
catalog <- if (file.exists(catalog_path))
    jsonlite::fromJSON(catalog_path, simplifyVector=FALSE) else NULL
files <- if (is.null(catalog)) sort(list.files(site,'^P8149_Lecture_.*[.]ipynb$')) else
    sort(vapply(catalog$notebooks, `[[`, '', 'file'))
if (!length(files)) stop('No maintained notebooks found')
if (anyDuplicated(files) || any(basename(files) != files)) stop('Invalid notebook catalog')
original <- setNames(lapply(file.path(site,files),readBook),files)
inputs <- unique(unlist(lapply(original,function(b)names(b$metadata$notebook_input_sha256))))
input_hashes <- setNames(lapply(inputs,function(path) {
    resolved <- normalizePath(file.path(site,path),mustWork=TRUE)
    if (!startsWith(resolved,paste0(site,.Platform$file.sep))) stop('Input outside site: ',path)
    sha(resolved)
}),inputs)
assets <- list.files(file.path(site,'widgets'), recursive=TRUE, full.names=FALSE)
assets <- if (length(assets)) file.path('widgets', assets) else character()
staged_paths <- unique(c(files, inputs, assets,
    if (!is.null(catalog)) vapply(catalog$datasets, `[[`, '', 'file'),
    if (!is.null(catalog)) 'catalog.json'))
for (path in staged_paths) {
    dest <- file.path(stage,path)
    dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    if (!file.copy(file.path(site,path),dest)) stop('Could not stage ',path)
}
validation_library <- file.path(stage, 'library')
dir.create(validation_library)
run(r,c('CMD','INSTALL',paste0('--library=',shQuote(validation_library)),shQuote(root)),
    file.path(stage,'install.log'))
version <- read.dcf(file.path(root,'DESCRIPTION'))[1L,'Version']

# Figure files are generated only when explicitly requested by the lecturer.
variables <- c('P8149_FIGURE_DIR','R_HOME','R_LIBS','PATH','JUPYTER_PATH')
old_environment <- Sys.getenv(variables, unset=NA_character_)
oldwd <- getwd()
on.exit({
    setwd(oldwd)
    for (name in variables) {
        if (is.na(old_environment[[name]])) Sys.unsetenv(name) else
            do.call(Sys.setenv, setNames(list(old_environment[[name]]), name))
    }
}, add=TRUE)
Sys.unsetenv(c('P8149_FIGURE_DIR','R_HOME'))
if (!is.null(figure_dir)) Sys.setenv(P8149_FIGURE_DIR=figure_dir)
kernel_root <- file.path(stage, 'jupyter')
dir.create(file.path(kernel_root, 'kernels/ir'), recursive=TRUE)
jsonlite::write_json(list(argv=c(r,'--slave','--vanilla','-e',
    'options(width=120L); IRkernel::main()',
    '--args','{connection_file}'), display_name='R',language='R'),
    file.path(kernel_root,'kernels/ir/kernel.json'),auto_unbox=TRUE,pretty=TRUE)
Sys.setenv(R_LIBS=validation_library,
    PATH=paste(dirname(r),Sys.getenv('PATH'),sep=.Platform$path.sep),
    JUPYTER_PATH=paste(c(kernel_root,if (!is.na(old_environment[['JUPYTER_PATH']]))
        old_environment[['JUPYTER_PATH']]),collapse=.Platform$path.sep))
setwd(stage)
for (file in files) {
    executable <- Filter(function(c) c$cell_type == 'code' &&
        nzchar(trimws(paste(c$source,collapse=''))), original[[file]]$cells)
    if (!length(executable)) {
        cat('Check instructions ',file,'\n',sep=''); next
    }
    cat('Execute ',file,'\n',sep=''); flush.console()
    run(jupyter,c('nbconvert','--to','notebook','--execute','--inplace',
        '--ExecutePreprocessor.kernel_name=ir','--ExecutePreprocessor.timeout=600',
        shQuote(file)),file.path(stage,paste0(file,'.log')))
}
run(rscript,c(shQuote(file.path(root,'tools/rebuild_hmm_widget.R')),shQuote(stage)),
    file.path(stage,'widgets.log'))
for (file in files) {
    book <- readBook(file.path(stage,file))
    if (!identical(sourceCells(original[[file]]),sourceCells(book)))
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
    executable <- Filter(function(c) nzchar(trimws(paste(c$source,collapse=''))), cells)
    book$metadata$hpgR_validation <- list(status=if (!length(executable))
        'instructions_checked' else 'fresh_R_session_pass',
        version=version)
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
for (file in files) if (!identical(sourceCells(readBook(file.path(site,file))),
                                  sourceCells(original[[file]])))
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
    run(rscript,shQuote(file.path(root,'tools/build_notebook_site.R')),
        file.path(stage,'site-build.log'))
}
cat('PASS: ',length(files),' notebooks; ',length(html_files),
    ' standalone HTML exports; hpgR ',version,
    if(build_site)'; local website rebuilt' else '',
    '. Logs: ',stage,'\n',sep='')
}
main()

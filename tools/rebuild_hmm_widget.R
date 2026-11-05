# Render-data rebuild. Scientific posterior calculations are in hpgR.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1L) stop('Supply the course notebook directory')
site <- normalizePath(args[1L],mustWork=TRUE)
library(hpgR)
book <- jsonlite::fromJSON(file.path(site,'P8149_Lecture_2_HMM_Basics.ipynb'),
    simplifyVector=FALSE)
env <- new.env(parent=globalenv())
for (id in c('6c58bcb3','52d5307a','d4ba9bc1')) {
    cells <- Filter(function(c)identical(c$id,id),book$cells)
    stopifnot(length(cells)==1L)
    invisible(eval(parse(text=paste(cells[[1L]]$source,collapse='')),env))
}
frames <- function(data,values) lapply(values,function(value) {
    part <- data[as.numeric(as.character(data$setting))==value,]
    stopifnot(identical(as.integer(part$position),seq_along(env$observed)))
    part$posterior
})
demonstrations <- list(hidden=env$hidden_state,observed=env$observed,
    switch=list(values=env$switch_grid,posterior=frames(env$switch_frames,env$switch_grid)),
    error=list(values=env$error_grid,posterior=frames(env$error_frames,env$error_grid)))
path <- file.path(site,'widgets/hmm-basics.mjs')
text <- readLines(path,warn=FALSE)
take <- startsWith(text,'const demonstrations = ')
stopifnot(sum(take)==1L)
text[take] <- paste0('const demonstrations = ',
    jsonlite::toJSON(demonstrations,auto_unbox=TRUE,digits=NA),';')
writeLines(text,path,useBytes=TRUE)
cat('Updated HMM display frames from the notebook and hpgR.\n')

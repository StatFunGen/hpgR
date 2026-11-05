#!/usr/bin/env Rscript
# Compile the supplied tskit C kernels in a minimal count-data harness.
# Usage: Rscript tools/verify_sequence_statistics.R POPGENTOOLS [FIXTURE_TSV]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 1:2)
library(hpgR)
source <- file.path(normalizePath(args[1L]), "tskitR/src/tskit/trees.c")
manifest <- subset(sourceProvenance(), id == "tskit-site-statistics")
hash <- strsplit(system2("shasum", c("-a", "256", shQuote(source)), stdout = TRUE),
                 "[[:space:]]+")[[1L]][1L]
stopifnot(identical(hash, manifest$sha256))
lines <- readLines(source)
extract <- function(name, type = "int") {
    start <- grep(paste0("^", name, "\\("), lines)
    stopifnot(length(start) == 1L)
    end <- start + which(lines[seq.int(start, length(lines))] == "}")[1L] - 1L
    c(paste("static", type), lines[start:end])
}
header <- c(
    "#include <stdio.h>", "#include <stdlib.h>", "#include <stdbool.h>",
    "#include <assert.h>", "#include <math.h>",
    "typedef size_t tsk_size_t; typedef unsigned tsk_flags_t;",
    "typedef struct { const tsk_size_t *sample_set_sizes; } sample_count_stat_params_t;",
    "typedef struct { tsk_size_t num_samples; } tsk_treeseq_t;",
    "typedef struct { tsk_size_t n, k; } tsk_site_t;",
    "#define TSK_UNUSED(x) x", "#define GET_2D_ROW(x,K,j) ((x)+(K)*(j))",
    "#define TSK_STAT_POLARISED 1", "#define TSK_ERR_NO_MEMORY -1",
    "#define tsk_malloc malloc", "#define tsk_safe_free free",
    "#define tsk_trace_error(x) (x)", "#define tsk_bug_assert assert",
    # The original update function receives two complete-call allele weights.
    "static int get_allele_weights(const tsk_site_t *s, const double *counts,",
    " tsk_size_t K, const double *total, tsk_size_t *na, double **weights) {",
    " assert(K==2); *na=2; *weights=malloc(4*sizeof(double));",
    " (*weights)[0]=(*weights)[1]=s->n-s->k;",
    " (*weights)[2]=(*weights)[3]=s->k; return 0; }",
    "static void increment_nd_array_value(double *afs, tsk_size_t nd,",
    " const tsk_size_t *dims, const tsk_size_t *coordinate, double value) {",
    " assert(nd==1); afs[coordinate[0]]+=value; }")
main <- c(
    "int main(void) {",
    " puts(\"kind\\tn\\ta\\tb\\tc\\tvalue\");",
    " for (tsk_size_t n=2;n<=24;n++) for(tsk_size_t k=0;k<=n;k++) {",
    "  sample_count_stat_params_t p={&n}; double x=k,y=n-k,d1,d2;",
    "  diversity_summary_func(1,&x,1,&d1,&p);",
    "  diversity_summary_func(1,&y,1,&d2,&p);",
    "  printf(\"diversity\\t%zu\\t%zu\\t0\\t0\\t%.17g\\n\",n,k,d1+d2);",
    "  tsk_treeseq_t ts={n}; tsk_site_t site={n,k}; tsk_size_t dims[2]={n+1,n+1};",
    "  for(int polarised=0;polarised<=1;polarised++) {",
    "   double afs[25]={0};",
    "   tsk_treeseq_update_site_afs(&ts,&site,NULL,NULL,1,0,dims,polarised,afs);",
    "   for(tsk_size_t bin=1;bin<=(polarised?n-1:n/2);bin++)",
    "    printf(\"%s\\t%zu\\t%zu\\t%zu\\t0\\t%.17g\\n\",",
    "     polarised?\"polarised\":\"folded\",n,k,bin,afs[bin]);",
    "  }",
    " }",
    " for(tsk_size_t n=2;n<=12;n++) for(tsk_size_t a=0;a<=n;a++)",
    " for(tsk_size_t b=0;b<=n-a;b++) for(tsk_size_t c=0;c<=n-a-b;c++) {",
    "  sample_count_stat_params_t p={&n}; double state[3]={a,b,c}, d,r;",
    "  D_summary_func(1,state,1,&d,&p); r2_summary_func(1,state,1,&r,&p);",
    "  printf(\"D\\t%zu\\t%zu\\t%zu\\t%zu\\t%.17g\\n\",n,a,b,c,d);",
    "  printf(\"r2\\t%zu\\t%zu\\t%zu\\t%zu\\t%.17g\\n\",n,a,b,c,r);",
    " } return 0; }")
directory <- tempfile("hpgr-statistics-")
dir.create(directory)
cfile <- file.path(directory, "reference.c")
writeLines(c(header, extract("fold", "inline void"),
             extract("tsk_treeseq_update_site_afs"),
             extract("diversity_summary_func"), extract("D_summary_func"),
             extract("r2_summary_func"), main), cfile)
binary <- file.path(directory, "reference")
status <- system2(Sys.getenv("CC", "cc"),
    c("-std=c99", "-O2", shQuote(cfile), "-lm", "-o", shQuote(binary)))
stopifnot(status == 0L)
output <- system2(binary, stdout = TRUE)
reference <- read.delim(text = paste(output, collapse = "\n"), na.strings = c("nan", "-nan"))
actual <- vapply(seq_len(nrow(reference)), function(i) {
    r <- reference[i, ]
    switch(r$kind,
        diversity = hpgR:::tskitDiversity(r$a, r$n),
        D = hpgR:::tskitD(r$a/r$n, (r$a+r$b)/r$n, (r$a+r$c)/r$n),
        r2 = hpgR:::tskitR2(r$a/r$n, (r$a+r$b)/r$n, (r$a+r$c)/r$n),
        polarised = hpgR:::tskitSfs(r$a, r$n)[r$b],
        folded = hpgR:::tskitSfs(r$a, r$n, FALSE)[r$b])
}, numeric(1))
stopifnot(identical(is.na(actual), is.na(reference$value)))
error <- max(abs(actual - reference$value), na.rm = TRUE)
stopifnot(error < 2e-15)
if (length(args) == 2L) write.table(reference, args[2L], sep = "\t", quote = FALSE,
                                  row.names = FALSE, na = "NaN")
cat(sprintf("PASS: %d original-C statistic cases; max absolute difference %.3g\n",
            nrow(reference), error))

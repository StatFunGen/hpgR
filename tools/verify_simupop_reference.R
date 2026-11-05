# All margins and seed ranges are fixed before reading the reference values.
# Experimental units for population summaries are complete independent runs.
library(hpgR)
args <- commandArgs(trailingOnly=TRUE)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
root <- normalizePath(file.path(dirname(script), ".."), mustWork=TRUE)
source(file.path(root, "tests/testthat/helper-simupop-rng.R"))
fixture <- readRDS(file.path(root, "tests/testthat/fixtures/simupop-reference.rds"))
output <- if(length(args)) args[1L] else tempfile("simupop-reference-", fileext=".rds")
rng <- simupopTestRng(fixture$rng)
stopifnot(all(rng$uniform_equal), all(rng$integer_equal),
          all(rng$geometric_equal), all(rng$copy_equal))
native <- fixture$distribution
design <- native$results
blocks <- as.integer(unlist(design$blocks))
per_block <- design$per_block
rows <- list()
comparisons <- list()
z <- qnorm(.995)

record <- function(name, x, y, margin, variance=TRUE) {
    stopifnot(length(x)==length(y), all(is.finite(x)), all(is.finite(y)))
    n <- length(x)
    se <- sqrt(var(x)/n+var(y)/n)
    delta <- mean(y)-mean(x)
    mean_ok <- abs(delta)+z*se <= margin
    variance_delta <- var(y)-var(x)
    variance_se <- function(a) sqrt(max(0,
        (mean((a-mean(a))^4)-(length(a)-3)/(length(a)-1)*var(a)^2)/length(a)))
    vse <- sqrt(variance_se(x)^2+variance_se(y)^2)
    vmargin <- .35*mean(c(var(x),var(y)))
    variance_ok <- !variance || abs(variance_delta)+z*vse <= vmargin+1e-12
    grid <- sort(unique(c(x,y)))
    cdf_distance <- max(abs(ecdf(x)(grid)-ecdf(y)(grid)))
    dkw <- 2*sqrt(log(4000000)/(2*n))
    cdf_ok <- cdf_distance <= dkw
    rows[[length(rows)+1L]] <<- data.frame(name=name, n=n, native_mean=mean(x),
        R_mean=mean(y), mean_difference=delta, mean_ci_half=z*se, margin=margin,
        mean_ok=mean_ok, native_variance=var(x),R_variance=var(y),
        variance_ok=variance_ok,cdf_distance=cdf_distance,cdf_ok=cdf_ok)
    comparisons[[name]] <<- list(native=x,R=y)
    if(!mean_ok || !variance_ok || !cdf_ok)
        cat("CHECK FAILED",name,"mean",mean_ok,"variance",variance_ok,"CDF",cdf_ok,"\n")
}

for(j in seq_along(design$clones)) {
    d <- design$clones[[j]]
    seeds <- as.integer(unlist(d$seeds))+99991L
    initial <- matrix(c(rep(1L,d$N),rep(0L,d$N)),ncol=1)
    t0 <- proc.time()[3]
    y <- vapply(seeds,function(seed) {
        r <- simulateForward(initial,rep(d$N,d$steps+1L),seed)
        sum(tail(r$states,1L)[[1L]])
    },0)
    x <- as.numeric(unlist(d$values))
    record(paste0("copy_N",d$N,"_t",d$steps),x/(2*d$N),y/(2*d$N),.05)
    # Full finite-state law, including loss/fixation probabilities.
    m <- 2*d$N
    transition <- t(vapply(0:m,function(k) dbinom(0:m,m,k/m),numeric(m+1L)))
    probability <- replace(numeric(m+1L),d$N+1L,1)
    for(t in seq_len(d$steps)) probability <- as.vector(probability %*% transition)
    bound <- sqrt(log(4000000)/(2*length(x)))
    stopifnot(max(abs(cumsum(tabulate(x+1L,m+1L)/length(x))-cumsum(probability)))<bound,
              max(abs(cumsum(tabulate(y+1L,m+1L)/length(y))-cumsum(probability)))<bound)
    cat("R clone",d$N,d$steps,length(y),"runs",proc.time()[3]-t0,"seconds\n")
}

mutate <- getFromNamespace(".simupopMutate","hpgR")
meiosis <- getFromNamespace(".simupopMeiosis","hpgR")
sample_diploids <- getFromNamespace(".simupopSample","hpgR")
operator_checks <- list()
for(d in design$operators) {
    seed <- as.integer(d$seed)+99991L
    set.seed(seed)
    if(d$kind=="mutation") {
        h <- matrix(c(rep(0L,1000L),rep(1L,1000L)),nrow=2000L,ncol=100L)
        result <- mutate(h,d$u,d$v)
        observed <- c(forward=sum(result[h==0]),reverse=sum(1L-result[h==1]))
        reference <- c(forward=d$forward,reverse=d$reverse)
        expected <- 100000*c(d$u,d$v)
        sd <- sqrt(100000*c(d$u*(1-d$u),d$v*(1-d$v)))
        stopifnot(all(abs(observed-expected)<=6*sd+1),all(abs(reference-expected)<=6*sd+1))
        operator_checks[[length(operator_checks)+1L]] <- data.frame(kind="mutation",seed=d$seed,
            metric=c("forward","reverse"),native=reference/100000,R=observed/100000,
            expected=c(d$u,d$v),trials=100000)
    } else if(d$kind=="meiosis") {
        h <- meiosis(rbind(rep(0L,4),rep(1L,4)),rep(1L,40000),c(2L,2L),
                      list(c(.2,0),c(.35,0)))
        y <- c(sum(h[,1]),sum(h[,1]!=h[,2]),sum(h[,3]!=h[,4]),sum(h[,1]!=h[,3]))
        x <- as.numeric(unlist(d$metrics))
        expected <- c(.5,.2,.35,.5)
        stopifnot(all(abs(x/40000-expected)<.015),all(abs(y/40000-expected)<.015))
        operator_checks[[length(operator_checks)+1L]] <- data.frame(kind="meiosis",seed=d$seed,
            metric=c("homologue","r_chr1","r_chr2","assortment"),native=x/40000,
            R=y/40000,expected=expected,trials=40000)
        genotype <- h[seq(1,40000,2),1]+h[seq(2,40000,2),1]
        stopifnot(max(abs(tabulate(genotype+1,3)/20000-c(.25,.5,.25)))<.02,
                  max(abs(as.numeric(unlist(d$genotype))/20000-c(.25,.5,.25)))<.02)
    } else if(d$kind=="sampling") {
        h <- matrix(rep(c(rep(1L,24L),rep(0L,56L)),1),ncol=1)
        y <- vapply(seq_along(d$values),function(i) sum(sample_diploids(h,10)$haplotypes)/2,0)
        x <- as.numeric(unlist(d$values))
        record(paste0("diploid_sampling_seed",d$seed),x,y,.15)
        stopifnot(abs(mean(x)-3)<.1,abs(mean(y)-3)<.1,
            abs(var(x)-10*.3*.7*(40-10)/(40-1))<.15,
            abs(var(y)-10*.3*.7*(40-10)/(40-1))<.15)
    }
}
operators <- do.call(rbind,operator_checks)

metrics <- c("allele_frequency","heterozygosity","segregating_sites","adjacent_D2","genotype_heterozygosity")
margins <- c(.02,.015,.5,.0008,.02)
block_rows <- list()
for(i in seq_along(design$populations)) {
    dataset <- design$populations[[i]]
    d <- dataset$design
    sizes <- as.integer(unlist(d$sizes))
    genetic <- lapply(1:2,function(j) seq(0,by=-log1p(-2*d$r)/2,length.out=8))
    x <- do.call(rbind,lapply(dataset$values,unlist))
    y <- matrix(NA_real_,nrow(x),length(metrics))
    seeds <- as.integer(unlist(dataset$seeds))+99991L
    for(j in seq_len(nrow(y))) {
        result <- simulatePopulation(sizes,c(8L,8L),d$u,10L,geneticPositions=genetic,
            burnin=d$burnin,founderSize=d$founders,reverseMutationRate=d$v,seed=seeds[j])
        h <- result$observed$haplotypes
        p <- colMeans(h)
        a <- c(1:7,9:15)
        D <- colMeans(h[,a,drop=FALSE]*h[,a+1,drop=FALSE])-p[a]*p[a+1]
        y[j,] <- c(mean(p),mean(2*p*(1-p)),sum(p>0&p<1),mean(D^2),
            mean(h[seq(1,nrow(h),2),]!=h[seq(2,nrow(h),2),]))
        if(j%%per_block==0) cat("R population",d$name,"block",j/per_block,"\n")
    }
    for(k in seq_along(metrics)) record(paste(d$name,metrics[k],sep="/"),x[,k],y[,k],margins[k])
    generations <- d$burnin+length(sizes)-1L
    expected_frequency <- d$u/(d$u+d$v)*(1-(1-d$u-d$v)^generations)
    stopifnot(abs(mean(x[,1])-expected_frequency)<=6*sd(x[,1])/sqrt(nrow(x)),
              abs(mean(y[,1])-expected_frequency)<=6*sd(y[,1])/sqrt(nrow(y)))
    for(b in seq_along(blocks)) for(k in seq_along(metrics)) {
        keep <- seq_len(per_block)+(b-1)*per_block
        delta <- mean(y[keep,k])-mean(x[keep,k])
        se <- sqrt(var(x[keep,k])/per_block+var(y[keep,k])/per_block)
        standardized <- if(se>0)delta/se else 0
        block_rows[[length(block_rows)+1L]] <- data.frame(case=d$name,block=b,
            metric=metrics[k],native=mean(x[keep,k]),R=mean(y[keep,k]),z=standardized)
    }
}
summary <- do.call(rbind,rows)
block_summary <- do.call(rbind,block_rows)
passed <- all(summary$mean_ok & summary$variance_ok & summary$cdf_ok) && all(abs(block_summary$z)<5.5)
out <- list(native_engine=native$engine,summary=summary,operators=operators,
            block_summary=block_summary,raw=comparisons,passed=passed)
# A plausible incorrect mutation timing has the same mean but a wrong variance.
timing_rows <- list()
for(d in fixture$timing$results) {
    n <- 2*d$N
    x <- as.numeric(unlist(d$values))/n
    seeds <- as.integer(unlist(d$seeds))+99991L
    y <- vapply(seeds,function(seed)
        mean(simulatePopulation(c(d$N,d$N),1L,d$u,d$N,
             reverseMutationRate=0,seed=seed)$observed$haplotypes),0)
    probability <- vapply(0:n,function(k)
        sum(dbinom(0:n,n,d$u)*dbinom(k,n,(0:n)/n)),0)
    expected_variance <- d$u*(1-d$u)/n*(2-1/n)
    wrong_variance <- d$u*(1-d$u)/n
    bound <- sqrt(log(4000000)/(2*length(x)))
    for(a in list(x,y)) {
        empirical <- cumsum(tabulate(round(a*n)+1L,n+1L)/length(a))
        stopifnot(max(abs(empirical-cumsum(probability)))<bound)
        se_var <- sqrt(max(0,(mean((a-mean(a))^4)-
            (length(a)-3)/(length(a)-1)*var(a)^2)/length(a)))
        stopifnot(abs(var(a)-expected_variance)<6*se_var,
                  abs(var(a)-wrong_variance)>6*se_var)
    }
    timing_rows[[length(timing_rows)+1L]] <- data.frame(N=d$N,n=length(x),
        native_mean=mean(x),R_mean=mean(y),expected_mean=d$u,
        native_variance=var(x),R_variance=var(y),
        expected_variance=expected_variance,wrong_timing_variance=wrong_variance)
}
out$rng <- rng
out$timing <- do.call(rbind,timing_rows)
out$reference <- fixture$reference
out$R_version <- R.version.string
out$package_version <- as.character(packageVersion("hpgR"))
saveRDS(out,output)
write.csv(summary,paste0(output,".csv"),row.names=FALSE)
print(summary,row.names=FALSE,digits=5)
cat("MAX BLOCK Z",max(abs(block_summary$z)),"\n")
stopifnot(passed)
cat("PASS multi-seed native/R equivalence margins, variances, CDF and exact model laws\n")

print(rng,row.names=FALSE)
print(out$timing,row.names=FALSE)
cat("PASS mutation-timing negative control and source-aligned RNG/copying\n")
cat("Report",output,"\n")


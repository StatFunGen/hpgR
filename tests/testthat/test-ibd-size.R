test_that("fixed-word segment matching agrees with the original C++", {
    fixtures <- readRDS(test_path("fixtures", "ibd-match-reference.rds"))
    for (x in fixtures) {
        z <- detectIbdSegments(x$haplotypes, x$positions, wordSize=x$wordSize,
                               gap=x$gap, minimumCm=2)
        z <- z[do.call(order,z[c("hap1","hap2","fromSite","toSite")]),]
        rownames(z) <- NULL
        expect_identical(unname(as.matrix(z[c("hap1","hap2","fromSite","toSite")])),
                         unname(as.matrix(x$expected[1:4])))
        expect_equal(z$lengthCm,x$expected$length_cm,tolerance=3e-6)
    }
})

test_that("IBD length expectations agree with the original C++", {
    fixtures <- readRDS(test_path("fixtures", "ibd-expect-reference.rds"))
    for (x in fixtures)
        expect_equal(expectedIbdLength(x$effectiveSizes,x$minimumCm,
            x$chromosomeLengthsCm,x$haplotypePairs),x$expected,tolerance=1e-13)
    expect_equal(sum(expectedIbdLength(rep(20,10000),0,c(20,10))),30,tolerance=1e-12)
    expect_true(all(expectedIbdLength(rep(20,100),21,20)==0))
    expect_error(expectedIbdLength(0,8,20),"effectiveSizes")
})

test_that("segment fractions use all compared pairs and analyzed chromosomes", {
    x <- data.frame(lengthCm=c(8,9,10,14,20))
    z <- summarizeIbdSharing(x,10,c(20,20))
    expect_equal(z$count,c(2L,1L,2L))
    expect_equal(z$fraction,c(17,10,34)/400)
    doubled <- summarizeIbdSharing(x,20,c(20,20))
    expect_equal(doubled$fraction,z$fraction/2)
    empty <- estimateIbdSize(data.frame(lengthCm=numeric()),10,20)
    expect_true(is.na(empty$Ne))
    expect_error(summarizeIbdSharing(x,10,20,c(8,10,9,Inf)),"bins")
    expect_error(summarizeIbdSharing(x,10,20,c(8,Inf,Inf)),"bins")
})

test_that("constant-size fitting recovers a source-generated sharing expectation", {
    n <- 200; breaks <- c(8,10,14,Inf); pairs <- 10000
    tails <- vapply(breaks,function(c)if(is.finite(c))
        sum(expectedIbdLength(rep(n,1000),c,20)) else 0,0)
    total <- (head(tails,-1)-tail(tails,-1))*pairs
    # Represent each predicted bin total as equal-length segments within that bin.
    counts <- as.integer(ceiling(total/c(9,12,17)))
    segments <- data.frame(lengthCm=unlist(lapply(1:3,function(i)
        rep(total[i]/counts[i],counts[i]))))
    fit <- estimateIbdSize(segments,pairs,20)
    expect_equal(fit$Ne,n,tolerance=1e-4)
    expect_lt(fit$loss,1e-9)
})


test_that("IBD inputs reject invalid lengths before calculation", {
    segments <- data.frame(lengthCm = c(8, 10, 15))
    expect_error(summarizeIbdSharing(data.frame(lengthCm = "10"), 10, 20),
                 "lengthCm.*numeric")
    expect_error(summarizeIbdSharing(data.frame(lengthCm = NA_real_), 10, 20),
                 "lengthCm.*finite")
    expect_error(summarizeIbdSharing(segments, 10, 20, c("8", "Inf")), "bins")
    expect_error(summarizeIbdSharing(segments, 10, 0), "positive")
    expect_error(summarizeIbdSharing(segments, 10, 20, c(8, 10, -Inf)), "bins")
    expect_error(estimateIbdSize(segments, 10, 20, maxGeneration = NA_real_),
                 "maxGeneration")
})

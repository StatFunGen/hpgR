test_that("PCA projection reuses fitted frequencies and loadings", {
    G <- rbind(c(0,0,1,0), c(1,1,2,0), c(2,2,0,0), c(0,1,1,0))
    colnames(G) <- paste0("v", 1:4)
    pc <- genotypePCA(G, 2L)
    expect_identical(unname(pc$keep), c(TRUE,TRUE,TRUE,FALSE))
    expect_equal(projectGenotypePCA(G, pc), pc$scores, tolerance=1e-12)
    expected <- standardizeGenotypes(G[,pc$keep], pc$allele_frequency[pc$keep]) %*% pc$rotation
    expect_identical(projectGenotypePCA(G, pc), expected)
    expect_equal(projectGenotypePCA(G[1,,drop=FALSE], pc), pc$scores[1,,drop=FALSE])
    expect_error(projectGenotypePCA(G[,4:1], pc), "marker names")
    expect_error(projectGenotypePCA(G[,-1], pc), "order")
})

test_that("fixed-panel ancestry maximizes the existing binomial likelihood", {
    F <- cbind(c(.1,.3,.2,.4), c(.8,.6,.7,.9))
    G <- rbind(c(0,1,0,2), rep(0,4), rep(2,4))
    q <- fitSourceAncestry(G, F)
    expect_equal(q[2:3], c(0,1))
    grid <- seq(0,1,.001)
    ll <- vapply(grid, function(x) ancestryLikelihood(c(1-x,x),G[1,],F,TRUE), 0)
    expect_gte(ancestryLikelihood(c(1-q[1],q[1]),G[1,],F,TRUE), max(ll)-1e-10)
    expect_lt(abs(q[1]-grid[which.max(ll)]), .001)
    expect_error(fitSourceAncestry(G, cbind(F[,1],F[,1])), "identify")
    expect_error(fitSourceAncestry(G, F[-1,]), "per marker")
    expect_error(fitSourceAncestry(matrix(2L,1,1), matrix(0,1,2)), "differ")
    expect_error(fitSourceAncestry(matrix(c(2L,1L),1), rbind(c(0,0),c(.1,.9))), "zero probability")
})

test_that("ancestry factors reproduce exact path enumeration", {
    F <- cbind(A=c(.2,.2), B=c(.8,.8))
    fit <- ancestryHMM(c(1,1), F, c(0,1), -log(.8))
    expect_equal(unname(fit$posterior[2,2]),148/161,tolerance=1e-12)
    expect_equal(unname(ancestryHMM(1,F[1,,drop=FALSE],0)$posterior[1,2]), .8)
    expect_equal(ancestryHMM(c(NA,NA),F,c(0,1),initial=c(.3,.7))$posterior,
                 matrix(c(.3,.3,.7,.7),2,dimnames=list(NULL,c("A","B"))))
    expect_equal(ancestryHMM(c(1,0),F,c(0,1),0)$posterior[1,],
                 ancestryHMM(c(1,0),F,c(0,1),0)$posterior[2,])
    expect_error(ancestryHMM(c(0,1),F,c(1,0)), "nondecreasing")
    expect_error(ancestryHMM(c(0,2),F,c(0,1)), "binary")
    expect_error(ancestryHMM(c(0,1),F,c(0,1),initial=c(.4,.4)), "sum")
    expect_error(ancestryHMM(c(1,1),F,c(0,1),-1), "nonnegative")
})

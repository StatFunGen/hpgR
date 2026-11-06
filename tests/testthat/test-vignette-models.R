courseFixture <- function() readRDS(test_path('fixtures', 'vignette-calculations.rds'))

test_that('Gaussian EM retains the fixed sample, soft counts and likelihood trace', {
    expected <- courseFixture()[['2_EM_Basics']]
    set.seed(8149)
    x <- simulateGaussianMixture(80L, c(-1.5, 1.5))
    expect_equal(x$observed, expected$observed, tolerance=1e-12)
    expect_identical(x$group, expected$hidden_group)
    fit <- fitGaussianMixture(x$observed, list(pi=.5, mu=c(-.5,.5)))
    expect_equal(fit[names(expected$fit)], expected$fit, tolerance=1e-12)
    expect_true(fit$converged)
    expect_equal(fit$iterations, nrow(fit$history) - 1L)
    expect_true(all(diff(fit$history$loglik) >= -1e-10))
    r <- gaussianResponsibilities(c(-100,0,100), .5, c(-1,1))
    expect_equal(r[2], .5)
    expect_true(all(r >= 0 & r <= 1))
})

test_that('the genotype error event model and EM agree with the fixed reference', {
    expected <- courseFixture()[['2_MAF_MLE_EM']]
    set.seed(581)
    x <- simulateGenotypeCalls(10000,.3,.05)
    expect_identical(x$truth, expected$g_true)
    expect_identical(x$observed, expected$d_observed)
    expect_identical(x$error_event, expected$error_event)
    fit <- fitGenotypeError(x$observed)
    expect_equal(fit[names(expected$em_fit)], expected$em_fit, tolerance=1e-12)
    expect_true(fit$converged)
    expect_true(all(diff(fit$history$loglik) >= -1e-8))
    expect_equal(unname(rowSums(genotypeErrorPosterior(0:2,.3,.05))), rep(1,3))
    expect_equal(genotypeErrorMatrix(0), diag(3), ignore_attr=TRUE)
    expect_equal(genotypeErrorMatrix(1), matrix(1/3,3,3), ignore_attr=TRUE)
    g <- c(0,2,1,0,0)
    ml <- optimize(genotypeLogLik, c(1e-6,1-1e-6), maximum=TRUE, g=g)
    expect_equal(ml$maximum, mean(g)/2, tolerance=1e-4)
})

test_that('HMM messages reproduce original C++ probability-space recurrences', {
    fixture <- readRDS(test_path('fixtures','probability-hmm.rds'))
    expect_lte(fixture$max_absolute_error,1e-12)
    for (case in fixture$cases) {
        actual <- hmmPosterior(case$emission,case$transition,case$initial)
        expect_equal(actual[names(case$expected)],case$expected,tolerance=1e-12)
    }
})

test_that('HMM marginals and likelihood agree with exhaustive hidden paths', {
    emission <- rbind(c(.8,.2),c(.3,.7),c(1,1),c(.9,.1))
    tr <- matrix(c(.9,.1,.2,.8),2,byrow=TRUE)
    initial <- c(.6,.4)
    paths <- as.matrix(expand.grid(rep(list(1:2),4)))
    weight <- apply(paths,1,function(s) initial[s[1]] *
        prod(emission[cbind(1:4,s)]) * prod(tr[cbind(head(s,-1),tail(s,-1))]))
    exact <- sapply(1:2,function(k) colSums((paths==k)*weight)/sum(weight))
    fit <- hmmPosterior(emission,tr,initial)
    expect_equal(fit$posterior,unname(exact),tolerance=1e-12)
    expect_equal(fit$log_likelihood,log(sum(weight)),tolerance=1e-12)
    expect_equal(hmmPosterior(matrix(c(.2,.8),1),tr,initial)$posterior,
        matrix(initial*c(.2,.8)/sum(initial*c(.2,.8)),1))
    varying <- array(rep(tr,3),c(2,2,3))
    expect_equal(hmmPosterior(emission,varying,initial),fit)
    expect_equal(hmmPosterior(matrix(1,3,1),matrix(1,1,1))$posterior,matrix(1,3,1))
    expect_equal(hmmPosterior(matrix(1,3,1),array(1,c(1,1,2)))$posterior,matrix(1,3,1))
    expect_error(hmmPosterior(matrix(0,2,2),tr), 'zero probability')
    expect_error(hmmPosterior(emission,tr*2), 'sum to one')
    expect_error(hmmPosterior(emission,tr,c(-1,2)), 'initial')
})

test_that('copying simulation and imputation retain the course example', {
    expected <- courseFixture()[['2_Li_Stephens_HMM']]
    set.seed(534)
    frequency <- runif(100,.2,.8)
    H <- replicate(20,simulateReferenceHaplotype(frequency))
    expect_equal(unname(H),unname(expected$H))
    dimnames(H) <- dimnames(expected$H)
    target <- simulateCopyingPath(H,.05,.02)
    typed <- sort(sample.int(100,30))
    observed <- rep(NA_integer_,100)
    observed[typed] <- target$haplotype[typed]
    expect_identical(target$path,expected$true_path)
    expect_identical(target$haplotype,expected$true_target)
    expect_identical(observed,expected$observed)
    tr <- copyingTransition(.95,rep(1/20,20))
    e <- copyingEmission(H,observed,.02)
    fit <- hmmPosterior(e,tr)
    expect_equal(fit$posterior,expected$fit$posterior,tolerance=1e-12)
    expect_equal(fit$log_likelihood,expected$fit$log_likelihood,tolerance=1e-12)
    expect_true(all(e[is.na(observed),]==1))
    expect_equal(rowSums(tr),rep(1,20))
    expect_equal(copyingTransition(1,c(.2,.8)),diag(2))
})

test_that('linkage and kinship retain the same probability calculations', {
    fixture <- courseFixture()
    expected <- fixture[['3_Linkage_LOD']]$comparison
    actual <- rbind(fitLinkage(10,0),fitLinkage(9,1))
    expect_equal(actual,expected,tolerance=1e-12)
    expect_equal(linkageLOD(.5,10,0),0)
    expect_equal(linkageLOD(0,10,0),10*log10(2))
    x <- siblingSharingTest(c(2L,8L,10L))
    expected <- fixture[['3_Affected_Sibling_Linkage']]
    expect_equal(unname(x),unname(unlist(expected)),tolerance=1e-12)
    expected <- fixture[['3_Genome_Wide_Kinship']]
    set.seed(8149)
    f <- runif(5000,.1,.5)
    u <- t(replicate(100,compareKinship(simulateUnrelatedPair(f),f)))
    s <- t(replicate(100,compareKinship(simulateSiblingPair(f),f)))
    expect_equal(f,expected$allele_frequency,tolerance=1e-12)
    expect_equal(u,expected$unrelated,tolerance=1e-12)
    expect_equal(s,expected$siblings,tolerance=1e-12)
    expect_equal(u[,'PLINK'],u[,'frequency_moment'],tolerance=1e-12)
    expect_equal(s[,'PLINK'],s[,'frequency_moment'],tolerance=1e-12)
})

test_that('family traits preserve conditional covariance and variance estimates', {
    expected <- courseFixture()[['4_Family_Heritability']]
    set.seed(8149)
    x <- simulateSiblingTraits(8000,25,4,16,170)
    for (key in c('phenotype_1','phenotype_2','ibd_sharing'))
        expect_equal(x[[key]][1:20],expected[[key]],tolerance=1e-12)
    fit <- lm(I((x$phenotype_1-x$phenotype_2)^2) ~ I(1-x$ibd_sharing))
    expect_equal(unname(coef(fit)[2]/2),expected$V_A_hat)
    expect_equal(unname(coef(fit)[1]/2),expected$V_E_hat)
    expect_equal(unname(coef(fit)[2]/2)/var(c(x$phenotype_1,x$phenotype_2)),
        expected$h2_hat)
})

test_that('SNP likelihood and ancestry EM retain their numerical answers', {
    set.seed(8149)
    f <- runif(3000,.1,.5)
    G <- simulateGenotypes(350,f)
    Z <- standardizeGenotypes(G,f)
    K <- geneticRelationship(Z)
    trait <- simulatePolygenicTrait(Z,.5,.5)
    eig <- eigen(K,symmetric=TRUE)
    y <- as.vector(crossprod(eig$vectors,trait$phenotype))
    fit <- optimize(geneticProfileNLL,c(.001,.999),
        eigenvalues=pmax(eig$values,1e-8),rotated_phenotype=y)
    # Eigenvectors and the numerical minimizer vary slightly across BLAS builds.
    expect_equal(fit$minimum,courseFixture()[['4_SNP_Heritability']]$h2_likelihood,
        tolerance=sqrt(.Machine$double.eps))
    expect_equal(var(trait$genetic_value),.5)
    expected <- courseFixture()[['5_Homework_5_Starter']]
    pc <- genotypePCA(expected$toy_x)
    expect_equal(pc[names(expected$toy_pc)],expected$toy_pc,tolerance=1e-12)
    expect_equal(projectGenotypePCA(expected$toy_x, pc), pc$scores, tolerance=1e-12)
    set.seed(5)
    ancestry <- fitAdmixture(expected$toy_x,2,pc$scores)
    expect_equal(ancestry,expected$toy_fit,tolerance=1e-12)
    expect_true(ancestry$converged)
    expect_equal(rowSums(ancestry$q),rep(1,nrow(expected$toy_x)))
    expect_true(all(diff(ancestry$log_likelihood)>-1e-6))
    frequency <- cbind(c(.3,.8,.6,.6,.2),c(.6,.5,.3,.8,.6))
    genotype <- c(0,2,1,2,0)
    expect_equal(ancestryLikelihood(c(1,0),genotype,frequency),.03468165,
        tolerance=1e-7)
    expect_equal(log(ancestryLikelihood(c(.4,.6),genotype,frequency)),
        ancestryLikelihood(c(.4,.6),genotype,frequency,log=TRUE))
})

test_that('migration, pooling, selection and sampling uncertainty retain their laws', {
    fixture <- courseFixture()
    expected <- fixture[['10_Migration_and_Wahlund']]
    actual <- migrationFrequencies(.8,.2,.05,30)
    expect_equal(actual,expected$frequency)
    expect_equal(pairwiseFst(actual[,1],actual[,2]),expected$observed_fst)
    pool <- poolAlleleFrequencies(c(.8,.2))
    expect_equal(c(pool$H_S,pool$H_T,pool$F_ST),c(.32,.5,.36))
    expect_true(is.na(poolAlleleFrequencies(c(0,0))$F_ST))
    set.seed(1010)
    split <- replicate(100,simulateMigrationDrift(c(.5,.5),100,0,30),simplify=FALSE)
    mixed <- replicate(100,simulateMigrationDrift(c(.8,.2),100,.05,30),simplify=FALSE)
    expect_identical(split,fixture[['10_Split_Drift_Migration']]$split_pairs)
    expect_identical(mixed,fixture[['10_Split_Drift_Migration']]$migration_pairs)
    w <- c(.25,1,.92)
    eq <- balancedFrequency(w)
    expect_equal(unname(selectedFrequency(eq,w)),unname(eq))
    expect_equal(unname(selectionFrequencies(.02,w,150)),
        unname(fixture[['11_Sickle_Cell_Trait']]$with_malaria))
    expect_equal(selectedFrequency(.3,c(1,1,1)),.3)
    expect_error(balancedFrequency(c(1,1,1)))
    expect_equal(zeroCountUpperBound(c(3,10,30,100)),
        fixture[['12_Ancient_Dogs_Uncertainty']]$upper_if_zero)
    expect_equal((1-zeroCountUpperBound(3))^3,.05)
})

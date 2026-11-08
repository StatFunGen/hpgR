test_that("demographic SFS matches equilibrium and the integrated pair clock", {
    for(n in c(4L,12L)){
        expect_equal(expectedHistorySFS(n,c(100,rep(100,60)),2.5e-5,10000),
                     expectedSFS(n,100,2.5e-5,10000),tolerance=1e-9)
    }
    # Pair survival during a 20-generation bottleneck, then ancestral size100.
    meanPair <- 40 * (1 - exp(-20/40)) + exp(-20/40) * 200
    expect_equal(expectedHistorySFS(2,c(100,rep(20,20)),.001,100),
                 2 * .001 * 100 * meanPair)
    parents<-c(rep(20,20),seq(20,1000,length.out=40))
    pred<-expectedHistorySFS(4,c(100,parents),1,1)
    trees<-lapply(1:1000,function(i)pgKingmanTree(100,4,814900L+i,populationSizes=rev(parents)))
    exposure<-t(vapply(trees,function(t)pgTreeSummaries(t)$branch_generations,numeric(3)))
    expect_true(all(abs(colMeans(exposure)-pred)<5*apply(exposure,2,sd)/sqrt(1000)))
})

test_that("piecewise coalescent waiting follows the integrated pair hazard", {
    sizes<-c(rep(10,10),rep(100,20))
    time<-vapply(1:1200,function(i){
        t<-pgKingmanTree(100,2,112000L+i,populationSizes=sizes)
        pgTreeSummaries(t)$TMRCA
    },0)
    for(t in c(5,15,40)){
        hazard<-sum(head(c(sizes,rep(100,40)),t)^-1)/2
        probability<-1-exp(-hazard)
        expect_lt(abs(mean(time<=t)-probability),5*sqrt(probability*(1-probability)/length(time)))
    }
})

test_that("binary mutations reproduce the finite-site pair transition law", {
    ancestry<-structure(list(nodes=data.frame(id=0:2,time=c(0,0,10),flags=c(1L,1L,0L)),
        edges=data.frame(left=0,right=1,parent=2L,child=0:1),sampleIds=0:1,
        sequenceLength=1,timeUnits="generations"),class=c("hpgR_ancestry","list"))
    difference<-vapply(1:1200,function(i){
        out<-mutateAncestry(ancestry,.1,seed=31000L+i,model="binary")
        data<-ancestryData(out)$observed$haplotypes
        expect_true(all(out$sites$position==0))
        expect_true(all(out$mutations$derived_state %in% 0:1))
        sum(data[1,]!=data[2,])
    },0L)
    p<-(1-exp(-4*.1*10))/2
    expect_lt(abs(mean(difference)-p),5*sqrt(p*(1-p)/1200))
})

test_that("fast point trees agree with full parental records", {
    pop<-simulatePopulation(c(8L,4L,8L),c(5L,5L),.01,4L,
        geneticPositions=list(seq(0,.2,length.out=5),seq(0,.2,length.out=5)),
        burnin=200L,founderSize=4L,seed=921L,recordSites=c(1L,4L,8L))
    for(site in c(1L,4L,8L)){
        fast<-populationTree(pop,site,c(1,3,5,8))
        slow<-pgCompressTree(populationGenealogy(pop,site,c(1,3,5,8)),
            sample_ids=paste0("H",c(1,3,5,8)))
        expect_equal(pgTreeSummaries(fast),pgTreeSummaries(slow),ignore_attr=TRUE)
        expect_equal(recordedSubset(fast)$ages,recordedSubset(slow)$ages)
    }
})

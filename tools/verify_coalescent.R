# Executable checks of the model calculations used in the genealogy lecture.
# Pritchard 2.2; canonical Notes 6.7 and 8.10--8.13;
# Laruson and Reed (2021), pp.176--184. No empirical outcomes are simulated.
set.seed(8149)
Ne <- 1000
B <- 100000L
g <- rgeom(B,1/(2*Ne))+1
e <- rexp(B,1/(2*Ne))
stopifnot(abs(mean(g)/(2*Ne)-1)<0.02,
          abs(mean(e)/(2*Ne)-1)<0.02)
cat(sprintf('Geometric mean %.2f; exponential mean %.2f; expectation %.2f\n',
            mean(g),mean(e),2*Ne))
for (n in c(2,4,6,100)) {
  k <- 2:n
  waits <- 2*Ne/choose(k,2)
  stopifnot(abs(sum(waits)-4*Ne*(1-1/n))<1e-8)
  stopifnot(abs(sum(k*waits)-4*Ne*sum(1/(1:(n-1))))<1e-8)
}
H <- rbind(c(0,0,0),c(1,0,0),c(0,1,0),c(0,0,1))
n <- nrow(H)
counts <- colSums(H)
pairs <- combn(n,2)
pi.direct <- mean(apply(pairs,2,function(pair) sum(H[pair[1],]!=H[pair[2],])))
pi.counts <- sum(2*counts*(n-counts)/(n*(n-1)))
thetaW <- ncol(H)/sum(1/(1:(n-1)))
stopifnot(pi.direct==1.5,pi.counts==pi.direct,abs(thetaW-18/11)<1e-12)
cat(sprintf('Four-sequence example: S=%d; pi_L=%.2f; theta_W=%.2f\n',ncol(H),pi.direct,thetaW))
for (n in c(4,10,100)) {
  i <- 1:(n-1)
  stopifnot(abs(sum(2*i*(n-i)/(n*(n-1))/i)-1)<1e-12)
  # Two copies sampled without replacement for SNP discovery.
  retain <- i*(n-i)/choose(n,2)
  stopifnot(max(abs(retain/i-2*(n-i)/(n*(n-1))))<1e-12)
}
growth.ratio <- 1+1000*100/(4*20000)
stopifnot(growth.ratio==2.25)
cat('Instantaneous-growth singleton ratio:',growth.ratio,'\n')
cat('Pair-count, time-sum, branch-length, SFS weighting and ascertainment identities PASS\n')

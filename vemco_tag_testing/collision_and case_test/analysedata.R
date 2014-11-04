## Script for testing if two transmitters perform equally well
## 04.02.2013
rm(list=ls())

T <- 120*60 ## Time of a trial
N <- 2 ## Number of tags in test

pr <- rep(0,N) ## Ping rate of each tag
pr[1] <- 120 
pr[2] <- 120

n <- rep(0,N) ## Expected number of pings
n[1] <- T/pr[1]
n[2] <- T/pr[2]

x <- rep(0,N) ## Number of successes (received pings)
x[1] <- round(0.9*n[1])
x[2] <- round(0.8*n[2])

prop.test(x,n)


sampsize <- function(P1,P2,r,alpha,beta){
  delta <- abs(P1-P2)
  Pbar <- (P1+r*P2)/(r+1)
  Qbar <- 1-Pbar
  zalpha <- qnorm(1-alpha)
  zbeta <- qnorm(1-beta)
  mprime <- (zalpha*sqrt((r+1)*Pbar*Qbar) + zbeta*sqrt(r*P1*(1-P1)+P2*(1-P2)))^2/(r*delta^2)
  m <- mprime/4*(1+sqrt(1+ (2*(r+1))/(r*mprime*delta)))^2
  N <- (r+1)*m
  m
}

## Sample size calculator (http://www.cct.cuhk.edu.hk/stat/proportion/Casagrande.htm)
P1 <- c(0.95,0.9,0.8,0.7,0.6,0.5)
P2 <- seq(0.95,0.05,length=100)
r <- 1 ## r = 1 implies that the two sample sizes are identical
## Significance level (probability of rejecting the true null hypothesis)
alpha <- 0.05/2 ## Divide by 2 to get two sided results (P1=P2 versus P1 neq P2)
## Power of the test (probability of not rejecting the false null hypothesis)
beta <- 0.05

np1 <- length(P1)
np2 <- length(P2)
m <- matrix(0,np1,np2)
for(i in 1:np1){
  for(j in 1:np2){
    m[i,j] <- sampsize(P1[i],P2[j],r,alpha,beta)
  }
}

##pdf('samplesize.pdf')
plot(P2,m[1,],typ='l',ylim=c(0,360),xlim=c(0,1),ylab='Number of sent pings required',xlab='Tag 2 detection probability',main=paste('alpha:',alpha,'beta:',beta))
for(i in 2:np1){
  inds <- P2<P1[i]
  lines(P2[inds],m[i,inds],col=i)
}
legend('topleft',legend=P1,col=1:np1,lty=1)
##dev.off()

## Script for formatting acoustic data from Diamond head ranging experiment
## 18.01.2012
rm(list=ls())
library(fields)
library(Hmisc)
source('/home/mwp/work/R/basicfuns.R')

## Meta data
latr <- 21 + 14.154/60  # Lat receivers
lonr <- 157 + 48.881/60 # Lon receivers
latrb <- 21 + 14.174/60
lonrb <- 157 + 48.86/60
latt1 <- 21 + 14.214/60
lont1 <- 157 + 48.973/60
latt1b <- 21 + 14.186/60
lont1b <- 157 + 49.047/60
latt2 <- 21 + 14.317/60
lont2 <- 157 + 49.175/60

x <- matrix(0,5,2)
x[1,] <- c(lonr,latr)
x[2,] <- c(lonrb,latrb)
x[3,] <- c(lont1,latt1)
x[4,] <- c(lont1b,latt1b)
x[5,] <- c(lont2,latt2)
distmat <- rdist.earth(x,x) ## Distance matrix


## Read receiver data
filename <- 'diamondheadranging.csv'
datraw <- read.csv(filename)
timeraw <- as.POSIXct(strptime(as.character(datraw[,1]),format='%Y-%m-%d %H:%M:%S'),tz='UTC')

## Plot pings in a time interval
timestart <- as.POSIXct('2012-12-10 00:00:01 HST') ## Begin time for full data
##timestart <- as.POSIXct('2011-01-13 00:00:01 HST')
##timeend <- as.POSIXct('2011-01-13 20:00:01 HST')
timeend <- as.POSIXct('2012-12-23 00:00:01 HST') ## End time for full data
inds <- which(timeraw >= timestart & timeraw < timeend)
time <- as.POSIXlt(timeraw[inds],tz='HST')
dat <- datraw[inds,]
inds <- 1:length(time) ## All

## Receiver names
recs <- levels(dat$Receiver)
recnam <- c('60m','30m','3m')
## Transmitter names
tags <- levels(dat$Transmitter)
reftags <- tags[7:12]
tagnam <- c('1 3m','1 30m','1 60m','2 3m','2 30m','2 60m')

tagdata <- list()
for(i in 1:6){
  tagdata[[i]] <- list()
  for(j in 1:3){
    inds <- dat$Transmitter==reftags[i] & dat$Receiver==recs[j]
    tagdata[[i]][[j]] <- time[inds]
    ##plot(tagdata[[i]][[j]],rep(1,length(tagdata[[i]][[j]])))
  }
}

dt <- 3*60*60 ## Window width
pr <- 180 ## Average ping rate (between 110 and 250)
expected <- dt/pr
numintv <- as.numeric(difftime(timeend,timestart,units='secs'))/dt
windat <- array(0,dim=c(6,3,numintv))
windatmc <- array(0,dim=c(6,3,numintv))
means <- matrix(0,3,6)
colnames(means) <- c('t1 3 m','t1 30 m','t1 60 m','t2 3 m','t2 30 m','t2 60 m')
rownames(means) <- c('r 60 m','r 30 m','r 3 m')
stds <- means

ping <- c() ## Anova data
tagh <- c()
rech <- c()
loc <- c()
depths <- as.factor(c('3m','30m','60m','3m','30m','60m'))

do.plot <- TRUE
if(do.plot){
  graphics.off()
  dev.new(height=8,width=12)
  par(mfrow=c(3,6))
}
for(i in 1:6){
  for(j in 1:3){
    tt <- tagdata[[i]][[j]]
    for(k in 1:numintv){
      ts <- timestart + (k-1)*dt
      te <- timestart + k*dt
      inds <- which(tt >= ts & tt < te)
      windat[i,j,k] <- length(inds)/expected
      ##windat[i,j,k] <- length(inds)
    }
    windatmc[i,j,] <- windat[i,j,]-means[j,i]
    ## Logit transform data
    da <- logit(windat[i,j,seq(1,numintv,by=1)])
    da <- da[is.finite(da)]
    ## Calculate mean and standard deviation
    means[j,i] <- mean(da)
    stds[j,i] <- sqrt(var(da)/numintv)
    ## Do different check to see if data is uncorrelated (acf), approximately Gaussian (hist, qqnorm)
    if(do.plot) hist(da)
    ## Collect data in a data frame for easier anova'ing
    ping <- c(ping,da)
    tagh <- c(tagh,rep(i,length(da)))
    rech <- c(rech,rep(j,length(da)))
  }
}

pingdata <- as.data.frame(list(ping=ping,tagh=as.factor(tagh),rech=as.factor(rech)))

## mod <- lm(ping~tagh+rech,data=pingdata) ## Anova is not appropriate since the observations don't have equal variance

## Constant receiver height, variable transmitter height
aa <- c(1,1,2)
bb <- c(2,3,3)
pvals <- matrix(0,3,6)
colnames(pvals) <- c('3,30','3,60','30,60','3,30','3,60','30,60')
rownames(pvals) <- c('r 60 m','r 30 m','r 3 m')
pvals2 <- pvals
intv2 <- pvals
est2 <- pvals
for(j in 1:3){
  for(i in 1:3){
    ## Assuming data are normal
    mu <- abs(means[j,aa[i]]-means[j,bb[i]])
    std <- sqrt(stds[j,aa[i]]^2+stds[j,bb[i]]^2)
    pvals[j,i] <- pnorm(0,mu,std)
    mu <- abs(means[j,aa[i]+3]-means[j,bb[i]+3])
    std <- sqrt(stds[j,aa[i]+3]^2+stds[j,bb[i]+3]^2)
    pvals[j,i+3] <- pnorm(0,mu,std)
    ## Non-parametrix
    w1 <- wilcox.test(windat[aa[i],j,],windat[bb[i],j,],alternative='t',conf.int=TRUE)
    w2 <- wilcox.test(windat[aa[i]+3,j,],windat[bb[i]+3,j,],alternative='t',conf.int=TRUE)
    pvals2[j,i] <- w1$p.value
    pvals2[j,i+3] <- w2$p.value
    intv2[j,i] <- abs(diff(w1$conf.int))
    intv2[j,i+3] <- abs(diff(w2$conf.int))
    est2[j,i] <- w1$estimate
    est2[j,i+3] <- w2$estimate
  }
}
## Which means are different (TRUE means they are different)
pvals<0.025 ## 0.025 since it is two sided
pvals2<0.025 ## 0.025 since it is two sided

foo <- latex(round(pvals,digits=3),label='tab:pvals',caption='Testing for difference in means between transmitters (normality assumed).',file='pvals.tex')
foo <- latex(round(pvals2,digits=3),label='tab:pvals2',caption='Testing for difference in means between transmitters (non-parametric, Mann-whitney test).',file='pvals2.tex')

## Script for formatting acoustic data from Diamond head ranging experiment
## 18.01.2012
rm(list=ls())
library(fields)
library(Hmisc)
##library(mixtools) ## For contour ellipse
##source('../../simulation/ssmpaper/analysis/simfuns.R')
##source('../../acousticfuns.R')

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

dt <- 60*60 ## Window width
pr <- 180 ## Average ping rate (between 110 and 250)
expected <- dt/pr
numintv <- as.numeric(difftime(timeend,timestart,units='secs'))/dt
windat <- array(0,dim=c(6,3,numintv))
windatmc <- array(0,dim=c(6,3,numintv))
means <- matrix(0,6,3)

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
    means[i,j] <- mean(windat[i,j,])
    windatmc[i,j,] <- windat[i,j,]-means[i,j]
  }
}


## --- Fit a simple linear model
p <- 1
nt <- dim(windatmc)[3]
wd <- windat
I <- 3
nfit <- 200
aa <- c(1,1,2)
bb <- c(2,3,3)
coeff <- matrix(0,6,3)
colnames(coeff) <- c('30~60','3~60','3~30')
rownames(coeff) <- c('t1 3 m','t1 30 m','t1 60 m','t2 3 m','t2 30 m','t2 60 m')
coeff <- t(coeff)
rsquare <- coeff
resse <- coeff

for(I in 1:6){
  Y <- cbind(wd[I,1,(1+p):nt],wd[I,2,(1+p):nt],wd[I,3,(1+p):nt])
  X <- cbind(wd[I,1,1:(nt-p)],wd[I,2,1:(nt-p)],wd[I,3,1:(nt-p)]) 
  for(j in 1:3){
    a <- aa[j]
    b <- bb[j]
    x <- Y[1:nfit,a]
    a1 <- mean(x)
    x <- x-a1
    y <- Y[1:nfit,b]
    a2 <- mean(y)
    y <- y-a2
    theta <- solve(t(x)%*%x)%*%t(x)%*%y ## Normal equations
    sigma <- (y-theta*x) %*% (y-theta*x)/(nt-1)
    lmres <- lm(y~x)
    coeff[j,I] <- lmres$coefficients[2]
    rsquare[j,I] <- summary(lmres)$r.squared
    resse[j,I] <- summary(lmres)$sigma
  }
}

coeff <- round(coeff,digits=3)
rsquare <- round(rsquare,digits=3)
resse <- round(resse,digits=3)
foo <- latex(coeff,label='tab:coeff',caption='Values of $\\theta$ for fit of different transmitter combinations.')
foo <- latex(rsquare,label='tab:rsquare',caption='Values of $R^2$ for fit of different transmitter combinations.')
foo <- latex(resse,label='tab:resse',caption='Values of residual standard error for fit of different transmitter combinations.')



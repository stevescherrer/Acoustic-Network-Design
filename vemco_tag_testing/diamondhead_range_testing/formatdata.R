## Script for formatting acoustic data from Diamond head ranging experiment
## 18.01.2012
rm(list=ls())
library(fields)
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

##pdf('stationplot.pdf')
##plot(lonr,latr,xlim=c(157.813,157.82),ylim=c(21.235,21.239))
##points(lonrb,latrb,pch=2)
##points(lont1,latt1,col=2)
##points(lont1b,latt1b,col=2,pch=2)
##points(lont2,latt2,col=3)
##legend('topleft',legend=c('drop','ranging'),pch=c(1,2))
##dev.off()

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
##inds <- which(dat$Station.Name=='Palmyra_AA')
##plot(time[inds],dat$Station.Name[inds])
##par(las=1)
##plot(time[inds],dat$Station.Name[inds],yaxt='n',ylab='',main=paste('Tag ID: ',substr(as.character(dat$Transmitter[1]),10,20),', # data: ',length(time),sep=''),xlab='')
##lvls <- substr(levels(dat$Station.Name[1]),9,20)
##axis(2,at=1:length(lvls),labels=lvls)

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
  }
}

graphics.off()
##pdf('detectionpercentage.pdf',height=7,width=12)
##pdf('histdetectperc.pdf',height=6,width=10)
##pdf('density.pdf',height=6,width=10)
dev.new(height=7,width=12)
par(mfrow=c(3,6))
corrmat <- matrix(0,18,18)
ci <- 1
dists <- c('121m','368m') ## Distances from receivers
for(j in 1:3){ for(i in 1:6){
  ##plot(density(as.integer(tagdata[[i]][[j]]-tagdata[[i]][[j]][1])),xlab='',ylab='',main=paste('i:',i,'j:',j))
  mn <- ''
  yl=''
  ##if(j==1) mn <- tagnam[i]
  if(j==1){if(i<4){ mn <- paste(dists[1],',  ',tagnam[i],sep='') }else{ mn <- paste(dists[2],',  ',tagnam[i],sep='') }}
  if(i==1) yl <- recnam[j]
  xlb <- paste('[',i,',',j,'] ',ci,sep='')
  plot(windat[i,j,],typ='l',ylim=c(0,1.2),main=mn,ylab=yl,xlab=xlb)  ## Detection percentage
  abline(h=1,lty=2,col=2)
  ##hist(windat[i,j,],typ='l',xlim=c(0,1.2),main=mn,ylab=yl,xlab=xlb)  ## Histogram of detection percentage
  ##abline(v=1,lty=2,col=2)
  ##plot(density(windat[i,j,]),xlim=c(0,1),main=mn,ylab=yl) ## Density plot
  ##acf(windat[i,j,],lag.max=25) ## Autocorrelation function
  ##acf(windat[i,j,seq(1,dim(windat)[3],by=2)],lag.max=25) ## Autocorrelation function
  cj <- 1
  for(jj in 1:3){ for(ii in 1:6){
    corrmat[ci,cj] <- max(abs(ccf(windat[i,j,],windat[ii,jj,],plot=FALSE)$acf))
    cj <- cj+1
  }}
  ci <- ci+1
}}

dev.off()

## Plot cross correlation between data sets
X11()
image(1:18,1:18,corrmat[18:1,],main='Maximum cross-correlation between data sets')
grid()
## It seems from the plot that data coming from the same transmitter are correlated, whereas receiver seems to have little effect.


## Assuming that there is not temporal variability in the data we can apply standard statistical models such as GLM
ping <- as.vector(apply(windat,c(1,2),sum))
noping <- 312*30-ping
tagh <- as.factor(rep(c(3,30,60),6))
rech <- as.factor(c(rep(60,6),rep(30,6),rep(3,6)))
loc <- as.factor(rep(c(rep('A',3),rep('B',3)),3))

pingdata <- as.data.frame(list(tagh=tagh,rech=rech,loc=loc,ping=ping,noping=noping))
n <- pingdata$ping+pingdata$noping
pingdata$prop <- pingdata$ping/n
pdA <- pingdata[pingdata$loc=='A',]

mod <- glm(cbind(ping,noping)~tagh+rech+loc,family=binomial)
mod2 <- glm(cbind(ping,noping)~rech+tagh*loc,family=binomial)
mod3 <- glm(cbind(ping,noping)~rech*tagh*loc,family=binomial)

modA <- glm(cbind(ping,noping)~tagh+rech,family=binomial,data=pdA)
## The GLMs are not very succesful in finding an appropriate model (and it might not be appropriate at all to apply a GLM to these data since there is a large difference in the variation between the data sets [heteroskedasticity])


## --- Multivariate time series analysis (chapter 9 in HM's book)
nt <- dim(windat)[3]
X <- matrix(0,nt,3)
I <- 3
X[,1] <- windat[I,1,] 
X[,2] <- windat[I,2,] 
X[,3] <- windat[I,3,] 
Y <- matrix(0,nt,3)
J <- 1
Y[,1] <- windat[1,J,] 
Y[,2] <- windat[2,J,] 
Y[,3] <- windat[3,J,] 
##ccf(Y[,1],Y[,2])

pn2pw <- function(alpha,theta,sigma,Salpha,Stheta,Ssigma){
  c(alpha[which(Salpha==1)],theta[which(Stheta==1)],sigma[which(Ssigma==1)])
}
pw2pn <- function(pv,Salpha,Stheta,Ssigma){
  m <- dim(Stheta)[1]
  alpha <- rep(0,m)
  theta <- matrix(0,m,m)
  sigma <- matrix(0,m,m)
  inds <- which(Salpha==1)
  n1 <- length(inds)
  alpha[inds] <- pv[1:n1]
  inds <- which(Stheta==1)
  n2 <- length(inds)
  theta[inds] <- pv[(n1+1):(n1+n2)]
  inds <- which(Ssigma==1)
  n3 <- length(inds)
  sigma[inds] <- pv[(n1+n2+1):(n1+n2+n3)]
  sigma[lower.tri(sigma)] <- sigma[upper.tri(sigma)]
  list(alpha=alpha,theta=theta,sigma=sigma)
}

likfun <- function(pv,X,Salpha,Stheta,Ssigma){
  nt <- dim(X)[1]
  m <- dim(X)[2]
  pn <- pw2pn(pv,Salpha,Stheta,Ssigma)
  a <- matrix(rep(pn$alpha,nt-1),m,nt-1)
  Xhat <- pn$theta %*% (t(X[-nt,]) - a) ## Predictions
  Xtilde <- (t(X[-1,])-pn$alpha) - Xhat ## One-step prediction errors
  lik <- 0
  invsigma <- solve(pn$sigma)
  for(t in 1:(nt-1)){
    lik <- lik + Xtilde[,t] %*% invsigma %*% Xtilde[,t]
  }
  0.5*lik
}

## Estimate AR model for X data series
mod <- ar(X)
theta <- matrix(mod$ar,3,3)
alpha <- as.vector(mod$x.mean)
alpha <- c(0.8,0.8,0.8)
sigma <- matrix(mod$var.pred,3,3)

pn <- pw2pn(pv,Salpha,Stheta,Ssigma)

## Guesses
Gtheta <- diag(rep(1,3))
Gsigma <- 0.01*diag(rep(1,3))
Galpha <- rep(0.8,3)

## Parameter structures
Stheta <- matrix(0,3,3)
Stheta[1,] <- c(1,0,0)
Stheta[2,] <- c(0,1,0)
Stheta[3,] <- c(0,0,1)
Ssigma <- matrix(0,3,3)
Ssigma[1,] <- c(1,1,1)
Ssigma[2,] <- c(0,1,1)
Ssigma[3,] <- c(0,0,1)
Salpha <- c(1,1,1)

## pv is used as guess in optimization
pv <- pn2pw(alpha,theta,sigma,Salpha,Stheta,Ssigma)
##likfun(pv,X,Salpha,Stheta,Ssigma)

## Optimize parameters (be careful with parameters that imply non-stationarity, need to look more into this before estimating all parameters. Alternatively investigate models with similar sigma, and/or alpha and/or theta)
op <- optim(pv,likfun,gr=NULL,X,Salpha,Stheta,Ssigma,method='BFGS')


## A process is stationary if all the roots of the characteristic equation is within the unit circle p. 121 in HM book. The roots can be found using the below function. The characteristic equation is given by the difference equation as described in Appendix 1 in HM book.
polyroot(c(-0.9,2.8,-2.9,1))


likfun(op$par,X,Salpha,Stheta,Ssigma)

##plot(Xp[,1],typ='l')
##lines(X[-1,1],col=3)


## Below test the ar function and shows how the residuals are calculated
do.extra <- FALSE
if(do.extra){
  mod <- ar(Y)
  mod$order
  theta <- matrix(mod$ar,3,3)
  mu <- mod$x.mean
  XX <- t(X)
  X[2,]-mod$x.mean - as.vector(theta %*% X[1,] - theta %*% mod$x.mean )
  t(X[2,]-mod$x.mean) - theta %*% (X[1,]-mod$x.mean) 
  mod$resid[3,]
  acf(mod$resid[-1,])
}

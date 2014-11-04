## Script for calculating number of received pings
## Question: How accurate is the information from Vemco that 39xxx tags send pings uniformly random with a delay of 60-180 seconds? (it seems from exp 2 that either this information is incorrect or alternatively many pings are lost because of background noise?)
## 20.02.2013
rm(list=ls())

dat <- read.csv('VR100_10320_D2013.02.25T07.33.16.csv')

nn <- dim(dat)[1]
i <- 1
time <- as.POSIXct(paste(as.character(dat$Date[i]),as.character(dat$Time[i])),tz='UTC')
for(i in 2:nn) time[i] <- as.POSIXct(paste(as.character(dat$Date[i]),as.character(dat$Time[i])),tz='UTC')
##timeraw <- as.POSIXlt(time,tz='HST')
time <- as.POSIXlt(time,tz='UTC')
data <- dat

pr <- 120 ## Average pinging interval in seconds (between 60 and 180 sec)

strt <- as.POSIXct("2013-02-23 00:00:01",tz="UTC")
end <- as.POSIXct("2013-02-25 17:30:01",tz="UTC")
dt1 <- as.numeric(difftime(end,strt,unit='secs'))
expectedpings <- dt1/pr

inds <- time>strt & time<end
time1 <- time[inds]
data1 <- data[inds,c(7,17)] ## Tag IDs and their signal strenght in dB

ndp <- length(time1)

diffs <- as.numeric(difftime(time1[-ndp],time1[-1],units='sec'))

inds <- diffs > 63.5 & diffs < 183.5
diffs2 <- diffs[inds]

(res <- chisq.test(hst$counts))

breaks <- seq(63.5,183.5,length=13)
nbin <- length(breaks)-1

pdf('pingrate.pdf')
hist(diffs2,breaks=breaks,main=paste('Ping rate, tag 39193 (experiment 3). p-value:',round(res$p.value,3)),xlab='Time, sec',ylab='Number of received pings')
dev.off()

## Test if distribution of pings is truly uniform as claimed by Vemco













### Junk below this ------------------------------------


##m1 <- lm(diffs2~1) ## Data are not normal, so this is not appropriate
##(est <- cbind(Estimate = coef(m1), confint(m1)))
##plot(time1,rep(1,ndp),pch=3)


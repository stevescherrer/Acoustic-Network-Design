## Script for calculating number of received pings
## Question: Does the number of received pings per transmitter decrease when the number of transmitters increase?
## 16.02.2013
rm(list=ls())

dat <- read.csv('VR100_10320_D2013.02.20T07.51.36.csv')

nn <- dim(dat)[1]
i <- 1
time <- as.POSIXct(paste(as.character(dat$Date[i]),as.character(dat$Time[i])),tz='UTC')
for(i in 2:nn) time[i] <- as.POSIXct(paste(as.character(dat$Date[i]),as.character(dat$Time[i])),tz='UTC')
##timeraw <- as.POSIXlt(time,tz='HST')
time <- as.POSIXlt(time,tz='UTC')
data <- dat

pr <- 120 ## Average pinging interval in seconds (between 60 and 180 sec)

## --- Collision test, 4 tags
notag1 <- 4 ## Four tags in this test
strt <- as.POSIXct("2013-02-15 23:17:28",tz="UTC")
end <- as.POSIXct("2013-02-16 00:17:34",tz="UTC")
dt1 <- as.numeric(difftime(end,strt,unit='secs'))
expectedpings <- dt1/pr*notag1
inds <- time>strt & time<end
time1 <- time[inds]
data1 <- data[inds,c(7,17)] ## Tag IDs and their signal strenght in dB
receivedpings <- length(time1)
ids <- unique(data1[,1])
noids <- length(ids)
tagpings <- rep(0,noids)
meandiff1 <- rep(0,noids)
for(i in 1:noids){
  tagpings[i] <- length(which(data1[,1]==ids[i]))
  meandiff1[i] <- mean(diff(time1[data1[,1]==ids[i]]))
}
rat1 <- receivedpings/expectedpings
nocoll1 <- expectedpings-receivedpings

## --- Collision test, 8 tags
notag2 <- 8
strt <- as.POSIXct("2013-02-16 00:17:38",tz="UTC")
end <- as.POSIXct("2013-02-16 01:18:27",tz="UTC")
dt2 <- as.numeric(difftime(end,strt,unit='secs'))
expectedpings <- dt2/pr*notag2
inds <- time>strt & time<end
time2 <- time[inds]
data2 <- data[inds,c(7,17)] ## Tag IDs and their signal strenght in dB
receivedpings <- length(time2)
ids <- unique(data2[,1])
noids <- length(ids)
tagpings <- rep(0,noids)
meandiff2 <- rep(0,noids)
for(i in 1:noids){
  tagpings[i] <- length(which(data2[,1]==ids[i]))
  meandiff2[i] <- mean(diff(time2[data2[,1]==ids[i]]))
}
rat2 <- receivedpings/expectedpings
nocoll2 <- expectedpings-receivedpings

## --- Collision test, 16 tags
notag3 <- 16
strt <- as.POSIXct("2013-02-16 01:18:31",tz="UTC")
end <- as.POSIXct("2013-02-16 02:17:01",tz="UTC")
expectedpings <- dt2/pr*notag3
inds <- time>strt & time<end
time3 <- time[inds]
data3 <- data[inds,c(7,17)] ## Tag IDs and their signal strenght in dB
receivedpings <- length(time3)
ids <- unique(data3[,1])
noids <- length(ids)
tagpings <- rep(0,noids)
meandiff3 <- rep(0,noids)
for(i in 1:noids){
  tagpings[i] <- length(which(data3[,1]==ids[i]))
  meandiff3[i] <- mean(diff(time3[data3[,1]==ids[i]]))  
}
rat3 <- receivedpings/expectedpings
nocoll3 <- expectedpings-receivedpings

x <- c(4,8,16)
y <- c(rat1,rat2,rat3)
pdf('collresult.pdf')
plot(x,y,ylim=c(0,1),xlab='Number of tags in the water',ylab='Proportion received',main='Collision test (Experiment 2)')
lines(x,y,col=2)
dev.off()

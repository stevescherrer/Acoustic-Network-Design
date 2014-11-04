## Script for reading the VR100 data
## 11.02.2013
rm(list=ls())

dat <- read.csv('VR100_10320_D2013.02.11T08.42.22.csv')

nn <- dim(dat)[1]
i <- 1
time <- as.POSIXct(paste(as.character(dat$Date[i]),as.character(dat$Time[i])),tz='UTC')
for(i in 2:nn) time[i] <- as.POSIXct(paste(as.character(dat$Date[i]),as.character(dat$Time[i])),tz='UTC')
timeraw <- as.POSIXlt(time,tz='HST')

strt <- as.POSIXct("2013-02-08 23:00:01",tz='UTC') ## Starting time for observations
inds <- timeraw > strt
time <- timeraw[inds]
data <- dat[inds,]
##inds <- which(data$Chan==1)
##time <- time[inds]
##data <- data[inds,]

## --- Tag signal strenght testing
strt <- as.POSIXct("2013-02-09 00:15:01",tz="UTC")
end <- as.POSIXct("2013-02-09 00:46:01",tz="UTC")
inds <- time>strt & time<end
time1 <- time[inds]
data1 <- data[inds,c(7,17)] ## Tag IDs and their signal strenght in dB


## --- Signal strenght in different locations (1, 2, 3, 4), no tag casing used
## Location 1, Ho'okele
strt <- as.POSIXct("2013-02-09 00:51:01",tz="UTC")
end <- as.POSIXct("2013-02-09 00:57:01",tz="UTC")
inds <- time>strt & time<end
time2 <- time[inds]
data2 <- data[inds,c(7,17)] ## Tag IDs and their signal strenght in dB

## Location 2, floating dock, far
strt <- as.POSIXct("2013-02-09 00:57:01",tz="UTC")
end <- as.POSIXct("2013-02-09 01:02:01",tz="UTC")
inds <- time>strt & time<end
time3 <- time[inds]
data3 <- data[inds,c(7,17)] ## Tag IDs and their signal strenght in dB

## Location 3, floating dock, near
strt <- as.POSIXct("2013-02-09 01:02:01",tz="UTC")
end <- as.POSIXct("2013-02-09 01:12:01",tz="UTC")
inds <- time>strt & time<end
time4 <- time[inds]
data4 <- data[inds,c(7,17)] ## Tag IDs and their signal strenght in dB

## Location 4, red raft (Ho'okele usual parking space)
strt <- as.POSIXct("2013-02-09 01:12:01",tz="UTC")
end <- as.POSIXct("2013-02-09 01:18:01",tz="UTC")
inds <- time>strt & time<end
time5 <- time[inds]
data5 <- data[inds,c(7,17)] ## Tag IDs and their signal strenght in dB


## --- Signal strenght in different cases at location 1 (Ho'okele)
## Small case, no holes
strt <- as.POSIXct("2013-02-09 01:18:01",tz="UTC")
end <- as.POSIXct("2013-02-09 01:30:01",tz="UTC")
inds <- time>strt & time<end
time6 <- time[inds]
data6 <- data[inds,c(7,17)] ## Tag IDs and their signal strenght in dB

## Small case, with holes
strt <- as.POSIXct("2013-02-09 01:30:01",tz="UTC")
end <- as.POSIXct("2013-02-09 01:37:01",tz="UTC")
inds <- time>strt & time<end
time7 <- time[inds]
data7 <- data[inds,c(7,17)] ## Tag IDs and their signal strenght in dB

## Large case
strt <- as.POSIXct("2013-02-09 01:37:01",tz="UTC")
end <- as.POSIXct("2013-02-09 01:43:01",tz="UTC")
inds <- time>strt & time<end
time8 <- time[inds]
data8 <- data[inds,c(7,17)] ## Tag IDs and their signal strenght in dB


## --- Collision test, location 1 (Ho'okele), no casings used
## Two tags: 39193, 39194
strt <- as.POSIXct("2013-02-09 01:43:01",tz="UTC")
end <- as.POSIXct("2013-02-09 01:50:01",tz="UTC")
inds <- time>strt & time<end
time9 <- time[inds]
data9 <- data[inds,c(7,17)] ## Tag IDs and their signal strenght in dB

## Five tags: 39193-39197
strt <- as.POSIXct("2013-02-09 01:50:01",tz="UTC")
end <- as.POSIXct("2013-02-09 02:06:01",tz="UTC")
inds <- time>strt & time<end
time10 <- time[inds]
data10 <- data[inds,c(7,17)] ## Tag IDs and their signal strenght in dB

## --- Transmission times (time is takes to transmit an ID code) (recorded by Chase using a hand held stop watch)
ttimes <- c(2.84,2.7,2.97,2.67,2.94,2.89,3.04,3.17,3.14,2.64,2.91,2.9) + 0.15 ## Add 0.15 s for human delay

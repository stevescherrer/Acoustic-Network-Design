## Script for simulating tag pings
## 08.02.2013
rm(list=ls())

prmin <- 60
prmax <- 180
T <- 2*60*60 ## Time left in the water
n <- 200

NN <- rep(0,n)
for(i in 1:n){
  t <- runif(1,prmin,prmax)
  c <- 1
  while(max(t)<T){
    c <- c+1
    t[c] <- t[c-1] + runif(1,prmin,prmax)
  }
  NN[i] <- length(t)
}

mu <- mean(NN)
sd <- sqrt(var(NN))
cv <- sd/mu
print(cv)

plot(density(NN))

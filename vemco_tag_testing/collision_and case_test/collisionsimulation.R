## Script for simulating the number of collisions
## 08.02.2013
rm(list=ls())

ntag <- 8 ## Number of tags (must be at least 2 to get collisions, duh)
tt <- 3 ## Transmission time in seconds
prmin <- 60
prmax <- 180
hours <- 1
T <- hours*60*60 ## Time left in the water

tags <- 1:ntag
ntrial <- 1000
ncols <- rep(0,ntrial)
for(r in 1:ntrial){
  NN <- list()
  for(i in 1:ntag){
    t <- runif(1,prmin,prmax)
    c <- 1
    while(max(t)<T){
      c <- c+1
      t[c] <- t[c-1] + runif(1,prmin,prmax)
    }
    NN[[i]] <- t
  }
  ncol <- 0
  for(i in tags){ ## Loop over tags
    nping <- length(NN[[i]])
    for(k in tags[-i]){ ## Loop over other tags
      if(k>i){ ## Avoid checking 1 versus 2 when 2 versus 1 has already been checked
        for(j in 1:nping){ ## Loop over pings
          ts <- NN[[i]][j]
          te <- NN[[i]][j] + tt
          ncol <- ncol + length(which( NN[[k]]>ts & NN[[k]]<te ))
          ncol <- ncol + length(which( (NN[[k]]+tt)>ts & (NN[[k]]+tt)<te ))
          ##if(any(NN[[k]]>ts & NN[[k]]<te)) print(paste(i,k,j))
          ##if(any((NN[[k]]+tt)>ts & (NN[[k]]+tt)<te)) print(paste(i,k,j))
        }
      }
    }
  }
  ncols[r] <- ncol
}

do.save <- FALSE
if(do.save) pdf('colldens.pdf')
hist(ncols,xlim=c(0,80),freq=0,main=paste('No tag:',ntag,' Hours:',hours,' Ping rate:',mean(c(prmin,prmax)),'s  Trans time:',tt),xlab='Number of collisions')
if(do.save) dev.off()

plot.seg <- FALSE
if(plot.seg){
  I <- 1
  II <- rep(I,length(NN[[I]]))

  plot(NN[[1]][1],II[1],ylim=c(0.5,ntag+0.5),xlim=range(NN[[1]]))
  segments(NN[[1]],II,NN[[1]]+tt,II)

  for(i in 2:ntag){
    ii <- rep(i,length(NN[[i]]))
    points(NN[[i]],ii,col=i)
  }
}

rm(list=ls())
library(zoo)

err1 <- read.csv("myErr-1.csv")
err2 <- read.csv("myErr-2.csv")
err3 <- read.csv("myErr-3.csv")
# err1 <- err1[,-1]
# err2 <- err2[,-1]
# err3 <- err3[,-1]

myPlot <- matrix(NA,4,1000)
myPlot[1,] <- 0:999 / 100
tempPlot <- matrix(NA,2,1000)
tempPlot[1,] <- 0:999 / 100

tempPlot <- cbind(tempPlot,err1)
interpi <- as.data.frame(t(tempPlot))
interpi <- zoo(interpi)
index(interpi) <- interpi[,1]
interpi <- na.approx(interpi, na.rm=F)
tempPlot <- t(as.matrix(interpi))
tempPlot <- tempPlot[,1:1000]
colnames(tempPlot) <- rownames(tempPlot) <- NULL
myPlot[2,] <- tempPlot[2,]

tempPlot <- matrix(NA,2,1000)
tempPlot[1,] <- 0:999 / 100

tempPlot <- cbind(tempPlot,err2)
interpi <- as.data.frame(t(tempPlot))
interpi <- zoo(interpi)
index(interpi) <- interpi[,1]
interpi <- na.approx(interpi, na.rm=F)
tempPlot <- t(as.matrix(interpi))
tempPlot <- tempPlot[,1:1000]
colnames(tempPlot) <- rownames(tempPlot) <- NULL
myPlot[3,] <- tempPlot[2,]

tempPlot <- matrix(NA,2,1000)
tempPlot[1,] <- 0:999 / 100

tempPlot <- cbind(tempPlot,err3)
interpi <- as.data.frame(t(tempPlot))
interpi <- zoo(interpi)
index(interpi) <- interpi[,1]
interpi <- na.approx(interpi, na.rm=F)
tempPlot <- t(as.matrix(interpi))
tempPlot <- tempPlot[,1:1000]
colnames(tempPlot) <- rownames(tempPlot) <- NULL
myPlot[4,] <- tempPlot[2,]

par(mfrow=c(1,1))
plot(myPlot[1,], rep(0,1000), type="l", ylim=c(0, 4), 
     xlab="t value", ylab="predictive error")
title("average 0.632+ bootstrap estimate", 
      cex.main = 2,   font.main= 1)
lines(myPlot[1,], myPlot[2,], lwd=1)
lines(myPlot[1,], myPlot[3,], lwd=1)
lines(myPlot[1,], myPlot[4,], lwd=1)
lines(myPlot[1,], colMeans(myPlot[2:4,],na.rm=T), lwd=4, col="RED")
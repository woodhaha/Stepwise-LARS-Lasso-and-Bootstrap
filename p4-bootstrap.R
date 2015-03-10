# Clean memory, read data, include sources
# install.packages("lars")
# install.packages("zoo")
rm(list=ls())
library(lars)
library(zoo)
source("p4-util-myLasso.R")
source("p4-util-processData.R")

myData <- read.csv("myData-3.csv")
# Remove row names
myData <- myData[,2:dim(myData)[2]]
Y <- processData(myData)$dataY
X <- processData(myData)$dataX
dn <- dim(X)[1]
dm <- dim(X)[2]

# -------------------- SET UP BOOTSTRAP -------------------- #

# Set up the bootStrap
numBoot <- 25
bootIndex <- list()
bootSet <- matrix(F,dn,numBoot)
myBeta <- list()
k <- rep(0,numBoot)
# Calculate all the bootstraps and store the result
b<-1
for(b in 1:numBoot) {
  bootIndex[[b]] <- sample.int(29,  size=29, replace=T)
  bootData <- myData[bootIndex[[b]],]
  myY <- processData(bootData)$dataY
  myX <- processData(bootData)$dataX
  
  bootSet[,b][bootIndex[[b]]] <- T
  
  myBeta[[b]] <- myLasso(myX, myY, numUnique=(sum(bootSet[,b])-1))
#   myBeta[[b]] <- t(lars(myX, myY, type="lasso")$beta)
  k[b] <- dim(myBeta[[b]])[2]
}

# -------------------- COMPUTE ERROR -------------------- #

mySum <- list()
# Calculate the sum of Q's for each i
for(i in 1:dn) {
  # Trace all the errors for a i
  erri <- list()
  ti <- list()
  Ccomp <- seq(1:numBoot)[!bootSet[i,]]
  # Calculate all the Q's of different bootstrap for a i
  for(b in Ccomp) {
    errib <- rep(0,k[b])
    tib <- rep(0,k[b])
    for(tt in 1:k[b]) {
      yiHat <- X[i,] %*% myBeta[[b]][,tt]
      errib[tt] <- abs(Y[i] - yiHat)
      tib[tt] <- sum(abs(myBeta[[b]][,tt])) 
    }
    erri[[which(Ccomp==b)]] <- errib
    ti[[which(Ccomp==b)]] <- tib
  }
  # Arrange the Q's
  nt <- length(do.call("c",erri))
  sumi <- matrix(NA,length(Ccomp)+1,nt)
  sumi[1,] <- do.call("c",ti)
  start <- 1
  for(b in Ccomp) {
    sumi[which(Ccomp==b)+1, start:(start+k[b]-1)] <- 
      erri[[which(Ccomp==b)]]
    start <- start+k[b]
  }
  # Sort the Q's according to t
  lineOrder <- order(sumi[1,])
  sumi <- sumi[,lineOrder]
  # Fill the matrix
  interpi <- as.data.frame(t(sumi))
  interpi <- zoo(interpi)
  index(interpi) <- interpi[,1]
  interpi <- na.approx(interpi, na.rm=F)
  sumi <- t(as.matrix(interpi))
  colnames(sumi) <- rownames(sumi) <- NULL
  # Arrange the sum accordingly
  mySumi <- matrix(NA,dn+1,dim(sumi)[2])
  mySumi[1,] <- sumi[1,]
  mySumi[i+1,] <- colMeans(sumi[-1,], na.rm=T)
  mySum[[i]] <- mySumi
}
# Final arranges
mySum <- do.call("cbind",mySum)
lineOrder <- order(mySum[1,])
mySum <- mySum[,lineOrder]
interpi <- as.data.frame(t(mySum))
interpi <- zoo(interpi)
index(interpi) <- interpi[,1]
interpi <- na.approx(interpi, na.rm=F)
mySum <- t(as.matrix(interpi))
colnames(mySum) <- rownames(mySum) <- NULL

# -------------------- FULL MODEL ERROR -------------------- #

larModel <- lars(X, Y, type="lasso")
larBeta <- t(larModel$beta)
k <- dim(larBeta)[2]
t <- rep(0,k)
gamma <- rep(0,k)

i <- j <- tt <- 2

for(tt in 1:k) {
  for(i in 1:dn) {for(j in 1:dn) {
    gamma[tt] <- gamma[tt] + abs(Y[i] - (X[j,] %*% larBeta[,tt]))
  }}
  t[tt] <- sum(abs(larBeta[,tt]))
}
gamma <- gamma / (dn^2)
myGamma <- matrix(NA,2,dim(mySum)[2]+k)
myGamma[1,1:dim(mySum)[2]] <- mySum[1,]
myGamma[1,-(1:dim(mySum)[2])] <- t
myGamma[2,-(1:dim(mySum)[2])] <- gamma

interpi <- as.data.frame(t(myGamma))
interpi <- zoo(interpi)
index(interpi) <- interpi[,1]
interpi <- na.approx(interpi, na.rm=F)
myGamma <- t(as.matrix(interpi))
myGamma <- myGamma[2,1:dim(mySum)[2]]

# -------------------- FULL MODEL ERROR -------------------- #

larModel <- lars(X, Y, type="lasso")
larBeta <- t(larModel$beta)
k <- dim(larBeta)[2]
t <- rep(0,k)
err <- rep(0,k)
for(i in 1:k) {
  yHat <- X %*% larBeta[,i]
  t[i] <- sum(abs(larBeta[,i]))
  err[i] <- mean(abs(Y - yHat))
}
errBar <- matrix(NA,2,dim(mySum)[2]+k)
errBar[1,1:dim(mySum)[2]] <- mySum[1,]
errBar[1,-(1:dim(mySum)[2])] <- t
errBar[2,-(1:dim(mySum)[2])] <- err

interpi <- as.data.frame(t(errBar))
interpi <- zoo(interpi)
index(interpi) <- interpi[,1]
interpi <- na.approx(interpi, na.rm=F)
errBar <- t(as.matrix(interpi))
errBar <- errBar[2,1:dim(mySum)[2]]

# -------------------- 632+ ERROR TOWARDS T -------------------- #

errHat <- colMeans(mySum[-1,], na.rm=T)
errHat <- c(errHat)
errBar <- c(errBar)
myGamma <- c(myGamma)

rHat <- (errHat-errBar) / (myGamma-errBar)
omegaHat <- 0.632 / (1 - (0.368*rHat))
errFinal <- ((1-omegaHat) * errBar) + (errHat * omegaHat)

write.csv(rbind(t,errFinal),file="myErr-3.csv")

# -------------------- PLOT BOOTSTRAP ERROR -------------------- #

par(mfrow=c(2,2))
t <- c(mySum[1,])
k <- length(t)

# ---------------------------------------- #

plot(t, rep(0,k), type="l", ylim=c(0, 4), 
     xlab="t value", ylab="absolute error")
title("1 bootstrap estimate", 
      cex.main = 2,   font.main= 1)
for(i in 1:dn+1) { lines(t, mySum[i,], lwd=1) }
lines(t, colMeans(mySum[-1,], na.rm=T), lwd=4, col="RED")

# ---------------------------------------- #

for(i in 1:dn+1) { mySum[i,] <- 0.632*c(mySum[i,]) + 0.368*c(errBar) }
plot(t, rep(0,k), type="l", ylim=c(0, 3), 
     xlab="t value", ylab="absolute error")
title("0.632 bootstrap estimate", 
      cex.main = 2,   font.main= 1)
for(i in 1:dn+1) { lines(t, mySum[i,], lwd=1) }
lines(t, colMeans(mySum[-1,], na.rm=T), lwd=4, col="RED")

# ---------------------------------------- #

plot(t, rep(0,k), type="l", ylim=c(0, 4), 
     xlab="t value", ylab="overfitting rate")
title("relative rate of overfitting", 
      cex.main = 2,   font.main= 1)
lines(t, rHat, lwd=4, col="RED")

# ---------------------------------------- #

plot(t, rep(0,k), type="l", ylim=c(0, 4), 
     xlab="t value", ylab="predictive error")
title("0.632+ bootstrap estimate", 
      cex.main = 2,   font.main= 1)
lines(t, errFinal, lwd=4, col="RED")
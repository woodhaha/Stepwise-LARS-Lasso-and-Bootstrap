rm(list=ls())
myData <- read.csv("myData-3.csv")
myData <- myData[,2:dim(myData)[2]]

myX <- as.matrix(myData[,2:dim(myData)[2]])
myY <- myData[,1]
dn <- dim(myX)[1]
dm <- dim(myX)[2]
myX <- scale(myX)
myY <- scale(myY)

logicRemove <- F
stepLength <- 0.001
maxSteps <- 150
XA <- NULL
activeSet <- rep(F, dm)
mu <- rep(0,dn)
error <- myY-mu
B <- matrix(0,dm,maxSteps)
beta <- rep(0,dm)
BETA <- matrix(0,dm,maxSteps)
cHat <- t(myX) %*% error
logicNonSingularity <- T

SST <- t(myY) %*% myY
SSE <- rep(0,maxSteps)
Rsquare <- rep(0,maxSteps)
varSeq <- list()

myThre <- 10^-15
k=1
while(k<=maxSteps && logicNonSingularity) {
  error <- myY-mu
  if(!logicRemove) {
    errCor <- cor(myX, error)
    errCor[activeSet] = 0
    index <- which.max(abs(errCor))
    activeSet[index] <- T
    varSeq[k] <- paste("variable no.",index,"added:",
                       names(as.data.frame(myX))[index])
  } else {
    index <- which(rmVar==names(as.data.frame(myX)))
    activeSet[index] <- F
    varSeq[k] <- paste("variable no.",index,"droped:",
                       names(as.data.frame(myX))[index])
  }
  sign <- sign(cor(myX, error))
  signId <- diag(dm)
  for(i in 1:dm) {
    signId[i,i] <-sign[i]
  }  
  XA <- (myX %*% signId)[,activeSet]
  XA <- as.matrix(XA)
  one <- rep(1,sum(activeSet))
  cHat <- t(myX) %*% (error)
  cCap <-max(abs(cHat))
  a <- t(myX) %*% XA %*% solve(t(XA) %*% XA) %*% one
  
  gamma <- c(((cCap-cHat)/(1-a))[!activeSet], 
             ((cCap+cHat)/(1+a))[!activeSet])
  gamma <- min(gamma[gamma>myThre])
  # lasso
  bTilde <- solve(t(XA) %*% XA) %*% one
  cStar <- -(BETA[activeSet,k-1] / (sign[activeSet]*bTilde))
  cc=Inf
  if(length(cStar[cStar>myThre])==0) {
    cc=Inf
    indexStar <- 0
    rmVar <- NULL
  } else {
    cc <- min(cStar[cStar>myThre])
    indexStar <- which(cStar==cc)
    rmVar <- names(as.data.frame(myX))[activeSet][indexStar]
  }
  logicRemove <- cc < gamma
  if(logicRemove) { gamma <- cc }
  
  u <- gamma * (XA %*% solve(t(XA) %*% XA) %*% one)
  mu <- mu + u
  b <- gamma * (solve(t(XA) %*% XA) %*% one)
  beta[activeSet] <- beta[activeSet] + b*sign[activeSet]
  BETA[,k] <- beta
  error <- myY-mu
  cHat <- t(myX) %*% error
  
  SSE[k] <- t(error) %*% error
  Rsquare[k] <- 1 - (SSE[k]/SST)
  
  k <- k+1
  logicNonSingularity <- max(0,dim(XA)[2])<max(28,dim(XA)[1]-1)
  if(logicRemove) {logicNonSingularity <- T}
}
varSeq <- sapply(varSeq,function(v){return(v)})
k <- k-1
BETA <- BETA[,1:k]
t <- colSums(abs(BETA))

plot(t, rep(0,k), type="l", ylim=c(-1.2, 1.2), 
     xlab="t value", ylab="beta value")
title("Lasso Path No.1", 
      cex.main = 2,   font.main= 1)
for(i in 1:dm) { lines(t, BETA[i,], lwd=2, col=sample(colours(), 1)) }
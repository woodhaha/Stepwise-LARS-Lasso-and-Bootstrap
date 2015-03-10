rm(list=ls())
myData <- read.csv("myData-3.csv")
# Remove row names
myData <- myData[,2:dim(myData)[2]]

# -------------------- ITERATIONS -------------------- #

myX <- as.matrix(myData[,2:dim(myData)[2]])
myY <- myData[,1]
dn <- dim(myX)[1]
dm <- dim(myX)[2]
# Why it is different from colMeans(myX)?
myX <- scale(myX, drop(rep(1,dn) %*% myX)/dn, FALSE)
myX <- scale(myX, FALSE, sqrt(drop(rep( 1,dn) %*% (myX^2))))
myY <- drop(myY - mean(myY))

# -------------------- ITERATIONS -------------------- #

logicRemove <- F
stepLength <- 0.001
maxSteps <- 28
myThre <- 10^-5
k=1
XA <- NULL
activeSet <- rep(F, dm)
mu <- rep(0,dn)
error <- myY-mu
B <- matrix(0,dm,maxSteps)
beta <- rep(0,dm)
BETA <- matrix(0,dm,maxSteps)
cHat <- t(myX) %*% error
logicNonSingularity <- T
# -------------------- #
SST <- t(myY) %*% myY
SSE <- rep(0,maxSteps)
Rsquare <- rep(0,maxSteps)
# -------------------- #
varSeq <- list()
# -------------------- #

while(k<=maxSteps && logicNonSingularity) {
  error <- myY-mu
  errCor <- cor(myX, error)
  errCor[activeSet] = 0
  index <- which.max(abs(errCor))
  activeSet[index] <- T
  varSeq[k] <- paste("variable no.",index,"added:",
                     names(as.data.frame(myX))[index])

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
  
  u <- gamma * (XA %*% solve(t(XA) %*% XA) %*% one)
  mu <- mu + u
  b <- gamma * (solve(t(XA) %*% XA) %*% one)
  beta[activeSet] <- beta[activeSet] + b*sign[activeSet]
  BETA[,k] <- beta
  error <- myY-mu
  cHat <- t(myX) %*% error
  # -------------------- #
  SSE[k] <- t(error) %*% error
  Rsquare[k] <- 1 - (SSE[k]/SST)
  # -------------------- #
  k <- k+1
  logicNonSingularity <- max(0,dim(XA)[2])<max(28,dim(XA)[1]-1)
  if(logicRemove) {logicNonSingularity <- T}
}
t <- colSums(abs(BETA))
beta <- solve(t(XA) %*% XA) %*% t(XA) %*% myY
varSeq <- sapply(varSeq,function(var){return(var)})
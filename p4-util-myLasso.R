myLasso <- function(myX, myY, numUnique=29) {
  method <- "lasso"
  logicRemove <- F
  stepLength <- 0.001
  maxSteps <- 150
  myThre <- 10^-15
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
  
  SST <- t(myY) %*% myY
  SSE <- rep(0,maxSteps)
  Rsquare <- rep(0,maxSteps)
  varSeq <- list()
  
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
    if(method == "lasso") {
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
    }
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
    logicNonSingularity <- max(0,dim(XA)[2])<min(numUnique,dim(XA)[1]-1)
    if(logicRemove) {logicNonSingularity <- T}
  }
  k <- k-1
  BETA <- BETA[,1:k]
  varSeq <- sapply(varSeq,function(var){return(var)})
  
  return(BETA)
}
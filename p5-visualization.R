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
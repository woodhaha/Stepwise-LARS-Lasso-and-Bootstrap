processData <- function(myData) {
  myX <- as.matrix(myData[,2:dim(myData)[2]])
  myY <- myData[,1] 
  myX <- scale(myX)
  myY <- scale(myY)
  
  result <- list(myY, myX)
  names(result) <- c("dataY", "dataX")
  
  return(result)
}
removeUp <- function(tempData) {
  # Remove completely missing columns and rows with completely missing predictors
  logicMissingCol <- as.numeric(colSums(is.na(tempData))) != dim(tempData)[1]
  logicMissingRow <- as.numeric(rowSums(is.na(tempData))) != dim(tempData)[2]-1
  tempData<- tempData[logicMissingRow,logicMissingCol]
  
  # Remove columns with constant predictors
  colVar <- apply(tempData, 2, function(xx){
    xx <- xx[!is.na(xx)]
    return(var(xx))
  })
  logicConstCol <- round(as.numeric(colVar), 3) != 0
  tempData <- tempData[,logicConstCol]
  
  # Remove (second or later occurrance of) columns which are collinear (correlation 1) with another.
  numNa <- as.numeric(colSums(is.na(tempData)))
  tempData_1 <- tempData[,numNa==0]
  tempData_1 <- tempData_1[!is.na(rowSums(tempData_1)),]
  tempData_2 <- tempData[,numNa==1]
  tempData_2 <- tempData_2[!is.na(rowSums(tempData_2)),]
  tempData_3 <- tempData[,numNa==3]
  tempData_3 <- tempData_3[!is.na(rowSums(tempData_3)),]
  
  myCor <- matrix(0, dim(tempData)[2], dim(tempData)[2])
  p1 <- dim(tempData_1)[2]
  p2 <- dim(tempData_2)[2]
  p3 <- dim(tempData_3)[2]
  myCor[1:p1,1:p1] <- as.numeric(round(cor(tempData_1), 3))
  myCor[p1+1:p2,p1+1:p2] <- as.numeric(round(cor(tempData_2), 3))
  myCor[p1+p2+1:p3,p1+p2+1:p3] <- as.numeric(round(cor(tempData_3), 3))
  
  logicCollMatrix <- myCor == 1
  logicColl <- rep(T,dim(tempData)[2])
  
  for(i in dim(logicCollMatrix)[1]:2) {
    for(j in (i-1):1) {
      if(logicCollMatrix[i,j]) {
        logicColl[i] <- F
      }
    }
  }
  
  tempData <- tempData[,logicColl]
  return(tempData)
}
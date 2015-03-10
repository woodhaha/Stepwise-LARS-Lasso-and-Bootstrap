arrange <- function(xData,yData){
  # Rearrange data
  myXtemp <- data.matrix(xData[,6:dim(xData)[2]])
  patientIDsTemp <- sort(unique(c(xData$PatID,yData$Pat.No)))
  timePointsTemp <- sort(unique(as.character(xData$TimePoint)))
  
  PFS=matrix(NA,length(patientIDsTemp),1)
  predictors=matrix(NA,length(patientIDsTemp),(dim(myXtemp)[2])*length(timePointsTemp))
  counter=1
  
  for(pt in patientIDsTemp){
    index=which(yData$Pat.No==pt)
    if(length(index)==1){
      PFS[counter]=log(yData$PFS[index])
    }
    for(jj in 1:length(timePointsTemp)){
      index=which((myXData$TimePoint==timePointsTemp[jj])&
                    (myXData$PatID==pt))
      if(length(index)==1){
        predictors[counter,((jj-1)*dim(myXtemp)[2]+1):(jj*dim(myXtemp)[2])]=myXtemp[index,]
      }
    }
    counter=counter+1
  }
  
  predictorsStar=apply(predictors,2,function(xx){
    if(sum(is.na(xx))>=2){
      return((xx-mean(xx,na.rm=TRUE))/sd(xx,na.rm=TRUE))
    }else{
      return(xx)
    }
  })
  
  data=data.frame(PFS,predictorsStar)
  data$PFS=data$PFS-mean(data$PFS,na.rm=TRUE)
  
  # Label columns
  myNames=c("PFS",
            paste(attributes(myXData)$names[6:dim(myXData)[2]],"_baseline1",sep=""),
            paste(attributes(myXData)$names[6:dim(myXData)[2]],"_baseline2",sep=""),
            paste(attributes(myXData)$names[6:dim(myXData)[2]],"_followUp",sep=""))
  attributes(data)$names=myNames
  
  return(data)
}
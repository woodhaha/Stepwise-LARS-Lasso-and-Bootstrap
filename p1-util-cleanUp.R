cleanUp <- function(xData){
  # Restrict attention to particular variables
  myLogical <- sapply(attributes(xData)$names,function(varname){
    if(varname %in% c("PatID","TimePoint","PFS")){
      return(TRUE)
    }else{
      return(((length(grep(pattern="v1_",varname))>0))&
               (!((length(grep(pattern="Xa",varname))>0)|
                    (length(grep(pattern="DeltaX",varname))>0))))
    }
  })
  xData <- xData[,myLogical]
  return(xData)
}
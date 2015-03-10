rm(list=ls())

misStar <- read.csv("misStar-3.csv")
misStar <- misStar[,2:dim(misStar)[2]]
misNames <- names(misStar)

base1 <- grep(pattern="_baseline1",misNames)
base2 <- grep(pattern="_baseline2",misNames)
follow <- grep(pattern="_followUp",misNames)
misBase1 <- misStar[,base1]
misBase2 <- misStar[,base2]
misFollow <- misStar[,follow]

names(misBase1) <- lapply(names(misBase1),function(name){
  return (substring(name,1,nchar(name)-nchar("_baseline1"))) })
names(misBase2) <- lapply(names(misBase2),function(name){
  return (substring(name,1,nchar(name)-nchar("_baseline2"))) })
names(misFollow) <- lapply(names(misFollow),function(name){
  return (substring(name,1,nchar(name)-nchar("_followUp"))) })

base1 <- names(misBase1)
base2 <- names(misBase2)
follow <- names(misFollow)
misNames <- unique(c(base1,base2,follow), fromLast=F)

misNames %in% base1
misNames %in% base2
misNames %in% follow

misAve <- sapply(misNames,function(name){
  b1 <- b2 <- rep(NA,dim(misStar)[1])
  if(name %in% base1) b1 <- misBase1[which(base1==name)]
  if(name %in% base2) b2 <- misBase2[which(base2==name)]
  return(rowMeans(cbind(b1,b2),na.rm=T))
})
misDiff <- sapply(misNames,function(name){
  b1 <- b2 <- fl <- rep(NA,dim(misStar)[1])
  if(name %in% base1) b1 <- misBase1[which(base1==name)]
  if(name %in% base2) b2 <- misBase2[which(base2==name)]
  if(name %in% follow) fl <- misFollow[which(follow==name)]
  return(fl - rowMeans(cbind(b1,b2),na.rm=T))
})

misAve <- as.data.frame(misAve)
misDiff <- as.data.frame(misDiff)
names(misDiff) <- paste(names(misAve),"_difference",sep="")
names(misAve) <- paste(names(misAve),"_average",sep="")
myData <- cbind(misStar[1],misAve,misDiff)
write.csv(myData,"myData-3.csv")
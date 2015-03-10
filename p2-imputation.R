# Multiple imputation
install.packages("mi")
library(mi)
myInfo <- mi.info(myData)
myInfo <- update(myInfo,"type", 
                 as.list(rep("predictive-mean-matching",
                             dim(myData)[2])))
myMIs <- mi(myData,myInfo)
misStar <- mi.completed(myMIs)
# Save the data in files
dump(c("myInfo", "myMIs", "misStar"), file="miData.R")
write.csv(misStar[1], "misStar-1.csv")
write.csv(misStar[2], "misStar-2.csv")
write.csv(misStar[3], "misStar-3.csv")
# source("miData.R")
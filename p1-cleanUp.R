# Clean memory, include sources
rm(list=ls())
source("p1-util-cleanUp.R")
source("p1-util-arrange.R")
source("p1-util-removeUp.R")
# Read in texture and response data
myXData <- read.csv("Texture.csv",header=TRUE)
myYData <- read.csv("TextureResponse.csv",header=TRUE)
# Clean, arrange and manipulate the data
myXData <- cleanUp(myXData)
myData <- arrange(myXData, myYData)
myData <- removeUp(myData)
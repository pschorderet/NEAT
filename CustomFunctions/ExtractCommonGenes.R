# Set working directory
setwd("~/PepsRscripts/RScripts/MeiSheng/")
MS <- read.csv("Workbook5.csv")
head(MS)
# Create dataframe from .csv
df = as.data.frame(MS)
df1 = df[,1:2]
df2 = df[,3:4]
head(df1)
# Set paramters
OV = df1[df1[,1] %in% df2[,1],]

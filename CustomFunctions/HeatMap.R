#------------------------------------------------------------
# Load libraries
library(plotrix)
source("~/PepsRscripts/RScripts/PepsFunctions/RectPlotter.R")

#------------------------------------------------------------
# ColorBar
ColorBar <- function(lut, min, max=-min, nticks=5, ticks=seq(min, max, len=nticks), title='', key="bottom") {
  scale = (length(lut)-1)/(max-min)
  
  if(key=="bottom"){
    plot(c(min,max),c(0,1), axes=FALSE, xlab="",ylab="", main=title, col='white')
    axis(1, ticks, las=1, cex.axis=0.5, lwd=0.5)
    for (i in 1:(length(lut)-1)) {
      y = (i-1)/scale + min
      rect(y,0,y+(1/scale),0.05, col=lut[i], border=NA)
    }
  }
  if(key=="top"){
    plot(c(min,max),c(0,1), axes=FALSE, xlab="",ylab="", main=title, col='white')
    axis(3, ticks, las=1, cex.axis=0.5, lwd=0.5)
    for (i in 1:(length(lut)-1)) {
      y = (i-1)/scale + min
      rect(y,0.95,y+(1/scale),1, col=lut[i], border=NA)
    }
  }
  
  if(key=="left"){
    plot(c(0,1), c(min,max), axes=FALSE, xlab="",ylab="", main=title, col='white')
    axis(4, ticks, las=1, cex.axis=0.5, lwd=0.5)
    for (i in 1:(length(lut)-1)) {
      y = (i-1)/scale + min
      rect(0.95,y,1,y+1/scale, col=lut[i], border=NA)
    }  
  }
  
  if(key=="right"){
    plot(c(0,1), c(min,max), axes=FALSE, xlab="",ylab="", main=title, col='white')
    axis(2, ticks, las=1, cex.axis=0.5, lwd=0.5)
    for (i in 1:(length(lut)-1)) {
      y = (i-1)/scale + min
      rect(0,y,0.05,y+1/scale, col=lut[i], border=NA)
    }  
  }
}

#------------------------------------------------------------
# HeatMap
HeatMap <- function(mymatrix, datalim=c(0,0), signmatrix=NULL, colorgradient=c("blue", "red"), title="Main title", ycenter=seq(from=0,to=1,length=(length(mymatrix[1,]))), emptylines=0, key="bottom", cexsize=0.5){
  # mymatrix is a matrix as follows: start - end - value
  # colorlim will work if you do not supply datalim
  # if datalim is provided, colors will range from 
  
  # Set different parameters
  matrix3plus <- mymatrix[,-(1:2)]
  steps <- 1000
  rbPal <- colorRampPalette(colorgradient)
  gradcol <- rbPal(steps)
  # If datalim is not provided, set to min and max of mymatrix
  if(datalim[1]==0 & datalim[2]==0){
    datalim <- c(floor(100*min(matrix3plus)), ceiling(100*max(matrix3plus)))/100
  }
  mycolors <- cbind(mymatrix[,1:2], array(rbPal(steps)[cut(matrix3plus, seq(datalim[1],datalim[2],length=steps))], dim=dim(matrix3plus)))  
  # Plot a barplot as a control for the colors
  #barplot(mymatrix[,-(1:2)], beside=T, col=mycolors[,-(1:2)]) 
  
  # Plot the color key on the lower part
  ColorBar(lut=rbPal(steps), min=datalim[1], max=datalim[2], key=key)
  par(new=TRUE)  

  # Plot the actual heatmap
  plot(c(mymatrix[1,1], mymatrix[length(mymatrix[,1]),2]), c(0,1), col="White", xlab="", ylab="", main=title, axes=F)
  
  yheight <- 1/length(ycenter)
  trackheight <- ycenter[2]
  
  for(i in 3:length(mymatrix[1,])){
    RectPlotter(matrixval=mymatrix[,c(1,2,i)], ycoord=ycenter[length(ycenter)-emptylines-i+4], trackheight=ycenter[2], mycolors=mycolors[,i])
   # text(mymatrix[1,1], ycenter[length(ycenter)-emptylines-i+4], colnames(mymatrix)[i], cex=cexsize, pos=4)
    
    # Plot + when bin is significant
  #  if(length(which(signmatrix[,i]!=0))!=0){
  #    index <- which(signmatrix[,i]!=0)
  #    for(k in index){
  #       #text(mean(mymatrix[k,1:2]), ycenter[length(ycenter)-i+4-emptylines], "+", cex=cexsize)
  #       text(mean(mymatrix[k,1:2]), ycenter[length(ycenter)-i+4-emptylines], "+", cex=1)
  #    }
  #  }
  }
}

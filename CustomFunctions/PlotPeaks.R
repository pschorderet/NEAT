#------------------------------------------------------------
# Libraries
source("~/PepsRscripts/RScripts/PepsFunctions/RectPlotter.R")


#------------------------------------------------------------
# PlotPeaks takes a matrix of peaks of form chr - start - end - val
# Plots peaks. If no value is provided, height of peak = trachheight
PlotPeaks  <- function(peaks=get(peak_chr3R), ypos=ycenter[length(ycenter)], trackheight=trackheight, xpostitle=get(peak_chr3R)[1,2], relativeheight="FALSE", title="Peaks", peakcolor="Black"){
  if(length(peaks[1,])!=4){
    peaks  <- cbind(peaks, trackheight)
  }
  if(relativeheight=="TRUE"){
    peaks[,4] <- (peaks[,4]/max(peaks[,4]))*trackheight
  }
  if(relativeheight=="FALSE"){
    peaks[,4] <- trackheight
  }
  
  for(k in 1:length(peaks[,1])){
    #RectPlotter(matrixval=peaks[k,2:3], ycoord=ypos, trackheight=trackheight, mycolors=peakcolor)
    RectPlotter(matrixval=peaks[k,2:3], ycoord=ypos, trackheight=trackheight, mycolors=peakcolor, relrectheight=peaks[k,4])
  }
  text(xpostitle, ypos, title, cex=0.5, pos=4)
}
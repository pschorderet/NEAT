#------------------------------------------------------------
# PlotIABregions plots the iab regions in colors on line=ycoord

#------------------------------------------------------------
# libraries
source("~/PepsRscripts/RScripts/PepsFunctions/CurvedArrows.R")
                                                        
PlotIABregions <- function(ycoord=0.5, trackheight=0.1, printlabels="TRUE"){
  annotations <- read.table(path_BXCannotation,as.is =T, header=T)
  iabcolors <- c('orange', 'red','lightblue','blue','lightgreen','green2','green3', 'green4','darkgreen')  
  regioncenters <- NULL
  pUbx <- 12559186; pabdA <- 12654643; pAbdB <- 12758754
  promotertargets <- c(pUbx, pUbx, pabdA, pabdA, pAbdB, pAbdB, pAbdB, pAbdB, pAbdB)
  
  # Check if labels are TRUE
  if(printlabels == "TRUE"){print(printlabels)
    for(i in 1:length(annotations[,1])){
      # Plots the regulatory regions
      RectPlotter(annotations[i,2:3], ycoord, trackheight, mycolors=iabcolors[i])
      text(mean(as.numeric(annotations[i,2:3])), ycoord, annotations[i,4], cex=0.5, srt=0, adj=0.5)
      
      # Plot the arrows
      regioncenters <- c(regioncenters,mean(as.numeric(annotations[i,2:3])))                                              
      CurvedArrows(x0=regioncenters[i], y0=(ycoord+(trackheight/2)), x1=promotertargets[i], y1=(ycoord+(trackheight/2)), 
                   arrowcolors=iabcolors[i], curveheight=(trackheight/30000), 
                   allabove=TRUE,  arrowheadwidth=1.5, arrowheadlength=0.2)
      
    } 
  }
  # Check if labels are TRUE
  if(printlabels == "FALSE"){print(printlabels) 
    for(i in 1:length(annotations[,1])){
      # Plots the regulatory regions
      RectPlotter(annotations[i,2:3], ycoord, trackheight, mycolors=iabcolors[i])     
      
      # Plot the arrows
      regioncenters <- c(regioncenters,mean(as.numeric(annotations[i,2:3])))                                              
      CurvedArrows(x0=regioncenters[i], y0=(ycoord+(trackheight/2)), x1=promotertargets[i], y1=(ycoord+(trackheight/2)), 
                   arrowcolors=iabcolors[i], curveheight=(trackheight/30000), 
                   allabove=TRUE,  arrowheadwidth=1.5, arrowheadlength=0.2)
      
    } 
  }
}
#------------------------------------------------------------
# RectPlotter
RectPlotter <- function(matrixval, ycoord, trackheight, mycolors="Black", relrectheight=1){
  
  # Plot peaks on entire trackheight
  if(relrectheight == 1){
    rect(matrixval[,1], (ycoord-(0.4*trackheight)), matrixval[,2], (ycoord+(0.4*trackheight)), col=mycolors, border = NA)    
  }
  
  # Plot peaks using specified height
  if(relrectheight != 1){
    rect(matrixval[,1], (ycoord-(0.4*trackheight)), matrixval[,2], ((ycoord-(0.4*trackheight))+relrectheight), col=mycolors, border = NA)    
  }
  
}
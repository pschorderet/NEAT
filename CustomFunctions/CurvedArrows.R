#------------------------------------------------------------
# CurvedArrows plots curved arrows

#------------------------------------------------------------
# libraries
library(igraph)
iArrows <- igraph:::igraph.Arrows

CurvedArrows <- function(x0, y0, x1, y1, arrowcolors='Black', curveheight=0.01, allabove=FALSE, allbelow=FALSE, arrowheadwidth=1, arrowheadlength=0.5){
  posfrac <- 1
  # If allpos==TRUE
  if(allabove==TRUE){
    # If the starting point is higher than the ending point, rectify the curve by multiplying by -1 (puts all to above)
    if(x0>x1){posfrac  <- -1}  
  }
  if(allbelow==TRUE){
    # If the starting point is higher than the ending point, rectify the curve by multiplying by -1 (puts all to above)
    if(x1>x0){posfrac  <- -1}  
  }
  iArrows(x0, y0, x1, y1,
          h.lwd=1, sh.lwd=1, sh.col=arrowcolors,
          curve=curveheight*posfrac, width=arrowheadwidth, size=arrowheadlength)
}
#------------------------------------------------------------
# OutlierValueInMatrix takes a matrix as input and calculates which of the values are outliers, either bycolumn (default) or by byrow
# Returns a similar dimensional matrix of 0s where outliers have original value
OutlierValueInMatrix  <- function(mymatrix, byrow=FALSE, IQR=3){
  # By default, calculate soutlier values by column
  
  
  if(byrow==TRUE){
    bp <- boxplot(t(mymatrix), col="grey", outcol="red", range=IQR, ces=0.5, las=2)
  }else{bp <- boxplot(mymatrix, col="grey", outcol="red", range=IQR, ces=0.5, las=2)}
  # Finds index of outliers and stores them in significant
  significant <- NULL
  for(k in 1:length(bp$out)){
    significant <- c(significant,which(mymatrix==bp$out[k]))
  }
  significant <- sort(significant)
  # Builds a matrix of 0s in which only outlier values are kept
  signmatrix <- matrix(0, nrow=dim(mymatrix)[1], ncol=dim(mymatrix)[2])
  signmatrix[significant]  <- mymatrix[significant]
  colnames(signmatrix) <- colnames(mymatrix)
  return(signmatrix)
}
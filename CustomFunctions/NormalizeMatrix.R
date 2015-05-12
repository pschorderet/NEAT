#------------------------------------------------------------
# NormalizeMatrix takes as input a matrix and the index of the column to normalize to
# returns the normalized version of the matrix
NormalizeMatrix  <- function(mymatrix, controlcolumnindex){
  mymatrixnorm <- (mymatrix/mymatrix[,controlcolumnindex])
  return(mymatrixnorm)
}

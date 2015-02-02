
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  LoadMartGRanges                               *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# LoadMartGRanges takes path to bed file (Mart object)
# return: GRanges object of Mart object
#         


LoadMartGRanges <- function(path2MartObject, binNumber, strand="+"){
  
  mymart <- read.table(path2MartObject, sep="\t", colClasses = "character", header=TRUE)
  head(mymart)
  cat(" \n There are ", dim(mymart[duplicated(mymart$mgi_symbol),])[1], " duplications in ", nameOfBed, ".", sep="")
  newMart <- RawMart2FormattedMart(mymart=mymart, strand=strand)
  
  # Store the GRanges object returned by Mart2GRanges in GRfeat
  binNumber <- binNumber
  GRfeat <- Mart2GRanges(martReply=newMart, binNumber=binNumber)
  head(width(GRfeat))
  # Return GRfeat
  return(GRfeat)
  
}


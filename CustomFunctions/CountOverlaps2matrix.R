
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  CountOverlaps2matrix                          *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# CountOverlaps2matrix takes 2 GRanges objects and computes their overlap. 
#          GRanges 1 = Features
#          GRanges 2 = Reads from experiment
# return: matrix of overlap

CountOverlaps2matrix <- function(GRanges1, GRanges2, normToLibrarySize=TRUE, binNumber=0, geneNamesSingle=0, centered=FALSE){  
  
  # Manage warnings
  suppressWarnings(warning("CountOverlaps2matrix"))
  
  # Text       
  cat(" \n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ", sep="")
        cat(" \n Counting overlaps", sep="")
  # Count overlaps
  tabChIP <- countOverlaps(GRanges1, GRanges2)
  
  # Create matrix of nrow=number of features and ncols=binNumber
        cat(" \n Creating count table \t\t\t(", length(geneNamesSingle), " features x ", binNumber, " bins)", sep="")
  tabChIP <- matrix(tabChIP, nrow=length(geneNamesSingle), ncol=binNumber, byrow=TRUE)  
  dim(tabChIP)
  head(tabChIP)
  which(tabChIP=='NaN'); which(tabChIP=='Inf')
  
  # Set parameters for row names
  if(geneNamesSingle[1]!=0){
    rownames(tabChIP)<- geneNamesSingle
  }
  # Text
        #cat(" \n GRanges of ", length(GRanges1), " on the ... strand", sep="")
  # Don't normalize for feature length if it's constant to save computational time
  if(length(unique(width(GRanges1)))==1){
      cat(" \n Did not normalize to features length \t\t(fixed at ", width(GRanges1[1]), "bp)", sep="")
      # Set parameters for column names
      colnames(tabChIP) <- round(seq(from=-0.5*binNumber*width(ranges(GRanges1)[1]), to=0.5*binNumber*width(ranges(GRanges1)[1]), length.out=binNumber),0)
  }
  # Normalize for the absolute length of each feature (per row) (it could be different for each gene, enhancer, etc.)
  if(length(unique(width(GRanges1)))!=1){
      cat(" \n Normalizing to features lengths (variable) of mean ", floor(mean(width(GRanges1))), " bp", sep="")  
     # WATCH OUT!! Need to transpose for the division to be correct!!
      tabChIP <- t(t(tabChIP)/width(GRanges1))
  }

  # Normalize by library size
  if(normToLibrarySize==TRUE){
        cat(" \n Normalizing to library size \t\t\t\t(", length(GRanges2), " reads)", sep="")
    tabChIP <- tabChIP/length(GRanges2)
  }
  # Manage 'NaN'
  cat(" \n Replacing by 0: \t\t\t\t\t(", length(which(tabChIP=='NaN')), " 'NaN'", sep="")
  if(as.numeric(length(which(tabChIP=='NaN')))>0){
    tabChIP[which(tabChIP=='NaN')] <- 0  
  }
  # Manage 'Inf'
  cat( "  -  ", length(which(tabChIP=='Inf')), " 'Inf')", sep="")
  if(as.numeric(length(which(tabChIP=='Inf')))>0){
    tabChIP[which(tabChIP=='Inf')] <- 0  
  }
  # Text
  cat(" \n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ", sep="")
  # Return tabChIP
  return(tabChIP)
}

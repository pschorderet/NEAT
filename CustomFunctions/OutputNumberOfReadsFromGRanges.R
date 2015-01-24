
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  OutputNumberOfReadsFromGRanges                *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# OutputNumberOfReadsFromGRanges takes a vector of GRanges objects and returns the number of reads in each one
# return: Outputs the number of reads in each GRanges object

OutputNumberOfReadsFromGRanges <- function(myGRsamples){
  cat(" \n **********************************************************", sep="")
  for(i in 1:length(myGRsamples)){
    # Text        
    cat(" \n Reads in ", myGRsamples[i], " :\t", length(get(myGRsamples[i])), sep="")
  }
  cat(" \n **********************************************************\n", sep="")
}
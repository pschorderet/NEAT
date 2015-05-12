
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  CountOverlaps2GRanges                         *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# CountOverlaps2GRanges takes 2 GRanges objects and computes their overlap. 
#          GRanges 1 = tiled genome (bins)
#          GRanges 2 = Reads from experiment
# return: GRanges

CountOverlaps2GRanges <- function(GRanges1, GRanges2, normFactRX, normToLibrarySize=FALSE){  
  
  # Text       
  cat(" \n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ", sep="")
        cat(" \n Counting overlaps", sep="")
  # Count overlaps
  tabChIP <- countOverlaps(GRanges1, GRanges2)
  # Add a metadata column counting number of counts per bin and normalized (over RX) counts per bin
  gr <- GRanges1

  #   length(which(tabChIP!=0))
  #   length(GRanges1)
  #   length(GRanges2)
  #   which(tabChIP!=0)

  gr$RawValue <- tabChIP
  gr$ValueRX <- tabChIP*normFactRX
  
	if(normToLibrarySize == "TRUE"){
		gr$ValueRX <- gr$ValueRX/length(GRanges2)
	}

  return(gr)

}

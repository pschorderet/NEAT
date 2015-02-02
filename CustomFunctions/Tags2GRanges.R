
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  Mart2GRanges		                          *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# Tags2GRanges takes a tag object (array of list)
# return: GRanges object of tag coordinates

# Load dependencies
library(GenomicRanges);


Tags2GRanges <- function(tags){
  tags.gr <- lapply(names(tags),function(n){
    pos <- tags[[n]]
    GRanges(seqnames=n, range=IRanges(start=pos,end=pos), strand="*")
  })
  tags.gr[1]
  # Concatenate
  return(do.call(c,tags.gr))
}

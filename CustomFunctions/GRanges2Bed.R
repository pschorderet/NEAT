

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  GRanges2bed                                   *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# Bam2GRanges takes path to a bam file
# return: Saves a GRanges object in samplespath folder
#         named samplespath.GRanges.RData countaining read starts




GRanges2Bed <- function(GRangesSample, path2bed){
  
  # Header
  
  
  # Format GRanges to bed
  GRbed <- data.frame(seqnames=seqnames(GRangesSample),
                      starts=start(GRangesSample),
                      ends=end(GRangesSample),
                      names=c(rep(".", length(GRangesSample))),
                      strands=strand(GRangesSample))
  GRbed
  head(GRbed)
  write.table(GRbed, file=path2bed, quote=F, sep="\t", row.names=F, col.names=F)

}

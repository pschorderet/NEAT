
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  Bed2GRanges                                   *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# Bed2GRanges takes a bed object (data.frame) of form: 
#            bed$BedChr - bed$BedStart - bed$BedEnd - bed$BedChr
# return: GRanges object

# Load dependencies
library(GenomicRanges);

Bed2GRanges <- function(bedReply, chromosomes=0){
  
  # Create a list of all starts and ends coordinates by chromosome
  chr <- tapply(bedReply$BedChr, as.vector(bedReply$BedChr), function(x) x)  
  starts <- tapply(bedReply$BedStart, as.vector(bedReply$BedChr), function(x) x)
  ends <- tapply(bedReply$BedEnd, as.vector(bedReply$BedChr), function(x) x)
  strand <- tapply(bedReply$BedStrand, as.vector(bedReply$BedChr), function(x) x)  
  
  # Set chromosome names if not supplied
  if(chromosomes[1]==0){
    chromosomes = c(paste("chr", seq(1,19),  sep = ""), "chrX", "chrY")
  }
  # Filter for correct chromosome names only
  chrfiltered <- chr[chromosomes]
  startsfiltered <- starts[chromosomes]
  endsfiltered <- ends[chromosomes]
  genenamesfiltered <- strand[chromosomes]
  
  chrvector <- startsvector <- endsvector <- geneStrand <- NULL
  for(i in 1:length(chromosomes)){
    chrvector <- c(chrvector, rep(names(startsfiltered[i]), length(startsfiltered[[i]])))
    startsvector <- c(startsvector, as.numeric(startsfiltered[[i]]))
    endsvector <- c(endsvector, as.numeric(endsfiltered[[i]]))
    geneStrand <- c(geneStrand, genenamesfiltered[[i]])
  }
  #GR <- GRanges(seqnames=bedReply$BedChr, ranges=IRanges(start=bedReply$BedStart, end=bedReply$BedEnd))
  length(chrvector)
  length(startsvector)
  length(startsvector)
  GRObject <- GRanges(seqnames=chrvector, range=IRanges(startsvector, endsvector))  
  
  # Return
  return(GRObject)
}

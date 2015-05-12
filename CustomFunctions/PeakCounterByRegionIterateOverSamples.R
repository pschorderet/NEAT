#------------------------------------------------------------
# Source functions and packages
source("~/PepsRscripts/RScripts/PepsFunctions/PeakCounterByRegion.R")
#------------------------------------------------------------
# PeakCounterByRegionIterator takes the names of bedgraph as input as well as regions
# get(samples[1]) should return a bedgraph
# returns a matrix of reads by region normalized
PeakCounterByRegionIterateOverSamples <- function(samples, regions, raw=FALSE){
  
  if(raw==TRUE){
    # This doesn't work because iteratively because all samples
    #  don't necessarily have the same number of 'bins'
    
    # regions <- PeakCounterByRegion(bedgraph=get(samples[1]), regions, raw=raw)[,1:2]
  }
  
  totfract <- regions
  for(i in samples){
    name <- paste("f",i, sep="")
    assign(name, PeakCounterByRegion(bedgraph=get(i), regions=regions, raw=raw)[,3])
    
    totfract <- cbind(totfract, get(name))
    colnames(totfract)[length(totfract[1,])] <- name
  }
  return(totfract)
}

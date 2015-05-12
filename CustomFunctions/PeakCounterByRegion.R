#------------------------------------------------------------
# Count number of reads in each region and normalize to total counts
# Requires a bedgraph as input (use ReadBedgraph function) and regions
# Returns a countdata matrix
PeakCounterByRegion <- function(bedgraph, regions, raw=FALSE){
  
  # If raw=TRUE, return raw data
  if(raw==TRUE){
    # This doesn't work because iteratively because all samples
    #  don't necessarily have the same number of 'bins'
    regionstart <- bedgraph[1,2]
    regionend <- bedgraph[length(bedgraph[,4]),3]
    bedgraph[,4] = bedgraph[,4]/sum(as.numeric(bedgraph[,4]))
    countdata <- bedgraph
    return(countdata)
  }
  
  # Calculate total peaks in region
  regionstart <- head(which(bedgraph[,2]>regions[1,1]))[1]
  regionend <- tail(which(bedgraph[,3]<=regions[length(regions[,1]),2]))[6]
  totalpeaksinregion <- sum(bedgraph[regionstart:regionend,4]) 
  totalpeaksinregion
  # Calculate relative fraction of peaks in each region
  rawcounts <- ratiocounts <- logcounts <- countsin <- NULL
  sumval <- 0
  for(i in 1:length(regions[,1])){
    start <- head(which(bedgraph[,2]>regions[i,1]))[1]
    end <- tail(which(bedgraph[,3]<=regions[i,2]))[6]
    #print(which(regions[i,1]<bedgraph[,2] & bedgraph[,3]<regions[i,2]))
    #print(paste(start, end, length(bedgraph[start:end,4]), sep= " - "))
    #countsin <- rbind(countsin, sum(bedgraph[start:end,4]))
    ratiocounts <- rbind(ratiocounts, sum(bedgraph[start:end,4])/totalpeaksinregion)  
    rawcounts <- rbind(rawcounts, sum(bedgraph[start:end,4]))
    logcounts <- rbind(logcounts, log(sum(bedgraph[start:end,4]),2))
  }
  #countdata <- cbind(regions,rawcounts)
  countdata <- cbind(regions,ratiocounts)
  return(countdata)
}


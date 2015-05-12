
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  Bam2GRanges                                   *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# Bam2GRanges takes path to a bam file
# return: Saves a GRanges object in samplespath folder
#         named samplespath.GRanges.RData countaining read starts

# Load dependencies
library(GenomicRanges);
# source Mart2GRanges.R

Bam2GRangesRData <- function(path2bam, chromosomes){
  
  if(DoesTheFileExist(path2bam)==FALSE){    
    stop(cat(" \n", path2bam, " could not be located", sep=""))
  }
  # Text 
      cat(" \n\t Reading .bam file \n    ", path2bam, sep="")
  # Read .bam file
  bamObj <- readGAlignmentsFromBam(path2bam)
  #table(as.vector(seqnames(bamObj)))
  
  # Create tag objects containing the start sit of each raw read
  #     (for all vectors in seqnames(bam), do x (-> grab them))
      cat(" \n\t Creating tags from bam", sep="")     
  tagsObj <- tapply(start(bamObj), as.vector(seqnames(bamObj)), function(x) x)
  # names(tagsObj)
  
  # Filter out random chromosomes
  if(chromosomes[1]==0){
    chromosomes = c(paste("chr", seq(1,19),  sep = ""), "chrX", "chrY")
  }
  #tagsObj <- tagsObj[chromosomes]
  
      cat(" \n\t Formatting tags to GRanges ", sep="")
  # This gives warnings()
  GRangesObj <- Tags2GRanges(tags=tagsObj)
  #table(as.vector(seqnames(GRangesObj)))
      cat(" \n\t Number of reads in GRanges         * ", length(GRangesObj), " *", sep="")
  # Save R object for further usage
      cat(" \n\t Saving GRanges object to:  \n\t    ", paste(sub("bam", "GRangesRData", path2bam), ".GRanges.RData", sep=""), sep="")
      cat(" \n\t Access GRanges object using load(GRangesObj)", sep="")    
      cat(" \n ----------------------------------------------------------", sep="")
  save(GRangesObj, file = paste(sub("bam", "GRangesRData", path2bam), ".GRanges.RData", sep=""))
  # #------------------------------------------------------------  
}





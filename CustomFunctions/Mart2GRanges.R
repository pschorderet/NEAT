
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  Mart2GRanges                                  *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# Mart2GRanges takes a Mart object (data.frame) of form: 
#            mart$MartChr - mart$MartStart - mart$MartEnd - mart$MartGeneName
# return: returned="geneNames" = GRanges object of intervals using start and end coordinates
# return: returned="GR" = Vector countaining geneNames

# Load dependencies
library(GenomicRanges);

Mart2GRanges <- function(martReply, chromosomes=0, binNumber){

    # Create a list of all TSSs coordinates by chromosome
    starts <- tapply(martReply$MartStart, as.vector(martReply$MartChr), function(x) x)
    ends <- tapply(martReply$MartEnd, as.vector(martReply$MartChr), function(x) x)
    genenames <- tapply(martReply$MartGeneName, as.vector(martReply$MartChr), function(x) x)  
    
    # Set chromosome names if not supplied
    if(chromosomes[1]==0){
      chromosomes = c(paste("chr", seq(1,19),  sep = ""), "chrX", "chrY")
    }
    # Filter for correct chromosome names only
    startsfiltered <- starts[names(starts)%in%chromosomes]
    endsfiltered <- ends[names(starts)%in%chromosomes]
    genenamesfiltered <- genenames[names(starts)%in%chromosomes]
    length(unlist(genenamesfiltered))
    # chrvector = Create a vector countaining 'chr...' of length equal to the number of TSS x number of bins per TSS. This is needed to build the GRanges objects later.
    # startsvector = Create a vector countaining the starts of the genes
    # endsvector = Create a vector countaining the ends of the genes
    # geneNames = Store geneNames for further use
    chrvector <- startsvector <- endsvector <- geneNames <- NULL
    for(i in 1:length(chromosomes)){
      chrvector <- c(chrvector, rep(names(startsfiltered[i]), length(startsfiltered[[i]])*binNumber))
      startsvector <- c(startsvector, as.numeric(startsfiltered[[i]]))
      endsvector <- c(endsvector, as.numeric(endsfiltered[[i]]))
      geneNames <- c(geneNames, genenamesfiltered[[i]])
    }
    length(startsvector)
    length(geneNames)
    
    if(binNumber != 1){
      # startsbinsStarts = For each row containing window start, fill in starts for all bins
      startsbinsStarts <- matrix(0, nrow=length(startsvector), ncol=(binNumber+1))
      startsbinsStarts[,1] <- startsvector
      startsbinsStarts[,ncol(startsbinsStarts)] <- endsvector
      # Define the window size for each row because each row has a potential to be different (f.ex. sizes of genes)
      binSizes <- (startsbinsStarts[,ncol(startsbinsStarts)] - startsbinsStarts[,1])/(binNumber)    
      head(binSizes,10)
      # Fill in the matrix column by column
      
      for(j in 2:binNumber){
        startsbinsStarts[,j] <- startsbinsStarts[,j-1] + binSizes
      }
    
      rownames(startsbinsStarts) <- geneNames
      
      # TSSbinsStarts / TSSbinsEnds = Store starts and ends into vector
      binsEnds <- as.vector(t(startsbinsStarts[,-1]-1))
      binsStarts <- as.vector(t(startsbinsStarts[,-length(startsbinsStarts[1,])]))
      # GRObject = Create a GRange object containing all bins of all TSS
      mystrand <- '*'
      if(martReply$Strand[1]=="-1"){mystrand <- '-'}
      if(martReply$Strand[1]=="1"){mystrand <- '+'}
      unique(chrvector)
      GRObject <- GRanges(seqnames=chrvector, range=IRanges(binsStarts, binsEnds), strand=mystrand)
      
      # Add metadata caontaining name of TSS each bin is overlapping
      geneNamesVector <- rep(geneNames, each=binNumber)
      elementMetadata(GRObject)$GeneNames <- geneNamesVector
      GRObject
      # Return
      return(GRObject)
    }
    
    # If binsize = 1, return start and end
    if(binNumber == 1){
      startsbinsStarts <- matrix(0, nrow=length(startsvector), ncol=(binNumber+1))
      startsbinsStarts[,1] <- startsvector
      startsbinsStarts[,ncol(startsbinsStarts)] <- endsvector
      rownames(startsbinsStarts) <- geneNames
      dim(startsbinsStarts)
      # TSSbinsStarts / TSSbinsEnds = Store starts and ends into vector
      binsEnds <- as.vector(t(startsbinsStarts[,-1]-1))
      binsStarts <- as.vector(t(startsbinsStarts[,-length(startsbinsStarts[1,])]))
      # GRObject = Create a GRange object containing all bins of all TSS
      mystrand <- '*'
      if(martReply$Strand[1]=="-1"){mystrand <- '-'}
      if(martReply$Strand[1]=="1"){mystrand <- '+'}
      
      length(chrvector)
      length(binsEnds)
      
      GRObject <- GRanges(seqnames=chrvector, range=IRanges(binsStarts, binsEnds), strand=mystrand)
      
      # Add metadata caontaining name of TSS each bin is overlapping
      geneNamesVector <- rep(geneNames, each=binNumber)
      elementMetadata(GRObject)$GeneNames <- geneNamesVector
      GRObject
      # Return
      return(GRObject)
    }
}

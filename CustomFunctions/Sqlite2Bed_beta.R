
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  Sqlite2Bed                                    *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# Sqlite2Bed takes a path to an sqlite and creates a bed 
#            dictionary containing gene symbols and gene names
# return: Stores annotations in bed file

#------------------------------------------------------------
# Dependencies
library(AnnotationDbi); library(Rsamtools); library(GenomicRanges); library(GenomicAlignments); 
source(paste(path2PepsRScripts, "DoesTheFileExist.R", sep=""))
source(paste(path2PepsRScripts, "ErrorOutput.R", sep=""))

#------------------------------------------------------------
# Function
Sqlite2Bed <- function(path2sqlitefolder, nameOfSqlite){
  
  cat(" \n\n **********************************************************", sep="")  
  cat(" \n *", sep="")
  cat(" \n *\tLoading Sqlite database for \t", nameOfSqlite, sep="")
  cat(" \n *\tPlease be patient while dictionary is being created.", sep="")
  cat(" \n *", sep="")
  cat(" \n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
  
  
  # Check for the existence of the corresponding .sqlite file
  if(DoesTheFileExist(path2file=paste(path2sqlitefolder, nameOfSqlite, ".sqlite", sep=""))!=TRUE){
    ErrorOutput(paste("The ", refGenome, ".sqlite file could not be found", sep=""))
  }  
  
  # Load refGene .sqlite file into GRanges object
  refGene <- loadDb(paste(path2sqlitefolder, nameOfSqlite, ".sqlite", sep=""))
  tx <- transcripts(refGene)
  txn <- values(tx)[,"tx_name"]
  # intronRanges <- intronsByTranscript(refGene, use.names=TRUE)
  # exonRanges <- exonsBy(refGene, "tx", use.names=TRUE)
  # length(exonRanges)
  #eByg <- exonsBy(refGene, by=c("gene"))
  #refGene
  #length(eByg)
    
  library(TxDb.Mmusculus.UCSC.mm9.knownGene)
  refGene <- TxDb.Mmusculus.UCSC.mm9.knownGene
  refGene
  tx <- transcripts(refGene)
  eByg <- exonsBy(refGene, by=c("gene"))
  length(eByg)
  eByg
  names(eByg)
  

  
  #---------------------------------------------------------------------------
  #----------------      This is a better way                 ----------------
  #---------------------------------------------------------------------------
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  eByg <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")
  symb <- select(org.Mm.eg.db, names(eByg), "SYMBOL")
  head(symb)
  sum(is.na(symb[,2]))
  any(duplicated(symb[,2]))
  all.equal(names(eByg), symb[,1])
  names(eByg) <- as.character(symb[,2])
  eByg
  #---------------------------------------------------------------------------
  #---------------------------------------------------------------------------
  #---------------------------------------------------------------------------
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # If the library bed file (geneDict.) does not exists, create it
  if(DoesTheFileExist(path2file=paste(path2sqlitefolder, nameOfSqlite, ".sqlite.bed", sep=""))!=TRUE){
    
    if(DoesTheFileExist(paste(path2sqlitefolder, nameOfSqlite, ".txt", sep=""))!=TRUE){
      ErrorOutput(paste("The ", refGenome, ".txt file could not be found", sep=""))  
    }
    
    cat(" \n\n -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+", sep="")  
    cat(" \n *", sep="")
    cat(" \n *\tCreating gene name to gene symbol dictionary for \t", nameOfSqlite, sep="")
    cat(" \n *\tPlease be patient while dictionary is being created.", sep="")
    cat(" \n *", sep="")
    cat(" \n\n -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+", sep="")  
    
    # Link Gene symbols to gene names and store in GS
    nameOfCorrespondingtxt <- 
    UC_Gene_Symbol <- read.table(paste(path2sqlitefolder, nameOfCorrespondingtxt, sep=""), sep="\t",colClasses = "character")
    head(exonRanges)
    head(exRanges)
    head(names(exonRanges))
    GeneSymbols <- unlist(names(exonRanges))
    GeneNames <- character(length=55419)
    for (i in 1:55419){ 
      x2 <- which(UC_Gene_Symbol$V1 == GeneSymbols[i])
      if(length(x2) > 0){
        GeneNames[i] <- UC_Gene_Symbol$V2[x2]
      }
    }  
    head(GeneNames)
    # Get rid of spaces in case they exists
    GeneNames <- gsub(" ", "", GeneNames)
    geneDict <- cbind(GeneSymbols, GeneNames)
    head(geneDict)
    head(geneDict)
    # Save count table to folder CountTables
    write.table(geneDict, paste(path2sqlitefolder, nameOfSqlite, ".sqlite.bed", sep=""), sep = "\t", row.names = FALSE, col.names=TRUE, quote=FALSE, na="")  
  }# End of if
  
  # Replace gene symbols by gene names
  geneDict <- read.table(paste(path2sqlitefolder, nameOfSqlite, ".sqlite.bed", sep=""), header=TRUE)
  names(exonRanges) <- geneDict$GeneNames
  
  cat(" \n *", sep="")
  cat(" \n *\tThe database ", nameOfSqlite, ".sqlite has been properly loaded", sep="")
  cat(" \n *", sep="")
  cat(" \n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")

  return(exonRanges)
  
}

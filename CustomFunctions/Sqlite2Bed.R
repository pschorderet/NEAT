
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
Sqlite2Bed <- function(TaxonDatabaseKG, TaxonDatabaseDict){
  
  cat(" \n\n **********************************************************", sep="")  
  cat(" \n *", sep="")
  cat(" \n *\tLoading database for \t", TaxonDatabaseKG, sep="")
  cat(" \n *", sep="")
  cat(" \n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
  
  # Load database
  #---------------------------------------------------------------------
  library(TaxonDatabaseKG, character.only = TRUE)    # TxDb.Mmusculus.UCSC.mm10.knownGene
  library(TaxonDatabaseDict, character.only = TRUE)  # org.Mm.eg.db
  
  exonRanges <- exonsBy(get(TaxonDatabaseKG), by = "gene")
  symb <- select(get(TaxonDatabaseDict), names(exonRanges), "SYMBOL")
#   head(symb)
#   sum(is.na(symb[,2]))
#   any(duplicated(symb[,2]))
#   all.equal(names(exonRanges), symb[,1])
  names(exonRanges) <- as.character(symb[,2])
  exonRanges
  
  return(exonRanges)
  
}

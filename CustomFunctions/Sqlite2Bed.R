
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

#------------------------------------------------------------
# Function
Sqlite2Bed <- function(TaxonDatabaseKG, TaxonDatabaseDict, Key){
  
  # Load database
  #---------------------------------------------------------------------
  library(TaxonDatabaseKG, character.only = TRUE)    # TxDb.Mmusculus.UCSC.mm10.knownGene
  library(TaxonDatabaseDict, character.only = TRUE)  # org.Mm.eg.db
  
  exonRanges <- exonsBy(get(TaxonDatabaseKG), by = "gene")
#  symb <- select(get(TaxonDatabaseDict), names(exonRanges), columns = "SYMBOL")
#  symb <- select(get(TaxonDatabaseDict), names(exonRanges), "SYMBOL", "FLYBASE")
#  symb <- select(get(TaxonDatabaseDict), names(exonRanges), "SYMBOL", "FLYBASE")
symb <- select(get(TaxonDatabaseDict), names(exonRanges), "SYMBOL", Key)


#    TaxonDatabaseKG <- "TxDb.Mmusculus.UCSC.mm9.knownGene"
#    TaxonDatabaseDict <- "org.Mm.eg.db"
#    library(TaxonDatabaseKG, character.only = TRUE)
#    library(TaxonDatabaseDict, character.only = TRUE)
#    exonRanges <- exonsBy(get(TaxonDatabaseKG), by = "gene")
#    symb <- select(get(TaxonDatabaseDict), names(exonRanges), columns = "SYMBOL")
#    mouseDB1 <- get("TxDb.Mmusculus.UCSC.mm9.knownGene")
#    keytypes(mouseDB1)
#    help('select')
#    mouseDB2 <- get("org.Mm.eg.db")
#    keytypes(mouseDB)
#    symb <- select(get(TaxonDatabaseDict), names(exonRanges), columns = "SYMBOL")
#   

  #  head(symb)
  #   sum(is.na(symb[,2]))
  #   any(duplicated(symb[,2]))
  #   all.equal(names(exonRanges), symb[,1])
  names(exonRanges) <- as.character(symb[,2])
  exonRanges
  # Finding missing values and renaming them
  dupl <- which(is.na(names(exonRanges)))
  names(exonRanges)[dupl] <- paste("missingValue_", dupl, sep="")
  
  cat(" \n\n **********************************************************", sep="")  
  cat(" \n *", sep="")
  cat(" \n *\tLoading database for \t", TaxonDatabaseKG, sep="")
  cat(" \n *", sep="")
  cat(" \n *\t\t",  length(dupl), " duplications renamed (", length(exonRanges), " genes tot)", sep="")
  cat(" \n *", sep="")
  cat(" \n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
  
  return(exonRanges)
  
}

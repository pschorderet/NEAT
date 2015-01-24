
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  Genes2GOTerm                                  *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# Genes2GOTerm takes a vector of gene names and a db
# return: GO term enrichment
#         

# Load dependencies
library(biomaRt)
library(GO.db)

Genes2GOTerm <- function(geneList, geneUniverse, GOdb="ensembl", orgLib){
    
  # See all databases available
  #listMarts() # To choose BioMart database
  # Choose db
  m <- useMart("GOdb"); listDatasets(m) 
  library(orgLib, character.only = TRUE)
  listAttributes(m) # Choose data types you want to download
  
  getGOLevel <- function(Node="GO:0008150", Children=GOBPCHILDREN, level) {
    for (i in seq_len(level-1)) {
      Node <- mget(Node, Children, ifnotfound=NA)
      Node <- unique(unlist(Node))
      Node <- as.vector(Node)
      Node <- Node[!is.na(Node)]
    }
    return(Node)
  }
  library(GO.db)
  # Selected only part of the hole GO graph; terms that are neither to small nor too big.
  go.lvl4 <- unique(unlist(lapply(4:10,function(i) getGOLevel(level=i)))) 
  head(go.lvl4)
  
  # Select gene universe (all analyzed genes)
  universe <- select(get(TaxonDatabaseDict), names(exonRanges), "SYMBOL")
  head(universe)
 
  # Select DEG
  selected.genes <- geneList[[1]]$genes
  head(selected.genes)
  selectedID <- NULL
  for(k in 1:length(geneList[[1]]$genes)){
    ID <- symb2ID[which(symb2ID$SYMBOL == geneList[[1]]$genes[k]),]
    selectedID <- rbind(selectedID, ID)
  }
  
  selectedID <- selectedID$ENTREZID
  universeID <- universe$ENTREZID
  
  #reading GO terms
  goterms <- unlist(eapply(GOTERM, function(x) x@Term))
  go2genes  <- as.list(org.Mm.egGO2ALLEGS)
  # Use only the IDs of universe
  go2genes2 <- lapply(go2genes, function(x) unique(x[x%in%universeID]) )
  # Select GO terms that have at least 10 genes
  selectedGOterms <- go2genes2[sapply(go2genes2,length)>10 & names(go2genes2)%in%go.lvl4] 
  head(selectedGOterms)
  
  # Iterate through all selectedGOterms 
  sigpval <-  sapply(selectedGOterms,function(xx){
    # Create a booleen vector for each gene: either it is or it is not in a given GO term
    #xx <- selectedGOterms[[]]
    universeGeneIsInGo <- universeID%in%xx
    #which(universeGeneIsInGo==TRUE)
    #universeID[universeGeneIsInGo]
    #   Do this for selected.genes
    selectedGeneIsInGo <- universeID%in%selectedID
    #which(selectedGeneIsInGo==TRUE)
    # Ensure that at least 40% of genes from GO term are in the of selected genes
    # This ensure GO terms that are very small do not become overrepresented all the time
    if(sum(xx%in%selectedID)>.2*length(xx)){
      fisher.test(table(universeGeneIsInGo,InMo),alternative="g")$p #fisher test to test association
    }else {1}
  })  
  quantile(sigpval)
  
  #Bonferroni correction of p-values
  #sigpval2 <- sort(sigpval[sigpval<(0.05/length(selectedGOterms))]) 
  sigpval2 <- sort(sigpval[sigpval<0.05]) 
  # Translate IDs to readable descriptions for the selected GO terms
  enrichedGOterms <-  strtrim(goterms[names(sigpval2)],70) 
  enrichedGOterms
  # Save table
  #write.table(overlap, paste(path2CountTables, nameCT, sep=""), sep = "\t", row.names = FALSE, col.names=TRUE, quote=FALSE, na="")    

  }
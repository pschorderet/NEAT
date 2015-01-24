
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  Features2Mart                                 *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# Features2Mart helps to create bed files using mart annotations
# return: Stores annotations in bed file

#------------------------------------------------------------
# Dependencies
library(biomaRt)

#------------------------------------------------------------
# Define paths to write annotation file to
path2Mart <- "~/PepsRscripts/RScripts/MartObjects/"
nameOfGenes <- "Transcripts"
alignment <- "mm9"

#------------------------------------------------------------
# Get feature from Biomart
# 2014-09-29: De facto db is mm10
# martDB <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# Access mm9
martDB <- useMart('ENSEMBL_MART_ENSEMBL',dataset='mmusculus_gene_ensembl', host="may2012.archive.ensembl.org")
# # See Datasets
# listMarts()
# listDatasets(martDB)
# # Find which attributes you want
attributes <- listAttributes(martDB)
attributes[grep(pattern="ensembl", x=attributes$description, ignore.case=T),]
# head(filters)

# Read in all known genes
# KG <- as.matrix(as.character(unique(martAns$mgi_symbol)))
#KG <- read.table(paste(path2Mart, alignment, "_KnownGenesUnique.bed", sep=""), colClasses = "character", header=TRUE)

# Get transcripts details
martAns <- getBM(attributes=c("chromosome_name", "transcript_start", "transcript_end", "strand", "ensembl_gene_id", "mgi_symbol"), mart=martDB)
# Get rid of duplicates
martAns <- martAns[!duplicated(martAns$mgi_symbol),]
dim(martAns)
length(KG[,1])

head(martAns)

# Fish out specific genes
#martAns <- martAns[grep(pattern=nameOfGenes, x=martAns$mgi_symbol, ignore.case=T),]
#martAns <- martAns[grep(pattern="^hox", x=martAns$mgi_symbol, ignore.case=T),] 
#martAns <-  martAns[which(martAns$mgi_symbol%in%Suz12$Gene==TRUE),]
  

# #------------------------------------------------------------
# # For TSS, use start +/- a given window
# winSize = 10000
# currentind <- martAns$strand==-1
# martAns[currentind,]$transcript_start <- as.numeric(martAns[currentind,]$transcript_end)-(winSize/2)
# martAns[currentind,]$transcript_end <- as.numeric(martAns[currentind,]$transcript_end)+(winSize/2)
# currentind <- martAns$strand==1
# martAns[currentind,]$transcript_end <- as.numeric(martAns[currentind,]$transcript_start)+(winSize/2)
# martAns[currentind,]$transcript_start <- as.numeric(martAns[currentind,]$transcript_start)-(winSize/2)
# head(martAns[currentind,])
# martAns$transcript_end-martAns$transcript_start

#------------------------------------------------------------
# Save feature in bed file
martAns$chromosome_name <- paste("chr", martAns$chromosome_name, sep = "")
head(martAns)
write.table(martAns, paste(path2Mart, paste(alignment, "_", nameOfGenes, ".bed", sep=""), sep=""), sep = "\t", row.names = FALSE, col.names=TRUE, quote=FALSE, na="")





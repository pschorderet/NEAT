
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  GeneList2Mart                                 *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# GeneList2Mart helps to create bed files using mart annotations from a gene list
# return: Stores annotations in bed file

#------------------------------------------------------------
# Dependencies
library(biomaRt)

#------------------------------------------------------------
# Define paths to write annotation file to
path2Mart           <- "~/PepsRscripts/RScripts/MartObjects/"
nameOfGeneListFile  <- "ES_PRC1_targets_SK.txt"
nameOfNewFile       <- "PRC1"
# DB parameters
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
# Get transcripts details
martAns <- getBM(attributes=c("chromosome_name", "transcript_start", "transcript_end", "strand", "ensembl_gene_id", "mgi_symbol"), mart=martDB)
# Get rid of duplicates
martAns <- martAns[!duplicated(martAns$mgi_symbol),]
dim(martAns)

# Read Gene list
geneList <- read.delim(paste(path2Mart, nameOfGeneListFile, sep=""))
geneListGenes <- geneList[,4]
length(geneListGenes)
length(which(duplicated(geneListGenes)))
head(geneListGenes)

# Find intersection between martAns and geneListGenes
martAns <-  martAns[which(martAns$mgi_symbol%in%geneListGenes==TRUE),]
dim(martAns)
head(martAns)


#------------------------------------------------------------
# Save feature in bed file
martAns$chromosome_name <- paste("chr", martAns$chromosome_name, sep = "")
head(martAns)
write.table(martAns, paste(path2Mart, paste(alignment, "_", nameOfNewFile, ".bed", sep=""), sep=""), sep = "\t", row.names = FALSE, col.names=TRUE, quote=FALSE, na="")







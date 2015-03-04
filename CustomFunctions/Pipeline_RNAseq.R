
#-----------------------------------------------------------------
# TEST LINES (DISREGARD)
# path2NEAT='/Users/patrick/Desktop/NEAT/'; path2MainFolder ='/Users/patrick/Desktop/EXAMPLE/'; topNgenes = 100; toHighlight = 10; 
# path2NEAT='/Users/patrick/NEAT/'; path2MainFolder ='~/Documents/Sciences/Kingston/DIPG/DIPG_2014-07-07_RNAseq/'; topNgenes = 100; toHighlight = 10; 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                RNAseq R pipeline for local                     *
#*                                                                *
#******************************************************************
#
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
#                                                                 #
# This script is intended to be run from bash/shell Mac terminal  #
#                                                                 #
# To check if the script works locally or to run the script       #
# through R, users will need to specify the following parameters: #
# parameters first:                                               #
#                                                                 #
#     path2NEAT <- "~/NEAT/"                                      #
#     MainFolder <- "MY_NEW_RNA_PROJECT/"; path2MainFolder <- paste("~/Desktop/", MainFolder, sep="")
#     topNgenes <- 2000; toHighlight <- 20;
#                                                                 #
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #

version <- "1.0.1 Jan 2015"

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*                Load and set parameters                         *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#--------------------------------------------------
# Load R and different packages
library(AnnotationDbi); library(Rsamtools); library(GenomicRanges); library(GenomicAlignments); library(edgeR); library(gtools); library(VennDiagram);

# Define paths to subfolders stored in RScripts
path2bam <- paste(path2MainFolder, "bam/", sep="")
path2GRangesRData <- paste(path2MainFolder, "GRangesRData/", sep="")
path2CountTables <- paste(path2MainFolder, "CountTables/", sep="")
path2Plots <- paste(path2MainFolder, "Plots/", sep="")
path2Logs <- paste(path2MainFolder, "Logs/", sep="")
path2ReferenceFiles <- paste(path2NEAT,"ReferenceFiles/", sep="")
path2CustFct <- paste(path2NEAT,"CustomFunctions/", sep="")
path2Targets <- paste(path2MainFolder, "DataStructure/Targets.txt", sep="")
path2chrlens <- paste(path2MainFolder, "DataStructure/chr_lens.dat", sep="")
path2Tophat <- paste(path2MainFolder, "Tophat/", sep="")

# Source pepsfunctions
source(paste(path2CustFct, "Bam2GRangesRData.R", sep=""))
source(paste(path2CustFct, "CountOverlaps2matrix.R", sep=""))
source(paste(path2CustFct, "DoesTheFileExist.R", sep=""))
source(paste(path2CustFct, "ErrorOutput.R", sep=""))
source(paste(path2CustFct, "LoadMartGRanges.R", sep=""))
source(paste(path2CustFct, "Mart2GRanges.R", sep=""))
source(paste(path2CustFct, "OutputNumberOfReadsFromGRanges.R", sep=""))
source(paste(path2CustFct, "RawMart2FormattedMart.R", sep=""))
source(paste(path2CustFct, "Sqlite2Bed.R", sep=""))
source(paste(path2CustFct, "Tags2GRanges.R", sep=""))

#------------------------------------------------------------
# Redirect output to log file
if(file.exists(path2Logs)==TRUE){cat(" \n Logs file exists.\n\n", sep="")}
if(file.exists(path2Logs)==FALSE){cat(" \n Creating\t", path2Logs, "\n\n", sep="");dir.create(path2Logs)}
sink(paste(path2Logs, Sys.Date(), "_", format(Sys.time(), "%X"), ".txt", sep=""), append=FALSE, split=FALSE)
#sink()

# Find Taxon Db name and dictionary
res <- readLines(path2Targets)
tDb <- res[grep("Proj_TaxonDatabase", res)]; tDb <- gsub("# ","",tDb); tDb <- gsub("\t","",tDb); tDb <- gsub("\"","",tDb)
TaxonDatabaseKG <- unlist(strsplit(tDb, split = "\\="))[2]
tDbDict <- res[grep("Proj_TaxonDatabaseDic", res)]; tDbDict <- gsub("# ","",tDbDict); tDbDict <- gsub("\t","",tDbDict); tDbDict <- gsub("\"","",tDbDict)
TaxonDatabaseDict <- unlist(strsplit(tDbDict, split = "\\="))[2]

# Find if run is SE or PE
SE <- res[grep("Paired_end_seq_run", res)]; SE <- gsub("# ","",SE); SE <- gsub("\t","",SE); SE <- gsub("\"","",SE)
SE <- unlist(strsplit(SE, split = "\\="))[2];
if(SE=="0"){ SE <- "TRUE" }
if(SE=="1"){ SE <- "FALSE" }

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*                Create proper arborescence                      *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Paste variables in log file
cat(" \n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  ", sep="")
cat(" \n ||", sep="")
cat(" \n ||\t RNAseq analysis version ", version, " (© Patrick Schorderet)", sep="")
cat(" \n ||", sep="")
cat(" \n ||\t Path to MainFolder : \t\t\t\t", path2MainFolder, sep="")
cat(" \n ||\t Path to CustomFunctions : \t\t\t", path2CustFct, sep="")
cat(" \n ||", sep="")
cat(" \n ||\t Number of DEG genes : \t\t\t\t", topNgenes, sep="")
cat(" \n ||\t Number of genes to highlight : \t", toHighlight, sep="")
cat(" \n ||\t Taxon database : \t\t\t\t\t", TaxonDatabaseKG, sep="")
cat(" \n ||\t Taxon database dict : \t\t\t\t", TaxonDatabaseDict, sep="")
cat(" \n ||", sep="")
cat(" \n ||\t Sequencing run type : \t\t\t\t", sep="")
if(SE=="TRUE"){ cat("Single-end", sep="") }
if(SE=="FALSE"){ cat("Paired-end", sep="") }
cat(" \n ||", sep="")
cat(paste(" \n ||\t Run date:\t\t\t\t\t\t\t", Sys.Date(), sep=""), sep="")
cat(" \n ||", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ======================================================================================", sep="")

# Create proper arborescence
cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Create / check arborescence of ", paste(path2MainFolder, sep=""), sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")
source(paste(path2CustFct, "RNAseqCreateArborescence.R", sep=""))
RNAseqCreateArborescence(path2MainFolder=path2MainFolder)

#------------------------------------------------------------
# Read the chr_lens.dat
#chromosomesFile <- read.delim(path2chrlens, comment.char="#")
chromosomesFile <- read.table(path2chrlens, comment.char="#")
chromosomes <- chromosomesFile$V1
#chromosomes = c(paste("chr", seq(1,19),  sep = ""), "chrX", "chrY")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                                 #
#                  Read Targets.txt                               #
#                                                                 #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Read Targets file", sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")
# Check if Targets file exists
if(file.exists(path2Targets)==FALSE){ErrorOutput(paste("No Targets.txt file in\t", path2Targets, sep="")) }

cat(" \n\n Targets file provided: \n\n", sep="")
Targets <- read.delim(path2Targets, comment.char="#")
print(Targets)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                                 #
#                  Load GRanges data                              #
#                                                                 #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Load datasets provided in Targets file
samples <- Targets$FileName
allsamples <- levels(Targets$FileName)

#------------------------------------------------------------
# Check existence of .bam.GRanges.RData or .bam files
cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Create / load GRanges object from bam files", sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")

# If .bam.GRanges.RData files DO NOT exists, create them
bamGRangesRDataPath <- paste(paste(path2GRangesRData, allsamples, sep=""), ".bam.GRanges.RData", sep="")
if(DoesTheFileExist(path2file=bamGRangesRDataPath)!=TRUE){
  # If .bam files cannot be found, abort
  bamPath <- paste(path2file=path2bam, allsamples, ".bam", sep="")
  if(DoesTheFileExist(path2file=bamPath)!=TRUE){
    ErrorOutput("Some .bam files could not be found")
  }
  # If .bam files are present, create .bam.GRanges.RData files
  if(DoesTheFileExist(path2file=bamPath)==TRUE){
    # Create GRangesRData files
    cat(" \n\n **********************************************************\n", sep="")
    cat(" \n Generating .bam.GRanges.RData files", sep="")
    for(k in 1:length(bamPath)){    
      # k=1
      cat(" \n\n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
      cat(" \n ", k, " / ", length(bamPath) ," :\t", bamPath[k], sep="")
      Bam2GRangesRData(path2bam=bamPath[k], chromosomes=chromosomes)
    }  
  }
}

# Check if .bam.GRanges.RData files exists, load them
GRangesSamples <- list.files(path2GRangesRData, pattern = "*.bam.GRanges.RData$")  
if(DoesTheFileExist(path2file=paste(paste(path2GRangesRData, allsamples, sep=""), ".bam.GRanges.RData", sep=""))==TRUE){
  # This doesn't work if put in a function!!  
  # Load each sample into its .GRanges named object
  for(i in 1:length(GRangesSamples)){
    path2load <- paste(path2GRangesRData, GRangesSamples[i], sep="")
    load(path2load)
    assign(GRangesSamples[i], GRangesObj)
  }
}

#------------------------------------------------------------
# Quality Control 
path2GRangesSamples <- paste(paste(path2GRangesRData, allsamples, sep=""), ".bam.GRanges.RData", sep="")
OutputNumberOfReadsFromGRanges(GRangesSamples)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                                 #
#                  Load reference files/genomes                   #
#                                                                 #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#------------------------------------------------------------
# Load .sqlite file
exonRanges <- Sqlite2Bed(TaxonDatabaseKG=TaxonDatabaseKG, TaxonDatabaseDict=TaxonDatabaseDict)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                                 #
#                      Count overlaps                             #
#                                                                 #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Compute overlaps for each sample over exons", sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")

# Store active countTable names in GRangesSamplesCT
GRangesSamplesCT <-NULL
prenameForCT <- "CountTable_OverExons_"
for(i in 1:length(GRangesSamples)){
  #i=1
  cat(" \n\n **********************************************************", sep="")  
  cat(" \n *", sep="")
  cat(" \n *\tCompute overlaps over exons for \t", GRangesSamples[i], sep="")
  cat(" \n *", sep="")
  # Compute overlap
  #    Also, summarizeOverlaps deals with GRangesList, so no need to 'consolidate' all exons to genes!
  overlap <- summarizeOverlaps(exonRanges, get(GRangesSamples[i]), mode="Union", ignore.strand=TRUE, inter.feature=FALSE, singleEnd=SE)  
  #table(as.vector(seqnames(get(GRangesSamples[i]))))  
  
  head(head(assays(overlap)$counts))
  assays(overlap)$counts[which(rownames(overlap)=="Ccnd1")]
  
  overlapDF <- as.data.frame(cbind(names(exonRanges), assays(overlap)$counts))
  colnames(overlapDF) <- c("GeneNames", "Counts")
  
  # Assign to nameCT
  nameCT <- paste(prenameForCT, GRangesSamples[i], ".bed", sep="")
  assign(nameCT, overlapDF)
  GRangesSamplesCT <- c(GRangesSamplesCT, nameCT)
  # Save to .bed file
  cat(" \n *\t\tStore countTable in : \t\t", paste(path2CountTables, nameCT, sep=""), sep="")    
  write.table(overlapDF, paste(path2CountTables, nameCT, sep=""), sep = "\t", row.names = FALSE, col.names=TRUE, quote=FALSE, na="")    

  cat(" \n *", sep="")
  cat(" \n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                                 #
#                      Compute DEG                                #
#                                                                 #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Compute RPKMs and consolidate into genes", sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")

# Generate DGEList and calculate RPKM
#CTallsamples <- list.files(path2CountTables, pattern = "*.bam.GRanges.RData.bed$")
CTallsamples <- paste(prenameForCT, allsamples, ".bam.GRanges.RData.bed", sep="")

# Store all countTables into CTall
# Start with first
cat(" \n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n", sep="")
cat(" \n \t Read countTable \t", CTallsamples[1], sep="")
nameCT <- paste(path2CountTables, CTallsamples[1], sep="")
CT <- read.table(nameCT, header=TRUE)
head(CT)
colnames(CT) <- c("GeneNames", CTallsamples[1])
CTall <-CT
# cbind all others
if(length(CTallsamples)>1){
  for(k in 2:length(CTallsamples)){
    # k=2
    cat(" \n \t Read countTable \t", CTallsamples[k], sep="")
    nameCT <- paste(path2CountTables, CTallsamples[k], sep="")
    CT <- read.table(nameCT, header=TRUE)
    CTall <- cbind(CTall, CT$Counts)
    colnames(CTall)[k+1] <- CTallsamples[k]
  }
}
# Save countTable containaing all samples
cat(" \n\n *\t\tStore countTable in : \t\t", paste(path2CountTables, "CountTable_OverExons_AllSamples.bed", sep=""), sep="")    
cat(" \n\n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n", sep="")
write.table(CTall, paste(path2CountTables, prenameForCT, "_AllSamples.bed", sep=""), sep = "\t", row.names = FALSE, col.names=TRUE, quote=FALSE, na="")    

# Format into matrix form
CTall <- as.matrix(CTall[2:(length(CTallsamples)+1)])
rownames(CTall) <- CT[,1]
#--------------------------------------------------
# Compute DGElist
#--------------------------------------------------
cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Compute DGEList and consolidate into genes", sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")

#--------------------------------------------------
# Create a DGEList object countaining number of reads per exon
cat(" \n \t Compute DGEList \t", sep="")
# !!! Deal with replicates here...

# Get library size
lib.size <- NULL
for(i in 1:length(GRangesSamples)){lib.size <- c(lib.size, length(get(GRangesSamples[i])))}

#d <- DGEList(counts=CTall, group=factor(Targets$Replicate), genes=rownames(CTall))
d <- DGEList(counts=CTall, lib.size=lib.size, group=Targets$Replicate, genes=rownames(CTall))
#d <- DGEList(counts=CTall, group=Targets$Replicate, genes=rownames(CTall))
#dim(CTall)
#colnames(CTall)
#--------------------------------------------------
# Estimate the different parameters
cat(" \n \t Compute Normalization factor", sep="")
d <- calcNormFactors(d)

# If there are no replicates
if(length(which(duplicated(Targets$Replicate)==FALSE))!=length(Targets$Replicate)){
  cat(" \n \t Estimate common dispersion", sep="")
  d <- estimateCommonDisp(d, verbose=TRUE)
  cat(" \n \t Estimate tagwise dispersion", sep="")
  d <- estimateTagwiseDisp(d, verbose=TRUE)
  bcv <- d$common.dispersion
  # Save the QC1 plot of common dispersion
  pdf(paste(path2Plots, Sys.Date(), "_QC1.pdf", sep=""), paper='USr', useDingbats=FALSE, width = 27, height = 21) 
  plotBCV(d, main="plotBCV(d)")
  cat("\n\n", sep="")
  dev.off()
}else{
  # Set common dispersion to 0.2 if there are no replicates
  cat(" \n \t\t Common dispersion set to 0.2 (no replicates available)", sep="")  
  bcv <- 0.2
  cat(" \n \t\t No replicates, tagwise dispersion cannot be computed", sep="")
}

#--------------------------------------------------
# Generate count table of RPKM
cat(" \n \t Compute RPKMs", sep="")
numBases <- sum(width(exonRanges))
RPKM <- rpkm(d,numBases)
RPKM <- cbind(CT[,1], RPKM)
colnames(RPKM) <- c("Gene Symbol", GRangesSamples)

#--------------------------------------------------
# Consolidated and save RPKM table
cat(" \n \t Consolidate RPKMs (this can take a few minutes...)", sep="")
RPKM2 <- tapply(1:nrow(RPKM), rownames(RPKM), function(x) {
  if(length(x)==1){ as.data.frame(RPKM)[x,-1]
  }else{colSums(as.data.frame(RPKM)[x,-1]) }})
RPKMf <- data.frame(matrix(unlist(RPKM2), nrow=length(RPKM2), byrow=T))
RPKMf <- cbind(names(RPKM2), RPKMf)
colnames(RPKMf) <- c("GeneNames", GRangesSamples)
# Write RPKM table to CounTables folder
path2CountTable <- paste(path2CountTables, "CountTable_OverGenes_RPKM.bed", sep="")
cat(" \n \t Store RPKMs in ", path2CountTable, sep="")
write.table(RPKMf, path2CountTable, sep = "\t", row.names = TRUE, col.names=TRUE, quote=FALSE)

#--------------------------------------------------
# Read in RPKM and counts data.frames generated above
RPKM = read.delim(path2CountTable, row.names=1)
# Make quality plots and save as pdf
pdf(paste(path2Plots, Sys.Date(), "_RPKM_byGene.pdf", sep=""), paper='USr', useDingbats=FALSE, width = 27, height = 21) 
Q_samples <- log2(RPKM[,-1])
pairs(Q_samples)
pairs(Q_samples, col=rgb(255,0,0,20,maxColorValue=255), pch=16)
cat("\n\n", sep="")
dev.off()

#--------------------------------------------------
# Start analysis DEG using EdgeR
cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Compute DEG using edgeR", sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")

# Open a pdf to store plotSmear plots
pdf(paste(path2Plots, Sys.Date(), "_plotSmears.pdf", sep=""), paper='USr', useDingbats=FALSE, width = 27, height = 21) 

#  ! ! ! ! ! ! ! ! This isn't suited for replicates ! ! ! ! ! ! ! !
# Find all combinations of samples and store in cmp
cmp <- combinations(n=length(d$samples$group), r=2, v=d$samples$group, repeats.allowed=FALSE)
#symb2ID <- select(get(TaxonDatabaseDict), names(exonRanges), "SYMBOL")

# Run exactTest for each comparison
for(k in 1:nrow(cmp)){  

    # Set the comparisons
    comp <- exactTest(d, pair=c(cmp[k,1], cmp[k,2]), dispersion = bcv^2)    
    
    # Set title
    name1 <- as.character(Targets$FileName[which(Targets$Replicate==cmp[k,1])])[1]
    name2 <- as.character(Targets$FileName[which(Targets$Replicate==cmp[k,2])])[1]
    title = paste(name1, "\n vs \n", name2, sep="")
    comp$comparison <- c(name1, name2)
    
    # Compute DEG
    fdr <- 0.05
    toptags <- topTags(comp, n=topNgenes, sort.by="logFC")
    tablediff <- summary(de <- decideTestsDGE(comp, p=fdr, adjust="BH"))
    
    # Save gene list for GO term analysis using GOrilla
    write.table(toptags[,1], paste(path2CountTables, "GOrilla_", name1, "-", name2, ".bed", sep=""), sep = "\t", row.names = FALSE, col.names=FALSE, quote=FALSE)
    
    # Text to logs
    cat(" \n\t- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n", sep="")
    cat(" \n\t Compararison ", k, "\t\t", comp$comparison[2], "\t -vs- \t", comp$comparison[1], "\n", sep="")
    cat(" \n\t\t Statistically significant changes: ", sep="")
    cat(" \n\t\t\t FDR \t\t\t\t\t ", fdr, sep="")
    cat(" \n\t\t\t Downregulated genes \t ", tablediff[1], sep="")
    cat(" \n\t\t\t Unchanged genes \t\t ", tablediff[2], sep="")
    cat(" \n\t\t\t Upregulated genes \t\t ", tablediff[3], sep="")
    cat(" \n\n\t Saving graphs in pdf format \t ", tablediff[3], sep="")
    cat(" \n\n\t Store top ", topNgenes," DEG list in \t ", paste(path2CountTables, "DEG_byExon_", name1, "-", name2, ".bed", sep=""), sep="")
    write.table(toptags, paste(path2CountTables, "DEG_byExon_", name1, "-", name2, ".bed", sep=""), sep = "\t", row.names = TRUE, col.names=TRUE, quote=FALSE)
    
    # plotSmears
        # statistically significant genes plotted in red
    detags <- rownames(d)[as.logical(de)]
    # plotSmear(comp, main=title, cex.main=0.6, col="Green")
    plotSmear(comp, de.tags=detags, main=title, cex.main=0.6, col=rgb(0,255,0,50,maxColorValue=255))
    xcoord <- c((par()$xaxp[1]+(par()$xaxp[3]*0.6)), (par()$xaxp[1]+(par()$xaxp[3]*0.8)))
    ycoord <- c(par()$yaxp[2], par()$yaxp[2])
    
    # Add name of genes
    comp.tags <- topTags(comp, n=topNgenes)$table
    
    #comp.tags
    updown <- comp.tags[1:toHighlight,]
    up <- comp.tags[which(updown$logFC>0), ]
    down <- comp.tags[which(updown$logFC<0), ]
    abline(h = c(-1, 1), col = "lightgrey")        
    
    # Highlight top toHighlight DEG in blue
    if(dim(down)[1]!=0){
      points(x=down$logCPM, y=down$logFC, cex.main=0.6, col="Blue", pch=21)
      text(x=down$logCPM, y=down$logFC, labels=down$genes, cex=0.5, pos=4, srt=-40)
      text(x=xcoord[1], y=ycoord[2], labels=paste(c("Downregulated\n", sort(down$genes)), sep="", collapse="\n"), cex=0.5, adj = c( 0, 1 ))    
    }
    if(dim(up)[1]!=0){
      points(x=up$logCPM, y=up$logFC, cex.main=0.6, col="Blue", pch=21)
      text(x=up$logCPM, y=up$logFC, labels=up$genes, cex=0.5, pos=4, srt=40)
      text(x=xcoord[2], y=ycoord[2], labels=paste(c("Upregulated\n", sort(up$genes)), sep="", collapse="\n"), cex=0.5, adj = c( 0, 1 ))    
    }
    
    # GO term analysis
    #Genes2GOTerm(geneList=toptags, geneUniverse=exonRanges, GOdb="ensembl", orgLib=TaxonDatabaseDict)
    
}
cat("\n\n", sep="")
dev.off() # "_plotSmears.pdf"

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                                 #
#                  Venn diagrams                                  #
#                                                                 #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Create Venn diagrams
cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Creating venn diagrams", sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")

# Load topNgenes for each sample and store in its file name
DEGlist <- list.files(path2CountTables)[grep("DEG_byExon_", list.files(path2CountTables))]
vennTable <- NULL
for(k in 1:length(DEGlist)){
    # k=1
    assign(DEGlist[k], read.delim(paste(path2CountTables, DEGlist[k], sep="")))  
    vennTable <- cbind(vennTable, as.character(get(DEGlist[k])$genes))
}
mycols <- c("Red", "blue", "green", "White", "Purple", "Grey", "Cyan")
vennList <- split(vennTable, rep(1:ncol(vennTable), each = nrow(vennTable)))


# Open a pdf to store plotSmear plots

if(length(vennList)>5){
  toKeep <- 1:5
  cat(" \n\t There are too many elements to compare in a Venn plot", sep="")
  cat(" \n\t Only the 5 first elements will be compared", sep="")
  cat(" \t\t ", DEGlist[toKeep], sep="\n\t\t")
  vennList <- vennList[toKeep]
  vennTable <- vennTable[,toKeep]
  DEGlist <- DEGlist[toKeep]
}

pdf(paste(path2Plots, Sys.Date(), "_Venn.pdf", sep=""), paper='USr', useDingbats=FALSE, width = 27, height = 21) 
venn.plot <- venn.diagram(vennList, NULL, fill=mycols[1:ncol(vennTable)],  alpha=rep(0.3, length(vennList)), cex = 2, cat.fontface=4, category.names=DEGlist, cat.cex = 0.8, main="Common DEG genes")
plot.new()
grid.draw(venn.plot)
dev.off() # "_Venn.pdf"

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                                 #
#                  GO terms analysis                              #
#                                                                 #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Create Venn diagrams
cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t GO terms analysis", sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                                 #
#                  Session info and closure                       #
#                                                                 #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Print sessionInfo()
cat(" \n\n\n || -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-   ", sep="")
cat(" \n\n", sep="")
print(sessionInfo())
cat(" \n\n", sep="")
cat(" \n || -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n\n   ", sep="")

cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  ", sep="")
cat(" \n ||", sep="")
cat(" \n ||\t RNAseq analysis version ", version," (© Patrick Schorderet 2014)", sep="")
cat(" \n ||", sep="")
cat(" \n ||\t Terminated with no known error.", sep="")
cat(" \n ||", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ======================================================================================\n\n", sep="")

# Close the R session when run from bash
quit(save = "no", status = 0, runLast = TRUE)

#}
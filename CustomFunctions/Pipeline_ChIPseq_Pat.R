
#Pipeline_ChIPseq <- function(path2NEAT, path2MainFolder, nameOfBed, binNumber, strand, runmeank, Venn, normInp ){

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                ChIPseq R pipeline for local                    *
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
#     
#     MainFolder <- "MY_NEW_CHIP_PROJECT/"; path2MainFolder <- paste("~/Desktop/", MainFolder, sep="")
#     path2NEAT='/Users/patrick/NEAT/'; path2MainFolder ='~/Documents/Sciences/Kingston/DIPG/DIPG_consolidated_ChIPseq/';
#     nameOfBed <- "mm9_TSS_10kb.bed" ; binNumber = 100 ; strand <- "+" ; Venn <- FALSE ; normInp <- FALSE; runmean=10
#     nameOfBed <- "mm9_Transcripts.bed" ;
#     nameOfBed <- "mm9_Enhancers" ;
#     runmeank <- 5
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
library(Rsamtools); library(GenomicRanges); library(GenomicAlignments); library(caTools);  library(VennDiagram);library(rtracklayer)

# Define paths to subfolders stored in RScripts
path2bam <- paste(path2MainFolder, "bam/", sep="")
path2GRangesRData <- paste(path2MainFolder, "GRangesRData/", sep="")
path2CountTables <- paste(path2MainFolder, "CountTables/", sep="")
path2Plots <- paste(path2MainFolder, "Plots/", sep="")
path2Logs <- paste(path2MainFolder, "Logs/", sep="")
path2Mart <- paste(path2NEAT, "MartObjects/", sep="")
path2ReferenceFiles <- paste(path2NEAT,"ReferenceFiles/", sep="")
path2CustFct <- paste(path2NEAT,"CustomFunctions/", sep="")
path2Targets <- paste(path2MainFolder, "DataStructure/Targets.txt", sep="")
path2wig <- paste(path2MainFolder, "/wig", sep="")

# Source CustomFunctions
source(paste(path2CustFct, "Bam2GRangesRData.R", sep=""))
source(paste(path2CustFct, "Bed2GRanges.R", sep=""))
source(paste(path2CustFct, "CountOverlaps2matrix.R", sep=""))
source(paste(path2CustFct, "CountOverlaps2GRanges.R", sep=""))
source(paste(path2CustFct, "DoesTheFileExist.R", sep=""))
source(paste(path2CustFct, "ErrorOutput.R", sep=""))
source(paste(path2CustFct, "GRanges2Bedwig.R", sep=""))
source(paste(path2CustFct, "LoadMartGRanges.R", sep=""))
source(paste(path2CustFct, "Mart2GRanges.R", sep=""))
source(paste(path2CustFct, "OutputNumberOfReadsFromGRanges.R", sep=""))
source(paste(path2CustFct, "RawMart2FormattedMart.R", sep=""))
source(paste(path2CustFct, "Tags2GRanges.R", sep=""))

#------------------------------------------------------------
# Redirect output to log file

if(file.exists(path2Logs)==TRUE){cat(" \n Logs file exists.\n\n", sep="")}
if(file.exists(path2Logs)==FALSE){cat(" \n Creating\t", path2Logs, "\n\n", sep="");dir.create(path2Logs)}
sink(paste(path2Logs, Sys.Date(), "_", format(Sys.time(), "%X"), ".txt", sep=""), append=FALSE, split=FALSE)
#sink()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*                Create proper arborescence                      *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Paste variables in log file
cat(" \n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  ", sep="")
cat(" \n ||", sep="")
cat(" \n ||\t ChIPseq analysis version ", version, " (© Patrick Schorderet)", sep="")
cat(" \n ||", sep="")
cat(" \n ||\t Path to MainFolder : \t\t", path2MainFolder, sep="")
cat(" \n ||\t Path to CustomFunctions : \t", path2NEAT, sep="")
cat(" \n ||\t Features to align to : \t", nameOfBed, sep="")
cat(" \n ||\t No of bins per feature : \t", binNumber, sep="")
cat(" \n ||\t Smoothing window : \t\t", runmeank, sep="")
cat(" \n ||\t Strand specificity : \t\t", strand, sep="")
cat(" \n ||\t Compute venn diagrams : \t", Venn, sep="")
cat(" \n ||\t Normalize to input : \t\t", normInp, sep="")
cat(" \n ||", sep="")
cat(paste(" \n ||\t Run date:\t\t\t\t\t", Sys.Date(), sep=""), sep="")
cat(" \n ||", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ======================================================================================", sep="")

# Create proper arborescence
cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Create / check arborescence of ", paste(path2MainFolder, sep=""), sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")
source(paste(path2CustFct, "ChIPseqCreateArborescence.R", sep=""))
ChIPseqCreateArborescence(path2MainFolder=path2MainFolder)

#--------------------------------------------------
# Define function to plot axes
PlotAxes <- function(mat, binrange=0, yrange=0, strand="+"){
  if(length(colnames(mat))==0){
    if(strand=="+"){axis(1, at=seq(from=1, to=ncol(mat), length.out=2), labels=c("TSS", "TES"))}
    if(strand=="-"){axis(1, at=seq(from=1, to=ncol(mat), length.out=2), labels=c("TES", "TSS"))}
  }
  if(length(colnames(mat))!=0){
    colnames(mat[binNumber/2])
    
    if(length(binrange)==1){}
    if(length(colnames(mat))!=1){
      myseq <- floor(seq(1, length(colnames(mat)), length.out=4))
      axis(1, at=myseq, labels=colnames(mat)[myseq])    
      abline(v=binrange[length(binrange)/2])
    }  
  }  
  # y axis
  if(length(yrange)>1){axis(2, at=yrange, labels=round(yrange,5))}  
}

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

#------------------------------------------------------------
# Read the Targets.txt file

cat(" \n\n Targets file provided: \n\n", sep="")
Targets <- read.delim(path2Targets, comment.char="#")
print(Targets)

res <- readLines(path2Targets)
for(i in 1:length(res)){
  newLine <- res[i]
  # Store some variables
  if(length(grep("Reference_genome\t", newLine))==1) { 
    currentline <- gsub("# ", "", newLine); currentline <- gsub("\t", "", currentline); currentline <- gsub("\"", "", currentline);
    refGenome <- unlist(strsplit(currentline, split = "\\="))[2]
  }    
}

#------------------------------------------------------------
# Read the chr_lens.dat
#chromosomesFile <- read.delim(path2chrlens, comment.char="#")
path2chrlens <- paste(path2MainFolder, "DataStructure/", refGenome, "/chr_lens.dat", sep="")
chromosomesFile <- read.table(path2chrlens, comment.char="#")
chromosomes <- chromosomesFile$V1
#chromosomes = c(paste("chr", seq(1,19),  sep = ""), "chrX", "chrY")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                                 #
#                  Load GRanges data                              #
#                                                                 #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Load datasets provided in Targets file
samples <- Targets$FileName
allsamples <- c(levels(Targets$FileName), levels(Targets$InpName))
# Remove '-' when no input is provided
if(length(which(allsamples=="-")) != 0){
  allsamples <- allsamples[-which(allsamples=="-")]    
}

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
      #Bam2GRangesRData(path2bam=bamPath[1], chromosomes=0)
    }  
  }
}

#dirname(bamGRangesRDataPath[1])

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
OutputNumberOfReadsFromGRanges(GRangesSamples)



path2bed <- paste(path2bedwig, "/", samples[1], ".bam.bedwig", sep="")

# Find Taxon Db name and dictionary
res <- readLines(path2Targets)
tDb <- res[grep("Proj_TaxonDatabase", res)]; tDb <- gsub("# ","",tDb); tDb <- gsub("\t","",tDb); tDb <- gsub("\"","",tDb)
TaxonDatabaseKG <- unlist(strsplit(tDb, split = "\\="))[2]
tDbDict <- res[grep("Proj_TaxonDatabaseDic", res)]; tDbDict <- gsub("# ","",tDbDict); tDbDict <- gsub("\t","",tDbDict); tDbDict <- gsub("\"","",tDbDict)
TaxonDatabaseDict <- unlist(strsplit(tDbDict, split = "\\="))[2]

TaxonDatabase = "BSgenome.Mmusculus.UCSC.mm9"
library(TaxonDatabase, character.only = TRUE)

#source(paste(path2CustFct, "GRanges2Bedwig.R", sep=""))




source(paste(path2CustFct, "CountOverlaps2GRanges.R", sep=""))
source(paste(path2CustFct, "CheckExistenceOfFolder.R", sep=""))

path2wig <- paste(path2MainFolder, "/wig", sep="")

CheckExistenceOfFolder(path2wig)

# Set variables and paths
path2chrlens <- paste(path2MainFolder, "DataStructure/", refGenome, "/chr_lens.dat", sep="")
chromosomesFile <- read.table(path2chrlens, comment.char="#")
seqlengths <- c(chromosomesFile[,2])
names(seqlengths) <- chromosomesFile[,1]
seqlengths

# Slice the genome into bins
bins <- tileGenome(seqlengths, tilewidth=500, cut.last.tile.in.chrom=TRUE)
bins


# Read in the statTable file countaining the ChIPrx normaization factors
path2statTable <- paste(path2MainFolder, "DataStructure/statTable.bed", sep="")
statTable <- read.table(path2statTable, header=TRUE, comment.char="#")

# Disable scientific notation to ensure no 'e' are found in the wig files
options(scipen=999)
# Set paths to tmp files

path2wigfilename_tmp <- paste(path2wig, "/tmp.wig", sep="")
path2wigfile_header <- paste(path2wig, "/header.wig", sep="")

for(h in 1:length(GRangesSamples)){
  # h=1
  gr <- CountOverlaps2GRanges(GRanges1=bins, GRanges2=get(GRangesSamples[h]), normFactRX=statTable$NormFact[h], normToLibrarySize=FALSE)
  df <- data.frame(seqnames=seqnames(gr), starts=start(gr)-1, ends=end(gr), value=elementMetadata(gr)$ValueRX)
  # Delete line that have no value
  dfnon0 <- df[which(df$value!=0),]
  # Set path to current wig file
  path2wigfilename <- paste(path2wig, "/", GRangesSamples[h], ".wig", sep="")
  # Copy dfnon0 to tmp file
  write.table(dfnon0, file=path2wigfilename_tmp, quote=F, sep="\t", row.names=F, col.names=F)
  # Add header to the file
  header = paste("track type=wiggle_0 name='", allsamples[h], "' description='coverata' visibility=dense color=0,100,200 priority=20", sep="")
  # Copy header to tmp file
  write.table(header, file=path2wigfile_header, quote=F, sep="\t", row.names=F, col.names=F)
  # Consolidate files
  mycode <- paste("`echo '", header, "' > ", path2wigfilename, "`", sep="")
  system(mycode)
  mycode <- paste("`cat ", path2wigfile_header, " ", path2wigfilename_tmp, " > ", path2wigfilename, "`", sep="")
  system(mycode)

}

# Delete temporary files
mycode <- paste("`rm ", path2wigfilename_tmp, "`", sep="")
system(mycode)
mycode <- paste("`rm ", path2wigfile_header, "`", sep="")
system(mycode)


print(statTable)
print(statTable)














#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*                Load reference files/genomes                    *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#------------------------------------------------------------
# Create GRanges objects contaninig bins around all feature
#     Use the Features2Mart script to create this object !!

# Load Mart object & store gene names
cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Create / load Mart file as GRanges object", sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")
# Load Mart object into GRanges
path2MartObject <- paste(path2Mart, nameOfBed, sep="")
GRfeat <- LoadMartGRanges(path2MartObject=path2MartObject, binNumber=binNumber, strand=strand)
geneNames <- GRfeat$GeneNames[seq(from=1, to=length(GRfeat), by=binNumber)]

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*                  Count overlaps                                *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# For each GRanges object:  - Compute overlap and store countTables
#                           - Plot average coverage and store in pdf
cat(" \n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Compute overlaps and generate count tables", sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")

# Open pdf to store the plots
pdf(paste(path2Plots, sub(".bed", "", nameOfBed), "_", strand,"_", runmeank, ".pdf", sep=""), paper='USr', useDingbats=FALSE, width = 27, height = 21) 

# Compute and store countTables, pdf, etc
cols <- c("Black", "Blue", "Red", "Green", "Yellow", "Purple", "Cyan", "Grey")

for(i in levels(Targets$Factor)){
  # i=levels(Targets$Factor)[1]
  cat(" \n\n **********************************************************", sep="")  
  cat(" \n *", sep="")
  cat(" \n *\tFactor: \t", i, sep="")
  cat(" \n *", sep="")
  cat(" \n * * * * * * * * * * * * * * * * * * * * * * * * * * * * *", sep="")
  # i = levels(Targets$Factor)[1]
  sameFactorSamples <- Targets$FileName[which(Targets$Factor==i)]
  
  # For each sample with the same Factor, compute overlap and plot
  activeCTNames <- activeCTNamesRep <- NULL
  for(j in sameFactorSamples){
    # Find which regions of HA are enriched
    # j =  sameFactorSamples[1]
    
    cat(" \n\n Sample: \t", as.character(j), "\n", sep="")    
    
    #--------------------------------------------------------
    #----------------   Replicates   ------------------------
    # If there is a Replicate corresponding to the current sample, average the countTable
    currentDupNum <- Targets$Replicate[which(Targets$FileName == j)]
    possibleDupIndex <- NULL
    for(y in sameFactorSamples){
      possibleDupIndex <- c(possibleDupIndex, Targets$Replicate[which(Targets$FileName == y)])
    }
    #DupNumSameFactorSamples <- Targets$Replicate[which(Targets$FileName == sameFactorSamples)]
    dupindex <- which(possibleDupIndex==currentDupNum)
    
    # 3 Replicates
    if(length(dupindex)==3){ ErrorOutput("Pipeline cannot deal with more than 2 replicates per condition. \nCheck your Target.txt file for details.") }
    if(length(dupindex)==2){      
      cat(" \n Replicate found:\n", sep="")
      for(l in sameFactorSamples[dupindex]){cat("\t\t Rep \t", l, "\n", sep="")}
    }
    
    if(Venn==TRUE){
      #--------------------------------------------------------
      # Load broadPeaks, narrowPeak
      #
      # broadPeak
      #
      path2currentSamplebroakPeak <- paste(path2MainFolder, "broadPeak/", j, ".broadPeak", sep="")
      if(file.exists(path2currentSamplebroakPeak)==FALSE){
        cat("\t", j, ".broadPeak file does not exist.\n", sep="")
      }
      if(file.exists(path2currentSamplebroakPeak)==TRUE){
        broadPeak <- read.table(path2currentSamplebroakPeak, header=FALSE)
        colnames(broadPeak) <- c("BedChr", "BedStart", "BedEnd", "-", "-", "-", "-", "BedStrand", "-")
        broadPeakGR <- Bed2GRanges(bedReply=broadPeak[, c("BedChr", "BedStart", "BedEnd", "BedStrand")])
        namebroadPeakGR <- paste(j, ".broadPeak.GRanges", sep="")
        assign(namebroadPeakGR, broadPeakGR)  
      }
      #
      # narrowPeak
      #
      path2currentSamplenarrowPeak <- paste(path2MainFolder, "narrowPeak/", j, ".narrowPeak", sep="")
      if(file.exists(path2currentSamplenarrowPeak)==FALSE){
        cat("\t", j, ".narrowPeak file does not exist.\n", sep="")
      }
      if(file.exists(path2currentSamplenarrowPeak)==TRUE){
        narrowPeak <- read.table(path2currentSamplenarrowPeak, header=FALSE)
        colnames(narrowPeak) <- c("BedChr", "BedStart", "BedEnd", "-", "-", "-", "-", "BedStrand", "-")
        narrowPeakGR <- Bed2GRanges(bedReply=narrowPeak[, c("BedChr", "BedStart", "BedEnd", "BedStrand")])
        namenarrowPeakGR <- paste(j, ".narrowPeak.GRanges", sep="")
        assign(namenarrowPeakGR, narrowPeakGR)  
      }
    }
    # end of if(Venn==TRUE)
    
    #--------------------------------------------------------
    # Replicates
    for(k in sameFactorSamples[dupindex]){

      nameCountTable <- paste("CountTable_", k, "_", nameOfBed, sep="")
      currentGR <- get(paste(k, ".bam.GRanges.RData", sep=""))
      countTable <- CountOverlaps2matrix(GRanges1=GRfeat, GRanges2=currentGR, binNumber=binNumber, geneNamesSingle=geneNames)
      head(countTable)
      unique(as.vector(seqnames(currentGR)))
      enrichments <- slice(currentGR)
      
      
      peaksHA <- import.bed(con="~/Documents/Sciences/Kingston/DIPG/DIPG_2014-07-07_ChIPseq/narrowPeak/PSa29-7_K27M1_Dox4d_HA_peaks.bed", asRangedData=F)
      cov <- coverage(peaksHA)      
      
      length(which(peaksHA$score>0))
      length(which(peaksHA$score>50))
      length(which(peaksHA$score>100))
      length(which(peaksHA$score>200))
      
      peaksHATop <- peaksHA[which(peaksHA$score>100)]
      countTable <- CountOverlaps2matrix(GRanges1=GRfeat, GRanges2=peaksHATop, binNumber=binNumber, geneNamesSingle=geneNames)
      
      mean(rowSums(countTable))
      TSSwHA <- countTable[which(rowSums(countTable)!=0),]
      head(TSSwHA)
      nrow(TSSwHA)
      nrow(countTable)
      nrow(TSSwHA)/nrow(countTable)
      
      # Assign countTable to the nameCountTable
      write.table(TSSwHA, paste(path2CountTables, "TSSwHAtop8pcPeaks_", Sys.Date(), ".bed", sep=""), sep = "\t", row.names = TRUE, col.names=TRUE, quote=FALSE, na="")
      
      # Create a Mart Object of these genes
      path2MartObject <- paste(path2Mart, nameOfBed, sep="")
      test <- LoadMartGRanges(path2MartObject=path2MartObject, binNumber=1, strand="-")
      geneNames <- GRfeat$GeneNames[seq(from=1, to=length(GRfeat), by=binNumber)]
      
      
      # Assign countTable to the nameCountTable
      assign(nameCountTable, countTable)
      # Save count table to folder CountTables
      cat(" \n Store countTable in : \t\t", paste(path2CountTables, sub(".bed", "", nameCountTable), "_bins", binNumber, "_", Sys.Date(), ".bed", sep=""), sep="")    
      cat(" \n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n", sep="") 
      write.table(get(nameCountTable), paste(path2CountTables, sub(".bed", "", nameCountTable), "_bins", binNumber, "_", Sys.Date(), ".bed", sep=""), sep = "\t", row.names = TRUE, col.names=TRUE, quote=FALSE, na="")
      # Store name of all active count tables in activeCTNames
      activeCTNames <- c(activeCTNames, nameCountTable)
    } # end of for(k in sameFactorSamples[dupindex])

    # Remove duplicated activeCTNames
    activeCTNames <- activeCTNames[!duplicated(activeCTNames)]
    activeCTNames
    
    #--------------------------------------------------------
    # Compute average countTables of 2 replicates
    #
    if(length(dupindex)==2){
      cat("\n\n +++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", sep="")
      cat(" \n Compute average of replicate countTables", sep="")    
      mylist <- list(get(activeCTNames[dupindex[1]]), get(activeCTNames[dupindex[2]]))
  #    str(mylist)
      dup1matrix <- simplify2array(mylist[[1]])
      dup2matrix <- simplify2array(mylist[[2]])
      # Compute the average between both duplicates
      averagedCountTable <- apply(simplify2array(mylist), 1:2, mean)
      # Plot the correlation per bin of each duplicate
      plotTitle <- paste("Correlation plot between replicate:\n", sameFactorSamples[dupindex[1]], "\n", sameFactorSamples[dupindex[2]], sep="")  
      plot(diag(cor(dup1matrix, dup2matrix)), xlab="bins", ylab="", type="b", ylim=c(-1,1), main=plotTitle, cex=0.5, cex.axis=0.8, cex.main=0.8)
      abline(h=0, lty=2)
      abline(h=mean(diag(cor(dup1matrix, dup2matrix))), lty=4, col="Red")
      
      # Assign countTable to the nameCountTable    
      nameCountTable <- paste("CountTable_AverageDup_", paste(sameFactorSamples[dupindex], collapse="--"), "_", nameOfBed, sep="")
      assign(nameCountTable, averagedCountTable)
      #activeCTNames <- c(activeCTNames, nameCountTable)
      activeCTNamesRep <- c(activeCTNamesRep, nameCountTable)
      # Save count table to folder CountTables        
      cat(" \n Store averaged countTable of replicate in : \t", nameCountTable, "_bins", binNumber, "_", Sys.Date(), ".bed", sep="")    
      write.table(get(nameCountTable), paste(path2CountTables, nameCountTable, sep=""), sep = "\t", row.names = TRUE, col.names=TRUE, quote=FALSE, na="")      
      
      cat("\n\n +++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n ", sep="")
    } # end of if
  } # end of j (go through samples that have the same factor)
  
  # If there is more than one one samples in sameFactorSamples (independently of whether there are replicates...)
  #     Consolidate all active countTable names
  if(length(sameFactorSamples)>1){
    currentTargets <- Targets[which(Targets$Factor==i),]
    df <- currentTargets$Replicate
    keepIndex <- which(currentTargets$Replicate == df[!(duplicated(df) | duplicated(df, fromLast = TRUE))])  
    activeCTNames <- c(activeCTNames[keepIndex], unique(activeCTNamesRep))
  }
  
  #--------------------------------------------------------
  # Prepare and plot metagenomic data
  #
  # Using countTable of the activeCTNames, find min and max range
  #
  rg <- range(colSums(get(activeCTNames[1])))
  miny <- rg[1]; maxy <- rg[2]
  for(k in 1:length(activeCTNames)){
    rg <- range(colSums(get(activeCTNames[k])))
    if(rg[1] < miny){miny <- rg[1]}
    if(rg[2] > maxy){maxy <- rg[2]}
  }

  # Clean activeCTNames to find FileShort of corresponding plot
  activeCTNames <- unique(activeCTNames)
  newActiveCTName <- sub("CountTable_", "", activeCTNames[1]); newActiveCTName <- sub("Norm2Inp_", "", newActiveCTName); newActiveCTName <- sub(paste("_", nameOfBed, sep=""), "", newActiveCTName); newActiveCTName <- sub(".bed", "", newActiveCTName)
  legends <- NULL
  # Smoothen data using runmean()
  plot(runmean(x=colSums(get(activeCTNames[1])), k=runmeank), type="l", col=cols[1], lwd=3, axes=F, xlab="Length [bp]", ylab="", ylim=c(miny, maxy), main=paste(i, nameOfBed," ; strand ", strand, sep= " "))
  #plot(colSums(get(activeCTNames[1])), type="l", col=cols[1], lwd=3, axes=F, xlab="Length [bp]", ylab="", ylim=c(miny, maxy), main=paste(i, nameOfBed," ; strand ", strand, sep= " "))
  for(h in 1:length(activeCTNames)){
    newActiveCTName <- sub("CountTable_", "", activeCTNames[h]); newActiveCTName <- sub("Norm2Inp_", "", newActiveCTName); newActiveCTName <- sub(paste("_", nameOfBed, sep=""), "", newActiveCTName); newActiveCTName <- sub(".bed", "", newActiveCTName); newActiveCTName <- sub("AverageDup_", "", newActiveCTName)
    newActiveCTName <- strsplit(newActiveCTName, "--")[[1]][1]
    legends <- c(legends,as.character(Targets$FileShort[which(Targets$FileName == newActiveCTName)]))
    head(get(activeCTNames[h]))
    points(runmean(x=colSums(get(activeCTNames[h])), k=runmeank), type="l", col=cols[h], lwd=3)  
    #points(colSums(get(activeCTNames[h])), type="l", col=cols[h], lwd=3)  
  }
  box()
  PlotAxes(mat=get(activeCTNames[1]), yrange=c(miny, maxy), strand=strand)  
  legend(1, maxy, legend=legends, lty=1, lwd=3, col=cols)  
  
  #--------------------------------------------------------
  # Prepare and plot Venn diagrams
  #
  # broadPeak venns  
  if(Venn==TRUE){
    
    broadPeakAll <- ls()[grep(".broadPeak.GRanges", ls())]
    if(length(broadPeakAll)!=0){
      path2VennsBroad <- paste(path2Plots, "broadPeakVenns/", sep="")
      cat("\n\n broadPeakAll: ", broadPeakAll, sep="\n")
      cat("\n Computing Venn diagrams of overlaps on broadPeak", sep="")
      mycols=c("Red", "Blue")
      for(v in 1:length(broadPeakAll)){
        for(w in 1:length(broadPeakAll)){
          if(v!=w){
            GRA <- get(broadPeakAll[v])
            GRB <- get(broadPeakAll[w])
            broadPeakNames <- c(broadPeakAll[v], broadPeakAll[w])
            AoverB <- sum(countOverlaps(GRA, GRB)!=0)
            venn.plot <- venn.diagram(list(1:length(GRA), (length(GRA)-AoverB):((length(GRA)-AoverB)+length(GRB))),  fill=mycols, alpha=rep(0.3, 2), cex = 2, cat.fontface=4, category.names=c("", ""), main=paste("Overlapping peaks", "\nRed = ", broadPeakNames[1], "\nBlue = ", broadPeakNames[2], sep=""), filename=paste(path2VennsBroad, broadPeakNames[1], "--", broadPeakNames[2], ".tiff", sep=""), cat.dist = rep(0.1, 2)) 
          } # end if
        } # end w
      } # end v  
      # Remove .broadPeak.GRanges objects
      rm(list=ls()[grep(".broadPeak.GRanges", ls())])
    } # end broadPeak venns  
    
    # narrowPeak venns  
    narrowPeakAll <- ls()[grep(".narrowPeak.GRanges", ls())]
    if(length(narrowPeakAll)!=0){
      path2VennsNarrow <- paste(path2Plots, "narrowPeakVenns/", sep="")
      cat("\n\n narrowPeakAll: ", narrowPeakAll, sep="\n")
      cat("\n Computing Venn diagrams of overlaps on narrowPeak", sep="")
      for(v in 1:length(narrowPeakAll)){
        for(w in 1:length(narrowPeakAll)){
          if(v!=w){
            GRA <- get(narrowPeakAll[v])
            GRB <- get(narrowPeakAll[w])
            narrowPeakNames <- c(narrowPeakAll[v], narrowPeakAll[w])
            AoverB <- sum(countOverlaps(GRA, GRB)!=0)  
            venn.plot <- venn.diagram(list(1:length(GRA), (length(GRA)-AoverB):((length(GRA)-AoverB)+length(GRB))),  fill=mycols, alpha=rep(0.3, 2), cex = 2, cat.fontface=4, category.names=c("", ""), main=paste("Overlapping peaks", "\nRed = ", broadPeakNames[1], "\nBlue = ", broadPeakNames[2], sep=""), filename=paste(path2VennsNarrow, narrowPeakNames[1], "--", narrowPeakNames[2], ".tiff", sep=""), cat.dist = rep(0.1, 2)) 
          } # end if
        } # end w
      } # end v  
      # Remove .broadPeak.GRanges objects
      rm(list=ls()[grep(".narrowPeak.GRanges", ls())])
    } # end narrowPeak venns  
  } # end Venn==TRUE
  
} # end of for(i in levels(Targets$Factor))
dev.off()

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
cat(" \n ||\t ChIPseq analysis version ", version," (© Patrick Schorderet 2014)", sep="")
cat(" \n ||", sep="")
cat(" \n ||\t Terminated with no known error.", sep="")
cat(" \n ||", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ======================================================================================\n\n", sep="")

# Close the R session when run from bash
quit(save = "no", status = 0, runLast = TRUE)

#}
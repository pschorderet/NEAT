#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  Bam2GRangesRemote.R                           *
#*                                                                *
#******************************************************************
#*                                                                *
#*----------------------------------------------------------------*
#*                                                                *
#* The below script creates GRanges objects from bam files        *
#*                                                                *
#*----------------------------------------------------------------*

cat("\n ----------------------------------------\n", sep="")
cat("\n Entered Bam2GRangesRemote.R \n", sep="")
cat("\n ----------------------------------------\n\n", sep="")

#*----------------------------------------------------------------*
# Initialize some parameters

args            <- commandArgs(TRUE)

path2expFolder  <- paste(args[1], "/", sep="")
path2bam	<- paste(args[2], "/", sep="")		# path to sample
path2bamRX	<- paste(args[2], "_RX/", sep="")
path2CustFct	<- paste(args[3], "/", sep="")		# path to CustomFunctions

path2MainFolder	<- path2expFolder
path2GRangesRData <- paste(path2MainFolder, "aligned/GRangesRData/", sep="")
path2GRangesRDataRX <- paste(path2MainFolder, "aligned/GRangesRData_RX/", sep="") 
path2statTable <- paste(path2MainFolder, "DataStructure/statTable.bed", sep="")

path2Targets	<- paste(path2MainFolder, "DataStructure/Targets.txt", sep="")
normConstant	<- 1e6


#-----------------------------------------------------------
# LOAD LIBRARY

library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(caTools)
library(rtracklayer)

source(paste(path2CustFct, "Bam2GRangesRData.R", sep=""))
source(paste(path2CustFct, "CountOverlaps2GRanges.R", sep=""))
source(paste(path2CustFct, "DoesTheFileExist.R", sep=""))
source(paste(path2CustFct, "ErrorOutput.R", sep=""))
source(paste(path2CustFct, "LoadMartGRanges.R", sep=""))
source(paste(path2CustFct, "Mart2GRanges.R", sep=""))
source(paste(path2CustFct, "Tags2GRanges.R", sep=""))
source(paste(path2CustFct, "OutputNumberOfReadsFromGRanges.R", sep=""))


#-----------------------------------------------------------
# Create folder
# wig
path2wig <- paste(path2MainFolder, "wig/", sep="")
checkExistenceOfFolder(path2CheckFolder=path2wig)

#------------------------------------------------------------
# Read Targets files

if(file.exists(path2Targets)==FALSE){ErrorOutput(paste("No Targets.txt file in\t", path2Targets, sep="")) }
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
  if(length(grep("Reference_genome_rx\t", newLine))==1) {
    currentline <- gsub("# ", "", newLine); currentline <- gsub("\t", "", currentline); currentline <- gsub("\"", "", currentline);
    refGenomeRX <- unlist(strsplit(currentline, split = "\\="))[2]
  }

}



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                                 #
#                Load reference genome parameters                 #
#                                                                 #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Read the chr_lens.dat
path2chrlens <- paste(path2MainFolder, "DataStructure/", refGenome, "/chr_lens.dat", sep="")
chromosomesFile <- read.table(path2chrlens, comment.char="#")
chromosomes <- chromosomesFile$V1

path2chrlensRX <- paste(path2MainFolder, "DataStructure/", refGenome, "/chr_lens.dat", sep="")
chromosomesFileRX <- read.table(path2chrlensRX, comment.char="#")
chromosomesRX <- chromosomesFileRX$V1


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
cat("\n Samples provided: ", allsamples, sep="\n")

#------------------------------------------------------------
# Check existence of .bam.GRanges.RData or .bam files
cat(" \n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Create / load GRanges object from bam files", sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")




#------------------------------------------------------------
# Create ChIPrx GRanges objects

# If .bam.GRanges.RData files DO NOT exists, create them
cat(" \n\n #########################################################", sep="")
cat("   \n #                                                       #", sep="")
cat("   \n #                Creating GRanges for RX                #", sep="")
cat("   \n #                                                       #", sep="")
cat("   \n ---------------------------------------------------------", sep="")

#------------------------------------------------------------
# Create GRanges objects
#
myChIPrx <- "TRUE"
countsRX <- NULL
bamGRangesRDataRXPath <- paste(paste(path2GRangesRDataRX, allsamples, sep=""), ".bam.GRanges.RData", sep="")
if(DoesTheFileExist(path2file=bamGRangesRDataRXPath)!=TRUE){
	# If .bam files cannot be found, assume RX is not used and set everything to 1
        bamPathRX <- paste(path2file=path2bamRX, allsamples, ".bam", sep="")
        if(DoesTheFileExist(path2file=bamPathRX)!=TRUE){
                # do nothing
		cat(" \n\n +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-", sep="")
                cat("   \n +-+-+-+-+-+-+-+-+-        NOTE       -+-+-+-+-+-+-+-+-+-+-\n", sep="")
		cat("   \n Some RX.bam files do not exists. \n Assume no ChIPrx normalization.\n", sep="")
		cat("   \n +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n\n\n", sep="")	
		myChIPrx <- "FALSE"		
        } else {
		# Create GRangesRData files
		cat(" \n\n **********************************************************\n", sep="")       
		cat(" \n Generating .bam.GRanges.RData files for RX", sep="")
		for(k in 1:length(bamPathRX)){
                        cat(" \n\n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
                        cat(" \n ", k, " / ", length(bamPathRX) ," :\t", bamPathRX[k], sep="")	
			Bam2GRangesRData(path2bam=bamPathRX[k], chromosomes=chromosomesRX)
		}
	}
}

#------------------------------------------------------------
# Load GRanges objects
#
if(myChIPrx == "TRUE"){

	GRangesSamplesRX <- list.files(path2GRangesRDataRX, pattern = "*.bam.GRanges.RData$")  
	if(DoesTheFileExist(path2file=paste(paste(path2GRangesRDataRX, allsamples, sep=""), ".bam.GRanges.RData", sep=""))==TRUE){
		# If .bam files cannot be found, assume RX is not used and set everything to 1
		bamPathRX <- paste(path2file=path2bamRX, allsamples, ".bam", sep="")
		if(DoesTheFileExist(path2file=bamPathRX)!=TRUE){
			countsRX <- c(countsRX, normConstant)
		}
		# Load each sample into its .GRanges named object
		for(i in 1:length(GRangesSamplesRX)){
			path2load <- paste(path2GRangesRDataRX, GRangesSamplesRX[i], sep="")
			load(path2load)
			assign(GRangesSamplesRX[i], GRangesObj)
			countsRX <- c(countsRX, length(get(GRangesSamplesRX[i])))
#			cat("\n", GRangesSamplesRX[i], "\t", length(get(GRangesSamplesRX[i])), sep="")
		}
	}

	OutputNumberOfReadsFromGRanges(GRangesSamplesRX)

}


cat(" \n\n #########################################################", sep="")
cat("   \n #                                                       #", sep="")
cat("   \n #                    Creating GRanges                   #", sep="")
cat("   \n #                                                       #", sep="")
cat("   \n ---------------------------------------------------------", sep="")

#------------------------------------------------------------
# Load GRanges objects
#
bamGRangesRDataPath <- paste(paste(path2GRangesRData, allsamples, sep=""), ".bam.GRanges.RData", sep="")

if(DoesTheFileExist(path2file=bamGRangesRDataPath)!=TRUE){
	# If .bam files cannot be found, abort
	bamPath <- paste(path2file=path2bam, allsamples, ".bam", sep="")
	if(DoesTheFileExist(path2file=bamPath)!=TRUE){
		ErrorOutput("Some .bam files could not be found")
	} else {
		# Create GRangesRData files
		cat(" \n\n **********************************************************\n", sep="")
		cat(" \n Generating .bam.GRanges.RData files", sep="")
		for(k in 1:length(bamPath)){    
			cat(" \n\n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
			cat(" \n ", k, " / ", length(bamPath) ," :\t", bamPath[k], sep="")
			Bam2GRangesRData(path2bam=bamPath[k], chromosomes=chromosomes)
		}  
	}
}

#------------------------------------------------------------
# Load GRanges objects
#
counts <- NULL
GRangesSamples <- list.files(path2GRangesRData, pattern = "*.bam.GRanges.RData$")
if(DoesTheFileExist(path2file=paste(paste(path2GRangesRData, allsamples, sep=""), ".bam.GRanges.RData", sep=""))==TRUE){
	# Load each sample into its .GRanges named object
	for(i in 1:length(GRangesSamples)){
		path2load <- paste(path2GRangesRData, GRangesSamples[i], sep="")
		load(path2load)
		assign(GRangesSamples[i], GRangesObj)
		counts <- c(counts, length(get(GRangesSamples[i])))
#		cat("\n", GRangesSamples[i], "\t", length(get(GRangesSamples[i])), sep="")
	}
}
OutputNumberOfReadsFromGRanges(GRangesSamples)



cat(" \n\n #########################################################", sep="")
cat("   \n #                                                       #", sep="")
cat("   \n #                 Creating Summary table                #", sep="")
cat("   \n #                                                       #", sep="")
cat("   \n ---------------------------------------------------------\n", sep="")

#------------------------------------------------------------
# Create summary table with number of reads and scaling factor
#
cat("\n statTable: \n\n", sep="")

if(myChIPrx == "TRUE"){
	# Compute normalization factor. This corresponds to alpha in the Orlando et al. Cell Reports 2014.
	# Multiplyin
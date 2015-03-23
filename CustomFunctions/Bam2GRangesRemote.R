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
path2GRanges	<- paste(args[3], "/", sep="") 
path2GRangesRX	<- paste(args[3], "_RX/", sep="")
path2CustFct	<- paste(args[4], "/", sep="")		# path to CustomFunctions


path2MainFolder	<- path2expFolder
path2Targets	<- paste(path2MainFolder, "DataStructure/Targets.txt", sep="")
path2statTable 	<- paste(path2MainFolder, "DataStructure/statTable.bed", sep="")
path2wig 	<- paste(path2MainFolder, "/wig", sep="")
normConstant	<- 1e6


#-----------------------------------------------------------
# LOAD LIBRARY

library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(caTools)
library(rtracklayer)

source(paste(path2CustFct, "Bam2GRangesRData.R", sep=""))
source(paste(path2CustFct, "CheckExistenceOfFolder.R", sep=""))
source(paste(path2CustFct, "CountOverlaps2GRanges.R", sep=""))
source(paste(path2CustFct, "DoesTheFileExist.R", sep=""))
source(paste(path2CustFct, "ErrorOutput.R", sep=""))
source(paste(path2CustFct, "LoadMartGRanges.R", sep=""))
source(paste(path2CustFct, "Mart2GRanges.R", sep=""))
source(paste(path2CustFct, "Tags2GRanges.R", sep=""))
source(paste(path2CustFct, "OutputNumberOfReadsFromGRanges.R", sep=""))


#------------------------------------------------------------
# Read Targets files

if(file.exists(path2Targets)==FALSE) { ErrorOutput(paste("No Targets.txt file in\t", path2Targets, sep="")) }
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
	if(length(grep("Proj_TaxonBSgenome\t", newLine))==1) {
                currentline <- gsub("# ", "", newLine); currentline <- gsub("\t", "", currentline); currentline <- gsub("\"", "", currentline);
                BSgenome <- unlist(strsplit(currentline, split = "\\="))[2]
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

#cat("\n\n", allsamples, sep="\n")

bamGRangesRXPath <- paste(paste(path2GRangesRX, allsamples, sep=""), ".bam.GRanges.RData", sep="")
if(DoesTheFileExist(path2file=bamGRangesRXPath)!=TRUE){
	# If .bam files cannot be found, assume RX is not used and set everything to 1
        bamPathRX <- paste(path2file=path2bamRX, allsamples, ".bam", sep="")
        
	if(DoesTheFileExist(path2file=bamPathRX)!=TRUE){
                # do nothing
		cat(" \n\n +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-", sep="")
                cat("   \n +-+-+-+-+-+-+-+-+-        NOTE       -+-+-+-+-+-+-+-+-+-+-\n", sep="")
		cat("   \n Some RX.bam files do not exists. \n Assume no ChIPrx normalization.\n", sep="")
		cat("   \n +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n\n\n", sep="")	
		myChIPrx <- "FALSE"		
		countsRX <- rep(normConstant, length(allsamples))
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

	GRangesSamplesRX <- list.files(path2GRangesRX, pattern = "*.bam.GRanges.RData$")  
	if(DoesTheFileExist(path2file=paste(paste(path2GRangesRX, allsamples, sep=""), ".bam.GRanges.RData", sep=""))==TRUE){
	
		# If .bam files cannot be found, assume RX is not used and set everything to 1
		bamPathRX <- paste(path2file=path2bamRX, allsamples, ".bam", sep="")
		if(DoesTheFileExist(path2file=bamPathRX)!=TRUE){
			countsRX <- c(countsRX, normConstant)
		}
		# Load each sample into its .GRanges named object
		for(i in 1:length(GRangesSamplesRX)){
			path2load <- paste(path2GRangesRX, GRangesSamplesRX[i], sep="")
			load(path2load)
			assign(GRangesSamplesRX[i], GRangesObj)
			countsRX <- c(countsRX, length(get(GRangesSamplesRX[i])))
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
bamGRangesPath <- paste(paste(path2GRanges, allsamples, sep=""), ".bam.GRanges.RData", sep="")

if(DoesTheFileExist(path2file=bamGRangesPath)!=TRUE){
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
GRangesSamples <- list.files(path2GRanges, pattern = "*.bam.GRanges.RData$")
if(DoesTheFileExist(path2file=paste(paste(path2GRanges, allsamples, sep=""), ".bam.GRanges.RData", sep=""))==TRUE){
	# Load each sample into its .GRanges named object
	for(i in 1:length(GRangesSamples)){
		path2load <- paste(path2GRanges, GRangesSamples[i], sep="")
		load(path2load)
		assign(GRangesSamples[i], GRangesObj)
		counts <- c(counts, length(get(GRangesSamples[i])))
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
NormFactRX <- round(normConstant/countsRX, 2)


statTable <- cbind(allsamples, counts, countsRX, NormFactRX)
colnames(statTable) <- c("FileName", "Counts", "CountsRX", "NormFactRX")

write.table(statTable, file=path2statTable, quote=F, sep="\t", row.names=F, col.names=T)

cat("\n statTable: \n\n", sep="")
print(statTable)



cat(" \n\n #########################################################", sep="")
cat("   \n #                                                       #", sep="")
cat("   \n #                   Creating wig files                  #", sep="")
cat("   \n #                                                       #", sep="")
cat("   \n ---------------------------------------------------------\n", sep="")

CheckExistenceOfFolder(path2wig)

cat("\n\n", sep="")

# Set variables and paths
seqlengths <- c(chromosomesFile[,2])
names(seqlengths) <- chromosomesFile[,1]
#seqlengths

# Slice the genome into bins
bins <- tileGenome(seqlengths, tilewidth=500, cut.last.tile.in.chrom=TRUE)

# Read in the statTable file countaining the ChIPrx normaization factors
path2statTable <- paste(path2MainFolder, "DataStructure/statTable.bed", sep="")
statTable <- read.table(path2statTable, header=TRUE, comment.char="#")

# Disable scientific notation to ensure no 'e' are found in the wig files
options(scipen=999)
# Set paths to tmp files

path2wigfilename_tmp <- paste(path2wig, "/tmp.wig", sep="")
path2wigfile_header <- paste(path2wig, "/header.wig", sep="")

for(h in 1:length(GRangesSamples)){

	cat("\n - - - - - - - - - - - - - - - - - - - - - - - - - - - - -", sep="")
	cat("\n Creating wig file:\t", GRangesSamples[h], "\n", sep="")

	# h=1
	gr <- CountOverlaps2GRanges(GRanges1=bins, GRanges2=get(GRangesSamples[h]), normFactRX=statTable$NormFactRX[h], normToLibrarySize=FALSE)
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

cat("\n\n", sep="")

# Delete temporary files
mycode <- paste("`rm ", path2wigfilename_tmp, "`", sep="")
system(mycode)
mycode <- paste("`rm ", path2wigfile_header, "`", sep="")
system(mycode)



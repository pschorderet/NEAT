#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*			QC for fastq                              *
#*                                                                *
#******************************************************************

cat("\n ----------------------------------------\n", sep="")
cat("\n Entered QC.R \n", sep="")
cat("\n ----------------------------------------\n\n", sep="")

#*----------------------------------------------------------------*
# Initialize some parameters
args		<- commandArgs(TRUE)
code		<- args[0] # path to peakcalling.R
path2expFolder	<- args[1] # path to the main experiment folder

#-----------------------------------------------------------
# LOAD LIBRARY
cat("\n Load systemPipeR library \n", sep="")
library(systemPipeR)

#*----------------------------------------------------------------*
# Define paths
path2Targets	<- paste( path2expFolder , "/DataStructure/Targets.txt", sep="")
path2fastq	<- paste( path2expFolder , "/fastq", sep="")
path2QC		<- paste( path2expFolder, "/QC/", sep="")

setwd(path2QC)

#--------------------------------------------------
# Read in Targets.txt
Targets <- read.delim(path2Targets, comment.char="#")

#--------------------------------------------------
# Create virtuel paths to all fastq
Targets			<- read.delim(path2Targets, comment.char="#")
samples			<- Targets$FileName
path2samplesfastq	<- paste(path2fastq, "/", samples, ".fastq", sep="")
names(path2samplesfastq) <- Targets$FileName

#--------------------------------------------------
# Load data in seeFastq object
cat("\n\n Load all fastq files in fqlist object \n\n", sep="")
fqlist			<- seeFastq(fastq=path2samplesfastq, batchsize=100000, klength=8)

#--------------------------------------------------
# Plot QC report
cat("\n\n Save QC report in \t", path2QC, "\n\n", sep="")
pdf(paste(path2QC, "/QCReport.pdf", sep=""), height=18, width=4*length(fqlist))
seeFastqPlot(fqlist)
dev.off()



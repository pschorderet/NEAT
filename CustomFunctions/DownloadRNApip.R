

DownloadRNApip <- function(LocalPath2NEAT, LocalPath2NewProject, RemotePath2MainFolderName){

#******************************************************************
  
version <- "1.0.1 Jan 2015"

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*                Load and set parameters                         *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

cat(" \n\n======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  ", sep="")
cat(" \n||", sep="")
cat(" \n||\t DownloadRNApip, version ", version," (© Patrick Schorderet 2014)", sep="")
cat(" \n||", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n======================================================================================\n\n\n", sep="")

#
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
#                                                                 #
# This script is intended to be run from bash/shell Mac terminal  #
#                                                                 #
# The script:                                                     #
#     - Creates the proper arborescence for RNAseq.R              #
#     - Downlaods filtered .bam files from remote Canute server   #
#     - Stores files in proper location (see below for arborescence)#
#                                                                 #
#*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*
#                                                                 #
# Experiment-specific parameters                                  #
#                                                                 #
# This needs to be the same name as the MainFolder in Canute!     #
# LocalPath2NEAT <- "~/NEAT/"
# Localpath2NewProject <- "~/Desktop/"
# RemotePath2MainFolderName <- "/data/schorderet/projects/RNA/"
# sshpath <- "schorderet@canute.mgh.harvard.edu"
#*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*                Load and set parameters                         *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

cat(paste("\n RemotePath2MainFolderName \t ", RemotePath2MainFolderName, "\n",sep=""))

NewMainFolderName <- tail(unlist(strsplit(RemotePath2MainFolderName, split="/")), 1)
if(NewMainFolderName==""){
  NewMainFolderName <- tail(unlist(strsplit(RemotePath2MainFolderName, split="/")), 2)[1]
}

#------------------------------------------------------------
# Read main folder

cat(paste(" LocalPath2NewProject \t\t ", LocalPath2NewProject, "\n",sep=""))
if(file.exists(LocalPath2NewProject)==TRUE){cat(" \n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n\n ", NewMainFolderName, " exists. \n Reading the Targets.txt file \n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n\n\n ", sep="");}
if(file.exists(LocalPath2NewProject)==FALSE){cat(" \n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n\n Creating\t", NewMainFolderName, "\n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n ", sep="");dir.create(LocalPath2NewProject)}


#--------------------------------------------------
# Define Paths using user's parameters
LocalPath2CustomFunctions <- paste( LocalPath2NEAT, "CustomFunctions", "/", sep="")
LocalPath2Targets <- paste( LocalPath2NewProject, "DataStructure/Targets.txt", sep="")
LocalPath2bam <- paste( LocalPath2NewProject, "bam/", sep="")
LocalPath2Tophat <- paste( LocalPath2NewProject, "Tophat/", sep="")

source(paste( LocalPath2CustomFunctions, "ErrorOutput.R", sep=""))

#------------------------------------------------------------
# Redirect output to log file
LocalPath2Logs <- paste(LocalPath2NewProject, "logs/", sep="") 
if(file.exists(LocalPath2Logs)==TRUE){cat(" \n Logs file exists.\n\n", sep="")}
if(file.exists(LocalPath2Logs)==FALSE){cat(" \n Creating\t", LocalPath2Logs, "\n\n", sep="");dir.create(LocalPath2Logs)}
sink(paste(LocalPath2Logs, Sys.Date(), "_", format(Sys.time(), "%X"), ".txt", sep=""), append=FALSE, split=FALSE)

#--------------------------------------------------
# Create proper arborescence
cat(" \n======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n||\t Create / check arborescence of ", LocalPath2NewProject, sep="")
cat(" \n|| .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")
source(paste(LocalPath2CustomFunctions, "RNAseqCreateArborescence.R", sep=""))
RNAseqCreateArborescence(path2MainFolder=LocalPath2NewProject)

#--------------------------------------------------
# Read Targets.txt file and find ssh parameters
res <- readLines(LocalPath2Targets)
for(i in 1:length(res)){
  newLine <- res[i]
  # Store some variables
  if(length(grep("My_personal_ssh", newLine))==1) { 
    currentline <- gsub("# ", "", newLine); currentline <- gsub("\t", "", currentline); currentline <- gsub("\"", "", currentline);
    sshpath <- unlist(strsplit(currentline, split = "\\="))[2]
  }    
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*             Transfer data from remote server                   *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#--------------------------------------------------
# Transfer file DataStructure from Canute to Datastructure/ folder
#
mycode <- paste("`scp -r ", sshpath, ":", RemotePath2MainFolderName, "DataStructure/ " , LocalPath2NewProject, "`", sep="")
cat(" \n\n\n======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n||\t Opening SSH connection (please enter your password in the terminal window...)", sep="")
cat(" \n|| .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n\n", sep="")
cat(" \n Downloading Targtes.txt file from remote server", sep="")
system(mycode)
cat(mycode)

#--------------------------------------------------
# Transfer QC bigwig, narrowPeak and broadPeak folders from *RemotePath2MainFolderName*
#
cat(" \n\n\n======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n||\t Downloading QC folder", sep="")
cat(" \n|| .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n\n", sep="")

# Download QC folder
cat(" \n Downloading QC folder from remote server", sep="")
mycode <- paste("`scp -r ", sshpath, ":", RemotePath2MainFolderName, "/QC/ " , LocalPath2NewProject, "`", sep="")
cat(mycode)
system(mycode)


#--------------------------------------------------
# Transfer Tophat folder (including filtered .bam files) found in *RemotePath2MainFolderName*
#
cat(" \n\n\n======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n||\t Downloading and opening Targets.txt files from remote server", sep="")
cat(" \n|| .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n\n", sep="")

if(file.exists(LocalPath2Targets)==FALSE){ErrorOutput(paste("No Targets.txt file in\t", LocalPath2Targets, sep="")) }
cat(" \n\n Targets file provided: \n\n", sep="")
Targets <- read.delim(LocalPath2Targets, comment.char="#")
print(Targets)

# Download Tophat folder
cat(" \n Downloading Tophat folder from remote server", sep="")
mycode <- paste("`scp -r ", sshpath, ":", RemotePath2MainFolderName, "/Tophat " , LocalPath2NewProject, "`", sep="")
system(mycode)

# Load datasets provided in Targets file
samples <- Targets$FileName
allsamples <- c(levels(Targets$FileName), levels(Targets$InpName))
# Remove '-' when no input is provided
if(length(which(allsamples=="-")) != 0){
  allsamples <- allsamples[-which(allsamples=="-")]    
}


#--------------------------------------------------
# At local level
#
cat(" \n\n\n ======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n||\t Downloading and reformating files and folders at local level", sep="")
cat(" \n|| .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n\n", sep="")
# Move .bam files before deleting entire folder
for(i in 1:length(allsamples)){
  # Download sample i to .bam folder
  # i=1
  cat("\n Sample: \t",  allsamples[i], "\n", sep="")
  # Move sample i to .bam folder
  cat("\n\t Moving \t\t ", allsamples[i], ".bam\t to\t" , LocalPath2bam, sep="")
  movecode <- paste("`mv ", LocalPath2Tophat, "/", allsamples[i], "/", allsamples[i], ".bam ", LocalPath2bam, "`", sep="")
  movecode
  system(movecode)
  cat("\n\n .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")
}


cat(" \n\n\n ======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n||\t Cleaning up local files and folders", sep="")
cat(" \n|| .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n\n", sep="")
cat("\n\n\t ", length(allsamples), " files have been added to \t\t", LocalPath2bam, sep="")
cat("\n\n .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                                 #
#                  Session info and closure                       #
#                                                                 #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Print sessionInfo()
cat(" \n\n\n|| -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-   ", sep="")
cat(" \n\n", sep="")
print(sessionInfo())
cat(" \n", sep="")
cat(" \n|| -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n\n   ", sep="")

cat(" \n\n\n======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  ", sep="")
cat(" \n||", sep="")
cat(" \n||\t Download RNApip, version ", version," (© Patrick Schorderet 2014)", sep="")
cat(" \n||", sep="")
cat(" \n||\t Terminated with no known error.", sep="")
cat(" \n||", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n======================================================================================\n\n", sep="")

# Close the R session when run from bash
quit(save = "no", status = 0, runLast = TRUE)

}


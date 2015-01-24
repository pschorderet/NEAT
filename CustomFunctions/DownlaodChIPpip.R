
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
cat(" \n||\t DownloadChIPpip, version ", version," (© Patrick Schorderet 2014)", sep="")
cat(" \n||", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n======================================================================================\n", sep="")


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*      Steps before running ChIPseq1.R pipeline for local        *
#*                                                                *
#******************************************************************
#
#
#
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
#                                                                 #
# This script is intended to be run from bash/shell Mac terminal  #
#                                                                 #
# The script:                                                     #
#     - Creates the proper arborescence for ChIPseq.R             #
#     - Downlaods required files from remote server               #
#     - Stores files in proper location (see arborescnce below)   #
#                                                                 #
#                                                                 #
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
#                                                                 #
# Setup ssh keys: Follow protocol here:                           #
#     http://www.howtogeek.com/66776/how-to-remotely-copy-files-over-ssh-without-entering-your-password/
#   as well as :                                                  #
#     https://github.com/markcarver/mac-ssh-askpass               #
#                                                                 #
#                                                                 #
#*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*
#                                                                 #
# Experiment-specific parameters                                  #
#                                                                 #
# This needs to be the same name as the MainFolder in Canute!     #
# NewMainFolderName <- "EXAMPLE/"
# LocalPath2NEAT <- "~/PepsRscripts/RScripts/"
# sshpath <- "schorderet@canute.mgh.harvard.edu"
# RemotePath2MainFolderName <- "/data/schorderet/projects/EXAMPLE/"
#*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arborescence of files local:                                    #
#                                                                 #
# *...* means the name is variable                                #
#                                                                 #
# ~/PepsRscripts/RScripts/                                        #
#                                                                 #
#     └─ *MainFolder*                                             #
#           └─ bam                                                #
#                 └─ *samples1*.bam                               #
#                 └─ ...                                          #
#                                                                 #
#           └─ DataStructure                                      #
#                 └─ Targets.txt                                  #
#                                                                 #
#           └─ GrangesRData                                       #
#                 └─ *samples1*.bam.GRangesRData                  #
#                 └─ ....                                         #
#                                                                 #
#           └─ CountTables                                        #
#                 └─ RPKM.bed                                     #
#                 └─ ....                                         #
#                                                                 #
#           └─ Plots                                              #
#                 └─ RPKM.pdf                                     #
#                 └─ ....                                         #
#                                                                 #
#     └─ MartObject                                               # 
#           └─ mm9_Trasncripts.bed                                #
#           └─ ...                                                #
#                                                                 #
#     └─ ReferenceFiles                                           #
#           └─ mm9KG.sqlite                                       #
#           └─ Mus_musculus.NCBIM37.67.allchr.gtf                 #
#           └─ UC_ID_to_GENE_symbol_mm9.txt                       #
#           └─ ...                                                #
#                                                                 #
#     └─ PepsFunctions                                            #
#           └─ ...                                                #
#                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                 #
# Arborescence of files remote                                    #
#                                                                 #
# *...* means the name is variable                                #
#                                                                 #
# /data/*path2ChIPseqScripts*   	                                #
#                                                                 #
# /data/*userFolder*                                              #
#                                                                 #
#                                                                 #
#     └─ * expFolder *	                                          #
#                                                                 #
#           └─ bwa_saf/                                           #
#	              └─ *samples1*     	                              #
#						        └─ *samples1*.bai     	                      #
#						        └─ *samples1*.bam     	                      #
#						        └─ splitbycrh	     	                          #
#							          └─ ...	 	              		              #
#  						      └─ ...    	 	              		              #
#                                                                 #
#           └─ bwa_sai/                                       	  #
#	              └─ *samples1*.sai   	                            #
#               └─ ...                                            #
#                                                                 #
#           └─ bwa_sam                                            #
#	              └─ *samples1*.sam   	                            #
#	              └─ *samples1*.u.sam   	                          #
#	              └─ *samples1*.u.sortedwpcr.bam     	              #
#	              └─ *samples1*.u.unsorted.bam     	                #
#               └─ ...                                            #
#                                                                 #
#           └─ DataStructure                                      #
#               └─ Targets.txt                                    #
#               └─ AdvancedSettings.txt                           #
#                                                                 #
#           └─ fastq                                              #
#	              └─ *samples1*.fastq                               #
#                                                                 #
#           └─ peakcalling                                        #
#               └─ bigwig/                                        #
#                   └─ *samples1*.bw                              #
#                   └─ ...                                        #
#                                                                 #
#               └─ broadPeak/                                     #
#	                  └─ *samples1*.broadPeak                       #
#	                  └─ ...             	                          #
#                                                                 #
#               └─ narrowPeak                                     #
#                   └─ *samples1*.narrowpeak                      #
#                   └─ ...                                        #
#                                                                 #
#               └─ *samples1*_fdr_*fdr*   	                      #
#                   └─ *samples1*.bindingpositions.txt            #
#                   └─ *samples1*.crosscorrelation.pdf            #
#                   └─ *samples1*.density.wig                     #
#                   └─ ...                                        #
#                                                                 #
#           └─ QC                                                 #
#               └─ QCreport.pdf                                   #
#                                                                 #
#           └─ scripts                                            #
#	              └─ *samples1*.bwa.sh	                            #
#               └─ *samples1*.filter.sh	                          #
#               └─ *samples1*.map.sh                              #
#               └─ *samples1*.unzip.sh                            #
#               └─ ...                                            #                                                                 #
#                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


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
#cat(paste("\n NewMainFolderName \t\t ", NewMainFolderName, "\n",sep=""))

#------------------------------------------------------------
# Create mainFolder 

LocalPath2NewProject <- paste( Localpath2NewProject, NewMainFolderName, "/", sep="")
#LocalPath2NewProject <- paste( LocalPath2NEAT, NewMainFolderName, "/", sep="")
cat(paste(" LocalPath2NewProject \t\t ", LocalPath2NewProject, "\n",sep=""))
if(file.exists(LocalPath2NewProject)==TRUE){cat(" \n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n\n ", NewMainFolderName, " already exists. \n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n\n\n ", sep=""); stop()}
if(file.exists(LocalPath2NewProject)==FALSE){cat(" \n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n\n Creating\t", NewMainFolderName, "\n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n ", sep="");dir.create(LocalPath2NewProject)}

#--------------------------------------------------
# Define Paths using user's parameters
LocalPath2CustomFunctions <- paste( LocalPath2NEAT, "CustomFunctions", "/", sep="")
LocalPath2Targets <- paste( LocalPath2NewProject, "DataStructure/Targets.txt", sep="")
LocalPath2bam <- paste( LocalPath2NewProject, "bam/", sep="")
LocalPath2saf <- paste( LocalPath2bam, "bwa_saf/", sep="" )

source(paste( LocalPath2CustomFunctions, "/ErrorOutput.R", sep=""))

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
source(paste(LocalPath2CustomFunctions, "ChIPseqCreateArborescence.R", sep=""))
ChIPseqCreateArborescence(path2MainFolder=LocalPath2NewProject)


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
# quikie-mart7

#--------------------------------------------------
# Transfer QC bigwig, narrowPeak and broadPeak folders from *RemotePath2MainFolderName*
#
cat(" \n\n\n======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n||\t Downloading QC, bigwig, narrowPeaks and broadPeaks folders", sep="")
cat(" \n|| .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n\n", sep="")

Sys.sleep(3)
# Download QC folder
cat(" \n Downloading QC folder from remote server", sep="")
mycode <- paste("`scp -r ", sshpath, ":", RemotePath2MainFolderName, "/QC/ " , LocalPath2NewProject, "`", sep="")
cat(mycode)
system(mycode)
# Download bigwig folder
cat(" \n Downloading QCReport.pdf from remote server", sep="")
mycode <- paste("`scp -r ", sshpath, ":", RemotePath2MainFolderName, "/peakcalling/bigwig " , LocalPath2NewProject, "`", sep="")
system(mycode)
# Download broadPeak folder
cat(" \n Downloading QCReport.pdf from remote server", sep="")
mycode <- paste("`scp -r ", sshpath, ":", RemotePath2MainFolderName, "/peakcalling/broadPeak " , LocalPath2NewProject, "`", sep="")
system(mycode)
# Download narrowPeak folder
cat(" \n Downloading QCReport.pdf from remote server", sep="")
mycode <- paste("`scp -r ", sshpath, ":", RemotePath2MainFolderName, "/peakcalling/narrowPeak " , LocalPath2NewProject, "`", sep="")
system(mycode)

#--------------------------------------------------
# Transfer filtered .bam files found in *RemotePath2MainFolderName*
#
cat(" \n\n\n======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n||\t Downloading and opening Targets.txt files from remote server", sep="")
cat(" \n|| .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n\n", sep="")

if(file.exists(LocalPath2Targets)==FALSE){ErrorOutput(paste("No Targets.txt file in\t", LocalPath2Targets, sep="")) }
cat(" \n\n Targets file provided: \n\n", sep="")
Targets <- read.delim(LocalPath2Targets, comment.char="#")
print(Targets)
cat(" \n Downloading .bam files from remote server", sep="")
mycode <- paste("`scp -r ", sshpath, ":", RemotePath2MainFolderName, "bwa_saf/ " , LocalPath2bam, "`", sep="")
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
  cat("\n\t Downloading \t", allsamples[i], "\t folder from remote server to local", sep="")
  # Move sample i to .bam folder
  cat("\n\t Moving \t\t ", allsamples[i], ".bam\t to\t" , LocalPath2bam, sep="")
  movecode <- paste("`mv ", LocalPath2saf, allsamples[i], "/", allsamples[i], ".bam ", LocalPath2bam, "`", sep="")
  movecode
  system(movecode)
  # Delete rest of folder
  #cat("\n\t Deleting \t\t ", LocalPath2saf, allsamples[i], sep="")
  #delcode <- paste("`rm -r ", LocalPath2saf, allsamples[i], "/`", sep="")
  #system(delcode)
  cat("\n\n .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")
}
# Delete rest of folder
cat("\n\t Deleting \t\t ", LocalPath2saf, sep="")
delcode <- paste("`rm -r ", LocalPath2saf, "`", sep="")
system(delcode)



cat(" \n\n\n ======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n||\t Cleaning up local files and folders", sep="")
cat(" \n|| .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n\n", sep="")
cat("\n\n\t ", length(allsamples), " files have been added to \t\t", LocalPath2bam, sep="")
cat("\n\t ", length(allsamples), " folders have been deleted from \t", LocalPath2bam, sep="")
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
cat(" \n||\t Remote2Local for ChIPseq, version ", version," (© Patrick Schorderet 2014)", sep="")
cat(" \n||", sep="")
cat(" \n||\t Terminated with no known error.", sep="")
cat(" \n||", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n======================================================================================\n\n", sep="")

# Close the R session when run from bash
quit(save = "no", status = 0, runLast = TRUE)

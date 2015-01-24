
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*      Steps before running RNAseq1.R pipeline for local         *
#*                                                                *
#******************************************************************
#
#
version <- "0.1.3 Oct 2014"
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
#   Connect to the ssh                                            #
#       ssh schorderet@canute.mgh.harvard.edu
#       quikie-mart7
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
#NewMainFolderName <- "DIPG_2014-07-07_RNAseq/";                  #
#path2RScripts <- "~/PepsRscripts/RScripts/"                      #
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
# Arborescence of files remote (Canute):                          #
#                                                                 #
# *...* means the name is variable                                #
#                                                                 #
# /data/*path2ChIPseqScripts*   	                                #
#                                                                 #
# /data/*userFolder*                                              #
#                                                                 #
#                                                                 #
#     └─ * expFolder *	                                          #
#           └─ DataStructure                                      #
#                 └─ Targets.txt                                  #
#                                                                 #
#           └─ fastq/                                             #
#                └─ *samples1*.fastq                              #
#                                                                 #
#           └─ HTseq/                                             #
#                                                                 #
#           └─ logs                                       	      #
#	              └─ *samples1*.bwa.logs                            #
#	              └─ *samples1*.filter.logs                         #
#                   └─ ....                                       #
#                                                                 #
#           └─ saf/                                               #
#                                                                 #
#           └─ scripts                                            #
#                └─ *samples1*.bwa.sh	                            #
#                                                                 #
#           └─ Tophat                                             #
#                └─ *samples1*                                    #
#                   └─ accepted_hits.bam                          #
#                   └─ deletions.bed                              #
#                   └─ insertions.bed                             #
#                   └─ junctions.bed                              #
#                   └─ logs                                       #
#                       └─ ....                                   #
#                   └─ prep_reads.info                            #
#                   └─ tmp                                        #
#                       └─ ....                                   #
#                   └─ unmapped.bam                               #
#                └─ ...                                           #
#                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*                Load and set parameters                         *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#--------------------------------------------------
# Define Paths
# - - - - - - - - - - - - - - - - - - - - - - - - #
path2expFolderCanute <- paste("/data/schorderet/projects/", NewMainFolderName, "/", sep="")
# - - - - - - - - - - - - - - - - - - - - - - - - #
path2MainFolder <- paste(path2RScripts, NewMainFolderName, sep="")
path2PepsRScripts <- paste(path2RScripts,"PepsFunctions/", sep="")
path2Targets <- paste(path2MainFolder, "DataStructure/Targets.txt", sep="")
path2bam <- paste(path2MainFolder, "bam/", sep="")
path2Logs <- paste(path2MainFolder, "Logs/", sep="")
path2PepsRScripts <- paste(path2RScripts,"PepsFunctions/", sep="")

source(paste(path2PepsRScripts, "ErrorOutput.R", sep=""))

#------------------------------------------------------------
# Create mainFolder 
if(file.exists(path2MainFolder)==TRUE){cat(" \n\n\n\n\n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n\n ", NewMainFolderName, " already exists. \n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n\n\n\n\n\n ", sep=""); stop()}
if(file.exists(path2MainFolder)==FALSE){cat(" \n\n\n\n\n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n\n Creating\t", path2MainFolder, "\n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n ", sep="");dir.create(path2MainFolder)}

#------------------------------------------------------------
# Redirect output to log file
if(file.exists(path2Logs)==TRUE){cat(" \n Logs file exists.\n\n", sep="")}
if(file.exists(path2Logs)==FALSE){cat(" \n Creating\t", path2Logs, "\n\n", sep="");dir.create(path2Logs)}
sink(paste(path2Logs, Sys.Date(), "_", format(Sys.time(), "%X"), ".txt", sep=""), append=FALSE, split=FALSE)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*                Load and set parameters                         *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

cat(" \n\n======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  ", sep="")
cat(" \n||", sep="")
cat(" \n||\t Canute2Local for RNAseq, version ", version," (© Patrick Schorderet 2014)", sep="")
cat(" \n||", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n======================================================================================\n\n\n", sep="")

#--------------------------------------------------
# Create proper arborescence
cat(" \n\n======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n||\t Create / check arborescence of ", path2MainFolder, sep="")
cat(" \n|| .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")
source(paste(path2PepsRScripts, "RNAseqCreateArborescence.R", sep=""))
RNAseqCreateArborescence(path2MainFolder=path2MainFolder)



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*             Transfer data from remote server                   *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#--------------------------------------------------
# Transfer file DataStructure from Canute to Datastructure/ folder
mycode <- paste("`scp -r schorderet@canute.mgh.harvard.edu:", path2expFolderCanute, "DataStructure/Targets.txt " , path2Targets, "`", sep="")
cat(" \n\n\n======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n||\t Opening SSH connection (please enter your password in the terminal window...)", sep="")
cat(" \n|| .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n\n", sep="")
system(mycode)
cat(mycode)
# quikie-mart7

#--------------------------------------------------
# Transfer filtered .bam files found in *path2expFolderCanute*
cat(" \n\n\n======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n||\t Downloading and opening Targets.txt files from remote server", sep="")
cat(" \n|| .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")

if(file.exists(path2Targets)==FALSE){ErrorOutput(paste("No Targets.txt file in\t", path2Targets, sep="")) }
cat(" \n\n Targets file provided: \n\n", sep="")
Targets <- read.delim(path2Targets, comment.char="#")
print(Targets)
mycode <- paste("`scp -r schorderet@canute.mgh.harvard.edu:", path2expFolderCanute, "TopHat/ " , path2bam, "`", sep="")
system(mycode)
cat(mycode)
# Load datasets provided in Targets file
samples <- Targets$FileName
allsamples <- c(levels(Targets$FileName))

#--------------------------------------------------
# At local level
#
cat(" \n\n\n======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n||\t Downloading and reformating files and folders at local level", sep="")
cat(" \n|| .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")

# Move .bam files before deleting entire folder
for(i in 1:length(allsamples)){
  # Rename sample i to *sample1*.bam  
  
  cat("\n Sample: \t",  allsamples[i], "\n", sep="")
  cat("\n\t Downloading \t accepted_hits.bam  from remote server to local", sep="")
  cat("\n\t Renaming \t\t accepted_hits.bam  to \t", allsamples[i], ".bam", sep="")
  renamecode <- paste("`mv ", path2bam, allsamples[i], "/accepted_hits.bam ", path2bam, allsamples[i], "/", allsamples[i], ".bam", "`", sep="")
  system(renamecode)
  # Move sample i to .bam folder
  cat("\n\t Moving \t\t ", allsamples[i], ".bam\t to\t" , path2bam, sep="")
  movecode <- paste("`mv ", path2bam, allsamples[i], "/", allsamples[i], ".bam ", path2bam, "`", sep="")
  system(movecode)
  # Delete rest of folder
  cat("\n\t Deleting \t\t ", path2bam, allsamples[i], sep="")
  delcode <- paste("`rm -r ", path2bam, allsamples[i], "/`", sep="")
  system(delcode)
  cat("\n\n .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")
}

cat(" \n\n\n ======================================================================================", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n||\t Cleaning up local files and folders", sep="")
cat(" \n|| .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n\n", sep="")
cat("\n\n\t ", length(allsamples), " files have been added to \t\t", path2bam, sep="")
cat("\n\t ", length(allsamples), " folders have been deleted from \t", path2bam, sep="")
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
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *||", sep="")
cat(" \n||                                                                                  ||", sep="")
cat(" \n||\t Canute2Local for RNAseq, version ", version," (© Patrick Schorderet 2014)||", sep="")
cat(" \n||                                                                                  ||", sep="")
cat(" \n||\t Terminated with no known error.                                            ||", sep="")
cat(" \n||                                                                                  ||", sep="")
cat(" \n|| * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *||", sep="")
cat(" \n======================================================================================\n\n", sep="")

# Close the R session when run from bash
quit(save = "no", status = 0, runLast = TRUE)

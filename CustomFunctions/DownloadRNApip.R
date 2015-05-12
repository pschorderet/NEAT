

DownloadRNApip <- function(LocalPath2NEAT, LocalPath2NewProject){

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
#                                                               
# LocalPath2NEAT <- "~/NEAT/"
# LocalPath2NewProject <- "~/Desktop/MY_NEW_RNA_PROJECT/"
#*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*                Load and set parameters                         *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#--------------------------------------------------
# Define Paths using user's parameters
LocalPath2CustomFunctions <- paste( LocalPath2NEAT, "CustomFunctions", "/", sep="")
LocalPath2Targets <- paste( LocalPath2NewProject, "DataStructure/Targets.txt", sep="")
LocalPath2bam <- paste( LocalPath2NewProject, "bam/", sep="")
LocalPath2aligned <- paste( LocalPath2NewProject, "Tophat/", sep="")

source(paste(LocalPath2CustomFunctions, "ErrorOutput.R", sep=""))
source(paste(LocalPath2CustomFunctions, "CheckExistenceOfFolder.R", sep=""))
source(paste(LocalPath2CustomFunctions, "CreateArborescenceRNAseq.R", sep=""))

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
  if(length(grep("My_project_title", newLine))==1) { 
    currentline <- gsub("# ", "", newLine); currentline <- gsub("\t", "", currentline); currentline <- gsub("\"", "", currentline);
    ProjectName <- unlist(strsplit(currentline, split = "\\="))[2]
  }
  if(length(grep("Remote_path_to_proj", newLine))==1) { 
    currentline <- gsub("# ", "", newLine); currentline <- gsub("\t", "", currentline); currentline <- gsub("\"", "", currentline);
    RemotePath2MainFolderName <- unlist(strsplit(currentline, split = "\\="))[2]
  }   
}

#------------------------------------------------------------
# Read main folder
cat(paste(" LocalPath2NewProject \t\t ", LocalPath2NewProject, "\n",sep=""))
if(file.exists(LocalPath2NewProject)==TRUE){cat(" \n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n\n ", ProjectName, " exists. \n Reading the Targets.txt file \n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n\n\n ", sep="");}
if(file.exists(LocalPath2NewProject)==FALSE){cat(" \n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n\n Creating\t", ProjectName, "\n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n ", sep="");dir.create(LocalPath2NewProject)}

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

CreateArborescenceRNAseq(path2MainFolder=LocalPath2NewProject)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*             Transfer data from remote server                   *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#--------------------------------------------------
# Transfer file DataStructure from Canute to Datastructure/ folder
#
mycode <- paste("`scp -r ", sshpath, ":", RemotePath2MainFolderName, "/", ProjectName, "/DataStructure/ " , LocalPath2NewProject, "`", sep="")
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
cat(" \n||\t Downloading QC and bam folders", sep="")
cat(" \n|| .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n\n", sep="")

# Download QC folder
cat(" \n Downloading QC folder from remote server", sep="")
mycode <- paste("`scp -r ", sshpath, ":", RemotePath2MainFolderName, "/", ProjectName, "/QC/ " , LocalPath2NewProject, "`", sep="")
system(mycode)
# Download bam folder
#cat(" \n Downloading bam folder from remote server", sep="")
#mycode <- paste("`scp -r ", sshpath, ":", RemotePath2MainFolderName, "/", ProjectName, "/aligned/bam/ ", LocalPath2NewProject, "`", sep="")
#system(mycode)

# Download GRanges folder
cat(" \n Downloading GRanges folder from remote server", sep="")
mycode <- paste("`scp -r ", sshpath, ":", RemotePath2MainFolderName, "/", ProjectName, "/GRangesRData/ ", LocalPath2NewProject, "`", sep="")
system(mycode)
# Download wig folder
cat(" \n Downloading wig folder from remote server", sep="")
mycode <- paste("`scp -r ", sshpath, ":", RemotePath2MainFolderName, "/", ProjectName, "/wig/ ", LocalPath2NewProject, "`", sep="")
system(mycode)




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


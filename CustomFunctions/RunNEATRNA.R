

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                NEAT run RNApip project                         *
#*                                                                *
#******************************************************************
#
# LocalPath2NEAT <- '~/NEAT/'
# LocalPath2Project <- '~/Desktop/MY_NEW_RNA_PROJECT'
#
#
#--------------------------------------------------
# Extrapolate and define paths
#
# Make sure the path is correct
LocalPath2Project <- paste(unlist(strsplit(LocalPath2Project, split="/")), sep="", collapse="/")
NewProjectName <- tail(unlist(strsplit(LocalPath2Project, split="/")), 1)
if(NewProjectName==""){
  NewProjectName <- tail(unlist(strsplit(LocalPath2Project, split="/")), 2)[1]
}

path2Targets <- paste(LocalPath2Project, "/DataStructure/Targets.txt", sep="")

#--------------------------------------------------
# Parse Targets.txt file and save variables
#
res <- readLines(path2Targets)
for(i in 2:length(res)){
 
  newLine <- res[i]
 
  # Store some variables
  if(length(grep("My_personal_ssh", newLine))==1) { 
    currentline <- gsub("# ", "", newLine); currentline <- gsub("\t", "", currentline); currentline <- gsub("\"", "", currentline);
    sshpath <- unlist(strsplit(currentline, split = "\\="))[2]
  }  
  if(length(grep("Remote_path_to_proj", newLine))==1) { 
    currentline <- gsub("# ", "", newLine); currentline <- gsub("\t", "", currentline); currentline <- gsub("\"", "", currentline);
    RemotePath2MainFolderName <- unlist(strsplit(currentline, split = "\\="))[2]
  }  
  if(length(grep("Remote_path_to_NEAT", newLine))==1) { 
    currentline <- gsub("# ", "", newLine); currentline <- gsub("\t", "", currentline); currentline <- gsub("\"", "", currentline);
    RemotePath2NEAT <- unlist(strsplit(currentline, split = "\\="))[2]
  }  
}

#--------------------------------------------------
# Connect to ssh an run pipeline
#
# Check for the existence of the project folder on remote server. If it does not exist, copy it
mycode <- paste("rsync -avz ", LocalPath2Project, " ", sshpath, ":", RemotePath2MainFolderName, "/", sep="")
cat("\n-------------------------------------------\n", sep="")
cat("\nCopying project folder to remote server", sep="")
cat("\n", mycode , "\n\n\n", sep="")
system(mycode)

# Run pipeline
mycode <- paste("ssh ", sshpath, " 'perl ", RemotePath2NEAT, "/RNApip/RNApipRunPipeline.pl ", RemotePath2MainFolderName, "/", NewProjectName, "'", sep="")
cat("\n-------------------------------------------\n", sep="")
cat("\n Running the pipeline", sep="")
cat("\n", mycode , "\n\n\n", sep="")
system(mycode)

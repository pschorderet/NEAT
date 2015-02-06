
CreateNewNEATRNA <- function(LocalPath2NEAT, LocalBasePath2Project, NewProjectName) {
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                NEAT create new RNApip project                  *
#*                                                                *
#******************************************************************
#
# LocalPath2NEAT <- '~/NEAT/'
# LocalBasePath2Project <- '~/Desktop/'
# NewProjectName <- 'TEST'
#
#
#--------------------------------------------------
# Copy the NewProject folder to the new path
#
mycode <- paste("`cp -r ", LocalPath2NEAT, "/RNApip/scripts/NewRNApipProject " , LocalBasePath2Project, "`", sep="")
system(mycode)
mycode <- paste("`mv ", LocalBasePath2Project, "NewRNApipProject " , LocalBasePath2Project, NewProjectName, "`", sep="")
system(mycode)


#--------------------------------------------------
# Extrapolate and define paths
#
path2Targets <- paste(LocalBasePath2Project, NewProjectName, "/DataStructure/Targets.txt", sep="")
path2TargetsTmp <- paste(LocalBasePath2Project, NewProjectName, "/DataStructure/Targets_tmp.txt", sep="")

#--------------------------------------------------
# Parse and modify the Targets.txt file
#
# Create an tmp Targets file
mycode <- paste("`echo '' >", path2TargetsTmp, "`", sep="")
system(mycode)

res <- readLines(path2Targets)
mycode <- paste("`echo '", res[1], "' > ", path2TargetsTmp, "`", sep="")
system(mycode)  
for(i in 2:length(res)){
  newLine <- res[i]
  # Test if line needs a replacement
  if(length(grep("<PROJECT_NAME>", newLine))==1) { newLine <- gsub("<PROJECT_NAME>", NewProjectName, newLine) }
  if(length(grep("<DATE>", newLine))==1) { newLine <- gsub("<DATE>", Sys.Date(), newLine) }
  if(length(grep("<LOCAL_PATH_TO_PROJECT>", newLine))==1) { newLine <- gsub("<LOCAL_PATH_TO_PROJECT>", LocalBasePath2Project, newLine) }
  if(length(grep("<LOCAL_PATH_TO_NEAT>", newLine))==1) { newLine <- gsub("<LOCAL_PATH_TO_NEAT>", LocalPath2NEAT, newLine) }
  
  # Store some variables
  if(length(grep("My_personal_ssh", newLine))==1) { 
    currentline <- gsub("# ", "", newLine); currentline <- gsub("\t", "", currentline); currentline <- gsub("\"", "", currentline);
    sshpath <- unlist(strsplit(currentline, split = "\\="))[2]
  }  
  
  # Re-write the line to the file
  mycode <- paste("`echo '", newLine, "' >> ", path2TargetsTmp, "`", sep="")
  system(mycode)  
#  cat("\n", mycode, sep="")

}
mycode <- paste("`mv ", path2TargetsTmp, " ", path2Targets, "`", sep="")
system(mycode)  

cat("\n\n Please fill in the Targets.txt file before going to the next step. \n\n", sep="")

}
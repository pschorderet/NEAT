InstallNEAT <- function(LocalPath2NEAT, sshPath){
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
#----------------------------------------------------------
#
#    NEAT requires the following R packages
#    To install them, copy paste this  text file into R
#
#----------------------------------------------------------
#
# LocalPath2NEAT <- '~/Desktop/NEAT/'
# sshPath <- "username@serveraddress.edu:~/home/"
#        
  path2Targets <- paste(LocalPath2NEAT, "ChIPpip/scripts/NewChIPpipProject/DataStructure/Targets.txt", sep="")
  path2TargetsTmp <- paste(LocalPath2NEAT, "ChIPpip/scripts/NewChIPpipProject/DataStructure/Targets_tmp.txt", sep="")
  
  #--------------------------------------------------
  RemotePath2NEAT <- strsplit (sshPath, ":")[[1]][2]

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
    # i=2
    newLine <- res[i]
    # Test if line needs a replacement
    if(length(grep("<DATE>", newLine))==1) { newLine <- gsub("<DATE>", Sys.Date(), newLine) }
    if(length(grep("<LOCAL_PATH_TO_NEAT>", newLine))==1) { newLine <- gsub("<LOCAL_PATH_TO_NEAT>", LocalPath2NEAT, newLine) }
    
    if(length(grep("<REMOTE_PATH_TO_NEAT>", newLine))==1) { newLine <- gsub("<REMOTE_PATH_TO_NEAT>", RemotePath2NEAT, newLine) }
    
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
  
  
  
  
  
  
  # Install packages
  cat("\n------------------------------------------------------------", sep='')
  cat("\n Installing various R packages required for NEAT \n\n", sep='')
  
	#install.packages('BiocInstaller')
  package <- c( 'AnnotationDbi', 'caTools', 'edgeR', 'GenomicRanges', 'GenomicAlignments', 'gtools', 'Rsamtools', 'VennDiagram' )  
  source("http://bioconductor.org/biocLite.R")
  BiocInstaller::biocLite(package)
 
 
}
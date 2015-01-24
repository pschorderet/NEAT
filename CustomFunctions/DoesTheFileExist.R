
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  DoesTheFileExist                              *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# DoesTheFileExist takes a path to a file
# return: Checks whether the file exists or not
#         TRUE if all files exists, FALSE if not

DoesTheFileExist <- function(path2file){
  ans <- file.exists(path2file)
  cat(" \n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ", sep="")
  cat("\n Checking existence of files...", sep="")
  TRUES <- 0
  for(i in 1:length(ans)){
    if(ans[i]==TRUE){
      TRUES <- TRUES+1
      # cat("\t", i, ":\t", path2file[i], " exists.\n", sep="")
    }
    if(ans[i]==FALSE){
      cat("\n\t", i, ":\t", path2file[i], " could not be found!", sep="")
    }
  }
  if(TRUES==length(path2file)){
    cat(" \n All files exist", sep="")
    cat(" \n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ", sep="")
    return(TRUE)
  }
  if(TRUES!=length(path2file)){
    cat(" \n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ", sep="")
    cat(" \n Some files could not be found!", sep="")
    cat(" \n **********************************************************\n", sep="")
    return(FALSE)
  }
}



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  CheckExistenceOfFolder                        *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# CheckExistenceOfFolder takes a path to a folder
# return: Checks whether the folder exists, creates folder if not
#         

CheckExistenceOfFolder <- function(path2CheckFolder){
  if(file.exists(path2CheckFolder)==TRUE){
    cat(" \n This file exists in: \t", path2CheckFolder, sep="")
  }
  if(file.exists(path2CheckFolder)==FALSE){
    cat(" \n Creating\t", path2CheckFolder, sep="")
    dir.create(path2CheckFolder)
  }
}
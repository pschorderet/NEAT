
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  RNAseqCreateArborescence                      *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# RNAseqCreateArborescence takes a path to a file
# return: Checks whether the files exists, creates them if not
#         

RNAseqCreateArborescence <- function(path2MainFolder){
  
#   CheckExistenceOfFolder <- function(path2CheckFolder){
#     if(file.exists(path2CheckFolder)==TRUE){
#       cat(" \n This file exists in: \t", path2CheckFolder, sep="")
#     }
#     if(file.exists(path2CheckFolder)==FALSE){
#       cat(" \n Creating\t", path2CheckFolder, sep="")
#       dir.create(path2CheckFolder)
#     }
#   }
  # bam
  path2curfolder <- paste(path2MainFolder, "bam/", sep="")
  CheckExistenceOfFolder(path2CheckFolder=path2curfolder)
  # CountTables
  path2curfolder <- paste(path2MainFolder, "CountTables/", sep="")
  CheckExistenceOfFolder(path2CheckFolder=path2curfolder)
  # DataStructure
  path2curfolder <- paste(path2MainFolder, "DataStructure/", sep="")
  CheckExistenceOfFolder(path2CheckFolder=path2curfolder)
  # GrangesRData
  path2curfolder <- paste(path2MainFolder, "GrangesRData/", sep="")
  CheckExistenceOfFolder(path2CheckFolder=path2curfolder)
  # Logs
  path2curfolder <- paste(path2MainFolder, "Logs/", sep="")
  CheckExistenceOfFolder(path2CheckFolder=path2curfolder)
  # Plots
  path2curfolder <- paste(path2MainFolder, "Plots/", sep="")
  CheckExistenceOfFolder(path2CheckFolder=path2curfolder)
}
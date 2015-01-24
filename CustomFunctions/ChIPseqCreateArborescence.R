
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  ChIPseqCreateArborescence                     *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# ChIPseqCreateArborescence takes a path to a file
# return: Checks whether the files exists, creates them if not
#         

ChIPseqCreateArborescence <- function(path2MainFolder){
  
  checkExistenceOfFolder <- function(path2CheckFolder){
    if(file.exists(path2curfolder)==TRUE){
      cat(" \n This file exists in: \t", path2curfolder, sep="")
    }
    if(file.exists(path2curfolder)==FALSE){
      cat(" \n Creating\t", path2curfolder, sep="")
      dir.create(path2curfolder)
    }
  }
  # bam
  path2curfolder <- paste(path2MainFolder, "bam/", sep="")
  checkExistenceOfFolder(path2CheckFolder=path2curfolder)
  # countTables
  path2curfolder <- paste(path2MainFolder, "CountTables/", sep="")
  checkExistenceOfFolder(path2CheckFolder=path2curfolder)
  # GrangesRData
  path2curfolder <- paste(path2MainFolder, "GrangesRData/", sep="")
  checkExistenceOfFolder(path2CheckFolder=path2curfolder)
  # Logs
  path2curfolder <- paste(path2MainFolder, "Logs/", sep="")
  checkExistenceOfFolder(path2CheckFolder=path2curfolder)
  # Plots
  path2curfolder <- paste(path2MainFolder, "Plots/", sep="")
  checkExistenceOfFolder(path2CheckFolder=path2curfolder)
  # narrowPeakVenns
  path2curfolder <- paste(path2MainFolder, "Plots/narrowPeakVenns/", sep="")
  checkExistenceOfFolder(path2CheckFolder=path2curfolder)
  # broadPeakVenns
  path2curfolder <- paste(path2MainFolder, "Plots/broadPeakVenns/", sep="")
  checkExistenceOfFolder(path2CheckFolder=path2curfolder)
}
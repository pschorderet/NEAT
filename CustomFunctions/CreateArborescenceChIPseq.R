
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  CreateArborescenceChIPseq                     *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# CreateArborescenceChIPseq takes a path to a file
# return: Checks whether the files exists, creates them if not
#         

CreateArborescenceChIPseq <- function(path2MainFolder){
  
  # bam
  path2curfolder <- paste(path2MainFolder, "bam/", sep="")
  CheckExistenceOfFolder(path2CheckFolder=path2curfolder)
  # countTables
  path2curfolder <- paste(path2MainFolder, "CountTables/", sep="")
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
  # narrowPeakVenns
  path2curfolder <- paste(path2MainFolder, "Plots/narrowPeakVenns/", sep="")
  CheckExistenceOfFolder(path2CheckFolder=path2curfolder)
  # broadPeakVenns
  path2curfolder <- paste(path2MainFolder, "Plots/broadPeakVenns/", sep="")
  CheckExistenceOfFolder(path2CheckFolder=path2curfolder)
}
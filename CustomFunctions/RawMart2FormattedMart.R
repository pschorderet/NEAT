
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                  RawMart2FormattedMart                         *
#*                                                                *
#******************************************************************

#------------------------------------------------------------
# RawMart2FormattedMart takes a raw mart object ()
# return: GRanges object of tag coordinates

RawMart2FormattedMart <- function(mymart, strand = "*"){
  head(mymart)
  mymart$MartChr <- mymart$chromosome_name; mymart$MartStart <- mymart$transcript_start; mymart$MartEnd <- mymart$transcript_end; mymart$MartGeneName <- mymart$mgi_symbol; mymart$Strand <- mymart$strand    
  newMart <- mymart[,c("MartChr", "MartStart", "MartEnd", "MartGeneName", "Strand")]
  head(newMart)
  if(strand=="*"){
          cat(" \n ----------------------------------------------------------", sep="")
          cat(" \n Bed file of features on both strands", sep="")
  }
  if(strand=="+"){
    newMart <- newMart[newMart$Strand=='1', ]
          cat(" \n ----------------------------------------------------------", sep="")
          cat(" \n Bed file of features on + strand only", sep="")               
  }
  if(strand=="-"){
    newMart <- newMart[newMart$Strand=='-1', ]
          cat(" \n ----------------------------------------------------------", sep="")
          cat(" \n Bed file of features on - strand only", sep="")           
  
  }
  cat(" \n Length of .bed: \t", length(newMart[,1]), sep="")
  cat(" \n ----------------------------------------------------------\n\n", sep="")
  head(newMart)
  return(newMart)
  
}



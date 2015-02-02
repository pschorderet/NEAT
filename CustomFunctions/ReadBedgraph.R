#------------------------------------------------
# Reads and stores bedgraphs using a path as input
# Returns a dataframe
ReadBedgraph <- function(file) {
  print(paste("File ", file, sep=""))
  dat <- scan(file=file, 
              what=list(character(),integer(),integer(),numeric()), 
              sep="\t", skip=1)
  dat <- data.frame(chr=dat[[1]], start=dat[[2]], 
                    end=dat[[3]], val=dat[[4]])
  return(dat)
}
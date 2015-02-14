#----------------------------------------------------------
#
#    NEAT requires the following R packages
#    To install them, copy paste this  text file into R
#
#----------------------------------------------------------

InstallNEAT <- function(){

  # Install packages
  cat("\n------------------------------------------------------------", sep='')
  cat("\n Installing various R packages required for NEAT \n\n", sep='')
  
	#install.packages('BiocInstaller')
  package <- c( 'AnnotationDbi', 'caTools', 'edgeR', 'GenomicRanges', 'GenomicAlignments', 'gtools', 'Rsamtools', 'VennDiagram' )  
  source("http://bioconductor.org/biocLite.R")
  BiocInstaller::biocLite(package)
 
 
}
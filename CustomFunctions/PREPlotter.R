# PREPlotter plots all annotations that are stored in a tab delimited bed file of form:
# PREPlotter takes the path to the file as input
# chr2	73709289	73709432	grey80	utr

#--------------------------------------------------------------------------------
# Other functions
#----------------------------------------------------------------------------------------
#***************************         box_plotter()               ************************

BoxPlotter <- function(coord, box_heigth_tot = 0.5, middle_gene_coord = 7, box_color_tot = "Black", fill_color = "white", box_lwd=0){
  Boxbox_lwd <- 0
  rect(coord[1], (middle_gene_coord-(0.5*box_heigth_tot)), coord[2], (middle_gene_coord+(0.5*box_heigth_tot)), col=fill_color, border=NA, lwd = box_lwd)
}


#--------------------------------------------------------------------------------
# PREPlotter
PREPlotter <- function(pathToGenes, ypos, boxwidth=0.2, colorPRE="Orange", inviewcoord){
  bedfile <- read.table(pathToGenes,as.is =T, header=T)
  #-----------------------------
  #  Different paramaters
  UTRw <- (boxwidth-(boxwidth*0.4))
  size_of_text <- 0.5
  #-----------------------------
  #  Find how many genes there are
  genenumber  <- which(bedfile[,1]==">" | bedfile[,1]=="<")
  #-----------------------------
  #  For each gene, parse file and plot types
  for(h in 1:nrow(bedfile)){  
    BoxPlotter(as.matrix(bedfile[h,2:3]), middle_gene_coord=ypos, fill_color=colorPRE, box_heigth_tot=UTRw)
  }
}

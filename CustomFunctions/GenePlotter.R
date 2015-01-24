# GenePlotter plots all genes that are stored in a tab delimited bed file of form:
# GenePlotter takes the path to the file as input
#  Gene	coordinates:
# >	73709289	73712084	-	Atp5g3
# chr2	73709289	73709432	grey80	utr
# chr2	73709433	73709547	grey70	exon
# chr2	73710023	73710218	grey70	exon
# chr2	73710703	73710781	grey70	exon
# chr2	73711794	73711832	grey70	exon
# chr2	73711833	73712084	grey80	utr

#--------------------------------------------------------------------------------
# Other functions
#----------------------------------------------------------------------------------------
#***************************         box_plotter()               ************************

BoxPlotter <- function(coord, box_heigth_tot = 0.5, middle_gene_coord = 7, box_color_tot = "Black", fill_color = "white", box_lwd=0){
  Boxbox_lwd <- 0
  rect(coord[1], (middle_gene_coord-(0.5*box_heigth_tot)), coord[2], (middle_gene_coord+(0.5*box_heigth_tot)), col=fill_color, border=NA, lwd = box_lwd)
}
PlotUTR <- function(geneline, yposition, color, heigth){
  BoxPlotter(coord=geneline[2:3], box_heigth_tot=heigth, middle_gene_coord=yposition, box_color_tot=color, fill_color=color, box_lwd=0.5)  
}
PlotExon <- function(geneline, yposition, color, heigth){
  BoxPlotter(coord=geneline[2:3], box_heigth_tot=heigth, middle_gene_coord=yposition, box_color_tot=color, fill_color=color, box_lwd=0.5)
}
PlotIntron <- function(geneline, yposition, color, heigth, strand){
  ngeneline  <- as.numeric(geneline[2:3])
  if(strand=="-"){heigth <- -(heigth)}
  segments(ngeneline[1], yposition, (ngeneline[1] + ((ngeneline[2]-ngeneline[1])/2)), (yposition + (0.5*heigth)), col= color, lwd=0.5)
  segments((ngeneline[1] + ((ngeneline[2]-ngeneline[1])/2)), (yposition + (0.5*heigth)), ngeneline[2], yposition, col= color, lwd=0.5)
}


#--------------------------------------------------------------------------------
# GenePlotter
GenePlotter <- function(pathToGenes, ypos, boxwidth=0.2, colorExon="grey30", colorUTR="grey60", inviewcoord){
  bedfile <- read.table(pathToGenes,as.is =T, header=F)
  #-----------------------------
  #  Different paramaters
  UTRw <- (boxwidth-(boxwidth*0.4))
  size_of_text <- 0.5
  #-----------------------------
  #  Find how many genes there are
  genenumber  <- which(bedfile[,1]==">" | bedfile[,1]=="<")
  #-----------------------------
  #  For each gene, parse file and plot types
  for(currentgene in 1:(length(genenumber)-1)){
      currentline <- genenumber[currentgene]
      currentstrand <- bedfile[currentline, 4]
    #  text(mean(as.numeric(bedfile[currentline,2:3])), (ypos+0.4*boxwidth), bedfile[currentline,5], cex=0.5, srt=0)
      while(currentline != genenumber[currentgene+1]){ 
      
        currentline <- currentline+1
        if(bedfile[currentline,5]=="utr"){PlotUTR(geneline=bedfile[currentline,], yposition=ypos, color=colorUTR, heigth=UTRw)}
        if(bedfile[currentline,5]=="exon"){PlotExon(geneline=bedfile[currentline,], yposition=ypos, color=colorExon, heigth=boxwidth)}                                                
        if(bedfile[currentline,5]=="intron"){PlotIntron(geneline=bedfile[currentline,], yposition=ypos, color=colorExon, heigth=boxwidth, strand=currentstrand)} 
      }
  }
  # Draw line between each gene
  genenumber <- genenumber[-(length(genenumber))]
  intergeneind <- cbind(c(inviewcoord[1],as.numeric(bedfile[genenumber,2])),c(as.numeric(bedfile[genenumber,3]),inviewcoord[2]))
  segments(intergeneind[,1], ypos, intergeneind[,2], ypos, col=colorUTR, cex=0.5)
}

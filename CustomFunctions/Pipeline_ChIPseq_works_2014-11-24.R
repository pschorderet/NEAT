
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                ChIPseq R pipeline for local                    *
#*                                                                *
#******************************************************************
#
# Connect to the ssh
# ssh -X ps607@ssh.research.partners.org
# Secretcode3
# ssh schorderet@canute.mgh.harvard.edu
# quikie-mart7

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
#                                                                 #
# This script is intended to be run from bash/shell Mac terminal  #
#                                                                 #
# To check if the script works locally, please  run these         #
# parameters first:                                               #
#     path2RScripts <- "~/PepsRscripts/RScripts/"                 #
#     MainFolder <- "ChIP_arborescence/"; path2MainFolder <- paste(path2RScripts, MainFolder, sep="")
#     MainFolder <- "DIPG_2014-07-07_ChIPseq/"; path2MainFolder <- paste(path2RScripts, MainFolder, sep="")
#     MainFolder <- "DIPG_2014-11-18_ChIPseq/"; path2MainFolder <- paste(path2RScripts, MainFolder, sep="")
#     MainFolder <- "DIPG_consolidated/"; path2MainFolder <- paste(path2RScripts, MainFolder, sep="")
#     nameOfBed <- "mm9_Hox.bed" ; binNumber = 10 ; strand <- "+" ; normInp <- FALSE
#     nameOfBed <- "mm9_TSS_10kb.bed" ; binNumber = 10 ; strand <- "+" ; normInp <- FALSE
#                                                                 #
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #

version <- "0.1.3 Oct 2014"

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arborescence of files:                                          #
#                                                                 #
# *...* means the name is variable                                #
#                                                                 #
# ~/PepsRscripts/RScripts/                                        #
#                                                                 #
#     └─ *MainFolder*                                             #
#           └─ bam                                                #
#                 └─ *samples1*.bam                               #
#                 └─ ...                                          #
#                                                                 #
#           └─ DataStructure                                      #
#                 └─ Targets.txt                                  #
#                                                                 #
#           └─ GrangesRData                                       #
#                 └─ *samples1*.bam.GRangesRData                  #
#                 └─ ....                                         #
#                                                                 #
#           └─ CountTables                                        #
#                 └─ CountTable_*samples1*_*martObject*_*date*.bed#
#                 └─ ....                                         #
#                                                                 #
#           └─ Plots                                              #
#                 └─ *FactorN*_*martObject*_*strand*.pdf          #
#                 └─ ....                                         #
#                                                                 #
#     └─ MartObject                                               # 
#           └─ mm9_Trasncripts.bed                                #
#           └─ ...                                                #
#                                                                 #
#     └─ ReferenceFiles                                           #
#           └─ mm9KG.sqlite                                       #
#           └─ Mus_musculus.NCBIM37.67.allchr.gtf                 #
#           └─ UC_ID_to_GENE_symbol_mm9.txt                       #
#           └─ ...                                                #
#                                                                 #
#     └─ PepsFunctions                                            #
#           └─ ...                                                #
#                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*                Load and set parameters                         *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#--------------------------------------------------
# Load R and different packages
library(Rsamtools); library(GenomicRanges); library(GenomicAlignments);  library(VennDiagram);

# Define paths to subfolders stored in RScripts
path2bam <- paste(path2MainFolder, "bam/", sep="")
path2GRangesRData <- paste(path2MainFolder, "GRangesRData/", sep="")
path2CountTables <- paste(path2MainFolder, "CountTables/", sep="")
path2Plots <- paste(path2MainFolder, "Plots/", sep="")
path2Logs <- paste(path2MainFolder, "Logs/", sep="")
path2Mart <- paste(path2RScripts, "MartObjects/", sep="")
path2ReferenceFiles <- paste(path2RScripts,"ReferenceFiles/", sep="")
path2PepsRScripts <- paste(path2RScripts,"PepsFunctions/", sep="")
path2Targets <- paste(path2MainFolder, "DataStructure/Targets.txt", sep="")

# Source pepsfunctions
source(paste(path2PepsRScripts, "Bam2GRangesRData.R", sep=""))
source(paste(path2PepsRScripts, "Bed2GRanges.R", sep=""))
source(paste(path2PepsRScripts, "CountOverlaps2matrix.R", sep=""))
source(paste(path2PepsRScripts, "DoesTheFileExist.R", sep=""))
source(paste(path2PepsRScripts, "ErrorOutput.R", sep=""))
source(paste(path2PepsRScripts, "LoadMartGRanges.R", sep=""))
source(paste(path2PepsRScripts, "Mart2GRanges.R", sep=""))
source(paste(path2PepsRScripts, "OutputNumberOfReadsFromGRanges.R", sep=""))
source(paste(path2PepsRScripts, "RawMart2FormattedMart.R", sep=""))
source(paste(path2PepsRScripts, "Tags2GRanges.R", sep=""))



#------------------------------------------------------------
# Redirect output to log file
if(file.exists(path2Logs)==TRUE){cat(" \n Logs file exists.\n\n", sep="")}
if(file.exists(path2Logs)==FALSE){cat(" \n Creating\t", path2Logs, "\n\n", sep="");dir.create(path2Logs)}
sink(paste(path2Logs, Sys.Date(), "_", format(Sys.time(), "%X"), ".txt", sep=""), append=FALSE, split=FALSE)
#sink()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*                Create proper arborescence                      *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Paste variables in log file
cat(" \n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  ", sep="")
cat(" \n ||", sep="")
cat(" \n ||\t ChIPseq analysis version ", version, " (© Patrick Schorderet)", sep="")
cat(" \n ||", sep="")
cat(" \n ||\t Path to MainFolder : \t\t", path2MainFolder, sep="")
cat(" \n ||\t Path to PepsRscripts : \t", path2RScripts, sep="")
cat(" \n ||\t Features to align to : \t", nameOfBed, sep="")
cat(" \n ||\t No of bins per feature : \t", binNumber, sep="")
cat(" \n ||\t Strand specificity : \t\t", strand, sep="")
cat(" \n ||\t Compute venn diagrams : \t", Venn, sep="")
cat(" \n ||\t Normalize to input : \t\t", normInp, sep="")
cat(" \n ||", sep="")
cat(paste(" \n ||\t Run date:\t\t\t\t\t", Sys.Date(), sep=""), sep="")
cat(" \n ||", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ======================================================================================", sep="")

# Create proper arborescence
cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Create / check arborescence of ", paste(path2MainFolder, sep=""), sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")
source(paste(path2PepsRScripts, "ChIPseqCreateArborescence.R", sep=""))
ChIPseqCreateArborescence(path2MainFolder=path2MainFolder)

#--------------------------------------------------
# Define function to plot axes
PlotAxes <- function(mat, binrange=0, yrange=0, strand="+"){
  if(length(colnames(mat))==0){
    if(strand=="+"){axis(1, at=seq(from=1, to=tail(binrange,1), length.out=2), labels=c("TSS", "TES"))}
    if(strand=="-"){axis(1, at=seq(from=1, to=tail(binrange,1), length.out=2), labels=c("TES", "TSS"))}
  }
  if(length(colnames(mat))!=0){
    colnames(mat[binNumber/2])
    
    if(length(binrange)==1){}
    if(length(colnames(mat))!=1){
      myseq <- floor(seq(1, length(colnames(mat)), length.out=4))
      axis(1, at=myseq, labels=colnames(mat)[myseq])    
      abline(v=binrange[length(binrange)/2])
    }  
  }  
  # y axis
  if(length(yrange)>1){axis(2, at=yrange, labels=round(yrange,5))}  
}

#------------------------------------------------------------
# Define parameters
chromosomes = c(paste("chr", seq(1,19),  sep = ""), "chrX", "chrY")


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                                 #
#                  Read Targets.txt                               #
#                                                                 #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Read Targets file", sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")
# Check if Targets file exists
if(file.exists(path2Targets)==FALSE){ErrorOutput(paste("No Targets.txt file in\t", path2Targets, sep="")) }

cat(" \n\n Targets file provided: \n\n", sep="")
Targets <- read.delim(path2Targets, comment.char="#")
print(Targets)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                                 #
#                  Load GRanges data                              #
#                                                                 #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Load datasets provided in Targets file
samples <- Targets$FileName
allsamples <- c(levels(Targets$FileName), levels(Targets$InpName))
# Remove '-' when no input is provided
if(length(which(allsamples=="-")) != 0){
  allsamples <- allsamples[-which(allsamples=="-")]    
}


#------------------------------------------------------------
# Check existence of .bam.GRanges.RData or .bam files
cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Create / load GRanges object from bam files", sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")

# If .bam.GRanges.RData files DO NOT exists, create them
bamGRangesRDataPath <- paste(paste(path2GRangesRData, allsamples, sep=""), ".bam.GRanges.RData", sep="")
if(DoesTheFileExist(path2file=bamGRangesRDataPath)!=TRUE){
  # If .bam files cannot be found, abort
  bamPath <- paste(path2file=path2bam, allsamples, ".bam", sep="")
  if(DoesTheFileExist(path2file=bamPath)!=TRUE){
    ErrorOutput("Some .bam files could not be found")
  }
  # If .bam files are present, create .bam.GRanges.RData files
  if(DoesTheFileExist(path2file=bamPath)==TRUE){
    # Create GRangesRData files
    cat(" \n\n **********************************************************\n", sep="")
    cat(" \n Generating .bam.GRanges.RData files", sep="")
    for(k in 1:length(bamPath)){    
      cat(" \n\n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
      cat(" \n ", k, " / ", length(bamPath) ," :\t", bamPath[k], sep="")
      Bam2GRangesRData(path2bam=bamPath[k])
    }  
  }
}

# Check if .bam.GRanges.RData files exists, load them
GRangesSamples <- list.files(path2GRangesRData, pattern = "*.bam.GRanges.RData$")  
if(DoesTheFileExist(path2file=paste(paste(path2GRangesRData, allsamples, sep=""), ".bam.GRanges.RData", sep=""))==TRUE){
  # This doesn't work if put in a function!!  
  # Load each sample into its .GRanges named object
  for(i in 1:length(GRangesSamples)){
    path2load <- paste(path2GRangesRData, GRangesSamples[i], sep="")
    load(path2load)
    assign(GRangesSamples[i], GRangesObj)
  }
}

#------------------------------------------------------------
# Quality Control 
GRallsamples <- sort(ls()[grep(ls(), pattern="bam.GRanges.RData$")]) 
myGRsamples <- GRallsamples[which(sub(".bam.GRanges.RData$", "", GRallsamples)%in%Targets$FileName)]
OutputNumberOfReadsFromGRanges(GRallsamples)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*                Load reference files/genomes                    *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#------------------------------------------------------------
# Create GRanges objects contaninig bins around all feature
#     Use the Features2Mart script to create this object !!

# Load Mart object into GRanges
path2MartObject <- paste(path2Mart, nameOfBed, sep="")
# Load Mart object & store gene names
cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Create / load Mart file as GRanges object", sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")

GRfeat <- LoadMartGRanges(path2MartObject=path2MartObject, binNumber=binNumber, strand=strand)
geneNames <- GRfeat$GeneNames[seq(from=1, to=length(GRfeat), by=binNumber)]

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#*                                                                *
#*                  Count overlaps                                *
#*                                                                *
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# For each GRanges object:  - Compute overlap and store countTables
#                           - Plot average coverage and store in pdf
cat(" \n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ||\t Compute overlaps and generate count tables", sep="")
cat(" \n || .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n", sep="")

# Open pdf to store the plots
pdf(paste(path2Plots, sub(".bed", "", nameOfBed), "_", strand, ".pdf", sep=""), paper='USr', useDingbats=FALSE, width = 27, height = 21) 

# Compute and store countTables, pdf, etc
cols <- c("Black", "Blue", "Red", "Green", "Yellow", "Purple", "Cyan", "Grey")
for(i in levels(Targets$Factor)){
  # (i=levels(Targets$Factor)[1])
  # (i=levels(Targets$Factor)[2])
  # (i=levels(Targets$Factor)[3])
  cat(" \n\n **********************************************************", sep="")  
  cat(" \n *", sep="")
  cat(" \n *\tFactor: \t", i, sep="")
  cat(" \n *", sep="")
  cat(" \n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
  
  sameFactorSamples <- Targets$FileName[which(Targets$Factor==i)]
  
  # For each sample with the same Factor, compute overlap and plot
  activeCTNames <- activeCTNamesRep <- NULL
  for(j in sameFactorSamples){
    # (j=sameFactorSamples[1])
    # (j=sameFactorSamples[2])
    # (j=sameFactorSamples[3])
    cat(" \n\n Sample: \t", j, "\n\n", sep="")    
    
    #--------------------------------------------------------
    #----------------   Replicates   ------------------------
    # If there is a Replicate corresponding to the current sample, average the countTable
    currentDupNum <- Targets$Replicate[which(Targets$FileName == j)]
    possibleDupIndex <- NULL
    for(y in sameFactorSamples){
      possibleDupIndex <- c(possibleDupIndex, Targets$Replicate[which(Targets$FileName == y)])
    }
    #DupNumSameFactorSamples <- Targets$Replicate[which(Targets$FileName == sameFactorSamples)]
    dupindex <- which(possibleDupIndex==currentDupNum)
    
    # 3 Replicates
    if(length(dupindex)==3){ ErrorOutput("Pipeline cannot deal with more than 2 replicates per condition.") }
    if(length(dupindex)==2){      
      cat(" \n Replicate found:\n", sep="")
      for(l in sameFactorSamples[dupindex]){cat("\t\t Rep \t", l, "\n", sep="")}
    }
    
    if(Venn==TRUE){
      #--------------------------------------------------------
      # Load broadPeaks, narrowPeak
      #
      # broadPeak
      #
      path2currentSamplebroakPeak <- paste(path2MainFolder, "broadPeak/", j, ".broadPeak", sep="")
      if(file.exists(path2currentSamplebroakPeak)==FALSE){
        cat("\t", j, ".broadPeak file does not exist.\n", sep="")
      }
      if(file.exists(path2currentSamplebroakPeak)==TRUE){
        broadPeak <- read.table(path2currentSamplebroakPeak, header=FALSE)
        colnames(broadPeak) <- c("BedChr", "BedStart", "BedEnd", "-", "-", "-", "-", "BedStrand", "-")
        broadPeakGR <- Bed2GRanges(bedReply=broadPeak[, c("BedChr", "BedStart", "BedEnd", "BedStrand")])
        namebroadPeakGR <- paste(j, ".broadPeak.GRanges", sep="")
        assign(namebroadPeakGR, broadPeakGR)  
      }
      #
      # narrowPeak
      #
      path2currentSamplenarrowPeak <- paste(path2MainFolder, "narrowPeak/", j, ".narrowPeak", sep="")
      if(file.exists(path2currentSamplenarrowPeak)==FALSE){
        cat("\t", j, ".narrowPeak file does not exist.\n", sep="")
      }
      if(file.exists(path2currentSamplenarrowPeak)==TRUE){
        narrowPeak <- read.table(path2currentSamplenarrowPeak, header=FALSE)
        colnames(narrowPeak) <- c("BedChr", "BedStart", "BedEnd", "-", "-", "-", "-", "BedStrand", "-")
        narrowPeakGR <- Bed2GRanges(bedReply=narrowPeak[, c("BedChr", "BedStart", "BedEnd", "BedStrand")])
        namenarrowPeakGR <- paste(j, ".narrowPeak.GRanges", sep="")
        assign(namenarrowPeakGR, narrowPeakGR)  
      }
    }
    # end of if(Venn==TRUE)
    
    #--------------------------------------------------------
    # Replicates
    for(k in sameFactorSamples[dupindex]){
      #k=sameFactorSamples[dupindex][1]
      #k=sameFactorSamples[dupindex][2]

      nameCountTable <- paste("CountTable_", k, "_", nameOfBed, sep="")
      currentGR <- get(paste(k, ".bam.GRanges.RData", sep=""))
      countTable <- CountOverlaps2matrix(GRanges1=GRfeat, GRanges2=currentGR, binNumber=binNumber, geneNamesSingle=geneNames)
      
      # Normalize to input here
      if(normInp=="TRUE"){
        cat(" \n\n Inp sample: \t", as.character(Targets$InpName[which(Targets$FileName==k)]), "\n", sep="")
        nameCountTable <- paste("CountTable_Norm2Inp_", k, "_", nameOfBed, sep="")
        # Find corresponding input
        currentGRInp <- get(paste(Targets$InpName[which(Targets$FileName==k)], ".bam.GRanges.RData", sep=""))
        countTableInp <- CountOverlaps2matrix(GRanges1=GRfeat, GRanges2=currentGRInp, binNumber=binNumber, geneNamesSingle=geneNames)        
        cat(" \n Normlizing to input sample", sep="")
        countTable <- countTable-countTableInp
      }
      
      # Assign countTable to the nameCountTable
      assign(nameCountTable, countTable)
      # Save count table to folder CountTables
      cat(" \n Store countTable in : \t\t", paste(path2CountTables, sub(".bed", "", nameCountTable), "_bins", binNumber, "_", Sys.Date(), ".bed", sep=""), sep="")    
      write.table(get(nameCountTable), paste(path2CountTables, sub(".bed", "", nameCountTable), "_bins", binNumber, "_", Sys.Date(), ".bed", sep=""), sep = "\t", row.names = TRUE, col.names=TRUE, quote=FALSE, na="")
      # Store name of all active count tables in activeCTNames
      activeCTNames <- c(activeCTNames, nameCountTable)
    } # end of for(k in sameFactorSamples[dupindex])

    # Remove duplicated activeCTNames
    activeCTNames <- activeCTNames[!duplicated(activeCTNames)]
    activeCTNames
    
    #--------------------------------------------------------
    # Compute average countTables of 2 replicates
    #
    if(length(dupindex)==2){
      cat("\n\n +++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", sep="")
      cat(" \n Compute average of replicate countTables", sep="")    
      mylist <- list(get(activeCTNames[dupindex[1]]), get(activeCTNames[dupindex[2]]))
  #    str(mylist)
      dup1matrix <- simplify2array(mylist[[1]])
      dup2matrix <- simplify2array(mylist[[2]])
      # Compute the average between both duplicates
      averagedCountTable <- apply(simplify2array(mylist), 1:2, mean)
      # Plot the correlation per bin of each duplicate
      plotTitle <- paste("Correlation plot between replicate:\n", sameFactorSamples[dupindex[1]], "\n", sameFactorSamples[dupindex[2]], sep="")  
      plot(diag(cor(dup1matrix, dup2matrix)), xlab="bins", ylab="", type="b", ylim=c(-1,1), main=plotTitle, cex=0.5, cex.axis=0.8, cex.main=0.8)
      abline(h=0, lty=2)
      abline(h=mean(diag(cor(dup1matrix, dup2matrix))), lty=4, col="Red")
      
      # Assign countTable to the nameCountTable    
      nameCountTable <- paste("CountTable_AverageDup_", paste(sameFactorSamples[dupindex], collapse="--"), "_", nameOfBed, sep="")
      assign(nameCountTable, averagedCountTable)
      #activeCTNames <- c(activeCTNames, nameCountTable)
      activeCTNamesRep <- c(activeCTNamesRep, nameCountTable)
      # Save count table to folder CountTables        
      cat(" \n Store averaged countTable of replicate in : \t", nameCountTable, "_bins", binNumber, "_", Sys.Date(), ".bed", sep="")    
      write.table(get(nameCountTable), paste(path2CountTables, nameCountTable, sep=""), sep = "\t", row.names = TRUE, col.names=TRUE, quote=FALSE, na="")      
      
      cat("\n\n +++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n ", sep="")
    } # end of if
  } # end of j (go through samples that have the same factor)
  
  # If there is more than one one samples in sameFactorSamples (independently of whether there are replicates...)
  #     Consolidate all active countTable names
  if(length(sameFactorSamples)>1){
    currentTargets <- Targets[which(Targets$Factor==i),]
    df <- currentTargets$Replicate
    keepIndex <- which(currentTargets$Replicate == df[!(duplicated(df) | duplicated(df, fromLast = TRUE))])  
    activeCTNames <- c(activeCTNames[keepIndex], unique(activeCTNamesRep))
  }
  
  #--------------------------------------------------------
  # Prepare and plot metagenomic data
  #
  # Using countTable of the activeCTNames, find min and max range
  #
  rg <- range(colSums(get(activeCTNames[1])))
  miny <- rg[1]; maxy <- rg[2]
  for(k in 1:length(activeCTNames)){
    rg <- range(colSums(get(activeCTNames[k])))
    if(rg[1] < miny){miny <- rg[1]}
    if(rg[2] > maxy){maxy <- rg[2]}
  }

  # Clean activeCTNames to find FileShort of corresponding plot
  activeCTNames <- unique(activeCTNames)
  newActiveCTName <- sub("CountTable_", "", activeCTNames[1]); newActiveCTName <- sub("Norm2Inp_", "", newActiveCTName); newActiveCTName <- sub(paste("_", nameOfBed, sep=""), "", newActiveCTName); newActiveCTName <- sub(".bed", "", newActiveCTName)
  legends <- NULL
  plot(colSums(get(activeCTNames[1])), type="l", col=cols[1], lwd=3, axes=F, xlab="Length [bp]", ylab="", ylim=c(miny, maxy), main=paste(i, nameOfBed," ; strand ", strand, sep= " "))
  for(h in 1:length(activeCTNames)){
    newActiveCTName <- sub("CountTable_", "", activeCTNames[h]); newActiveCTName <- sub("Norm2Inp_", "", newActiveCTName); newActiveCTName <- sub(paste("_", nameOfBed, sep=""), "", newActiveCTName); newActiveCTName <- sub(".bed", "", newActiveCTName); newActiveCTName <- sub("AverageDup_", "", newActiveCTName)
    newActiveCTName <- strsplit(newActiveCTName, "--")[[1]][1]
    legends <- c(legends,as.character(Targets$FileShort[which(Targets$FileName == newActiveCTName)]))
    points(colSums(get(activeCTNames[h])), type="l", col=cols[h], lwd=3)  
  }
  box()
  PlotAxes(mat=get(activeCTNames[1]), yrange=c(miny, maxy), strand=strand)  
  legend(1, maxy, legend=legends, lty=1, lwd=3, col=cols)
  
  
  #--------------------------------------------------------
  # Prepare and plot Venn diagrams
  #
  # broadPeak venns  
  if(Venn==TRUE){
    
    broadPeakAll <- ls()[grep(".broadPeak.GRanges", ls())]
    if(length(broadPeakAll)!=0){
      path2VennsBroad <- paste(path2Plots, "broadPeakVenns/", sep="")
      cat("\n\n broadPeakAll: ", broadPeakAll, sep="\n")
      cat("\n Computing Venn diagrams of overlaps on broadPeak", sep="")
      mycols=c("Red", "Blue")
      for(v in 1:length(broadPeakAll)){
        for(w in 1:length(broadPeakAll)){
          if(v!=w){
            GRA <- get(broadPeakAll[v])
            GRB <- get(broadPeakAll[w])
            broadPeakNames <- c(broadPeakAll[v], broadPeakAll[w])
            AoverB <- sum(countOverlaps(GRA, GRB)!=0)
            venn.plot <- venn.diagram(list(1:length(GRA), (length(GRA)-AoverB):((length(GRA)-AoverB)+length(GRB))),  fill=mycols, alpha=rep(0.3, 2), cex = 2, cat.fontface=4, category.names=c("", ""), main=paste("Overlapping peaks", "\nRed = ", broadPeakNames[1], "\nBlue = ", broadPeakNames[2], sep=""), filename=paste(path2VennsBroad, broadPeakNames[1], "--", broadPeakNames[2], ".tiff", sep=""), cat.dist = rep(0.1, 2)) 
          } # end if
        } # end w
      } # end v  
      # Remove .broadPeak.GRanges objects
      rm(list=ls()[grep(".broadPeak.GRanges", ls())])
    } # end broadPeak venns  
    
    # narrowPeak venns  
    narrowPeakAll <- ls()[grep(".narrowPeak.GRanges", ls())]
    if(length(narrowPeakAll)!=0){
      path2VennsNarrow <- paste(path2Plots, "narrowPeakVenns/", sep="")
      cat("\n\n narrowPeakAll: ", narrowPeakAll, sep="\n")
      cat("\n Computing Venn diagrams of overlaps on narrowPeak", sep="")
      for(v in 1:length(narrowPeakAll)){
        for(w in 1:length(narrowPeakAll)){
          if(v!=w){
            GRA <- get(narrowPeakAll[v])
            GRB <- get(narrowPeakAll[w])
            narrowPeakNames <- c(narrowPeakAll[v], narrowPeakAll[w])
            AoverB <- sum(countOverlaps(GRA, GRB)!=0)  
            venn.plot <- venn.diagram(list(1:length(GRA), (length(GRA)-AoverB):((length(GRA)-AoverB)+length(GRB))),  fill=mycols, alpha=rep(0.3, 2), cex = 2, cat.fontface=4, category.names=c("", ""), main=paste("Overlapping peaks", "\nRed = ", broadPeakNames[1], "\nBlue = ", broadPeakNames[2], sep=""), filename=paste(path2VennsNarrow, narrowPeakNames[1], "--", narrowPeakNames[2], ".tiff", sep=""), cat.dist = rep(0.1, 2)) 
          } # end if
        } # end w
      } # end v  
      # Remove .broadPeak.GRanges objects
      rm(list=ls()[grep(".narrowPeak.GRanges", ls())])
    } # end narrowPeak venns  
  } # end Venn==TRUE
  
} # end of for(i in levels(Targets$Factor))
dev.off()






#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                                 #
#                  Session info and closure                       #
#                                                                 #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Print sessionInfo()
cat(" \n\n\n || -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-   ", sep="")
cat(" \n\n", sep="")
print(sessionInfo())
cat(" \n\n", sep="")
cat(" \n || -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n\n   ", sep="")

cat(" \n\n\n ======================================================================================", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  ", sep="")
cat(" \n ||", sep="")
cat(" \n ||\t ChIPseq analysis version ", version," (© Patrick Schorderet 2014)", sep="")
cat(" \n ||", sep="")
cat(" \n ||\t Terminated with no known error.", sep="")
cat(" \n ||", sep="")
cat(" \n || * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ", sep="")
cat(" \n ======================================================================================\n\n", sep="")

# Close the R session when run from bash
quit(save = "no", status = 0, runLast = TRUE)

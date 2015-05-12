
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                      PeakCalling                               *
#*                                                                *
#******************************************************************
#*                                                                *
#*----------------------------------------------------------------*
#*                                                                *
#* The below peakcalling script runs MACS 	                  *
#*                                                                *
#*----------------------------------------------------------------*

cat("\n ----------------------------------------\n", sep="")
cat("\n Entered PeakCaller_MACS.R \n", sep="")
cat("\n ----------------------------------------\n\n", sep="")

#*----------------------------------------------------------------*
# Initialize some parameters

args		<- commandArgs(TRUE)
code		<- args[0]		# path to peakcalling.R
path2expFolder	<- args[1]	# path to the main experiment folder
path2sample	<- args[2]		# path to sample
path2input	<- args[3]		# path to input 
fdr		<- as.numeric(args[4]) 	# fdr 1e-09
makeunique	<- args[5]		# SPP specific
removepcr	<- args[6]		# SPP specific

#-----------------------------------------------------------
# LOAD TARGETS

Path2Targets <- paste(path2expFolder, "/DataStructure/Targets.txt", sep="")
cat("\n Path2Targets: \t", Path2Targets, "", sep="")
res <- readLines(Path2Targets)
for(i in 1:length(res)){
	newLine <- res[i]
	# Store some variables

	if(length(grep("Reference_genome", newLine))==1) { 
		currentline <- gsub("# ", "", newLine); currentline <- gsub("\t", "", currentline); currentline <- gsub("\"", "", currentline);
		genome <- unlist(strsplit(currentline, split = "\\="))[2]
		if( length( grep("hs", genome)!=0) ) { genomeMACS = "hs" }
		if( length( grep("mm", genome)!=0) ) { genomeMACS = "mm" }
		if( length( grep("dm", genome)!=0) ) { genomeMACS = "dm" }
		cat("\n genomeMACS: \t", genomeMACS, sep="")
	}
}

#-----------------------------------------------------------
# LOAD LIBRARY


#*----------------------------------------------------------------*
# Find name of sample and input from their resp. paths

mySample	<- tail(unlist(strsplit(path2sample, split="/")),2)[1]
myInput 	<- tail(unlist(strsplit(path2input, split="/")),2)[1]

#*----------------------------------------------------------------*
# Define paths 
path2peakcalling	<- paste( path2expFolder, "/peakcalling/", sep="" )
path2currSampleFolder	<- paste( path2peakcalling, "/", mySample, "/", sep="" )
path2currInputFolder	<- paste( path2peakcalling, "/", myInput, "/", sep="" )


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                    MACS run PeakCalling                        *
#*                                                                *
#******************************************************************

# Samples
#-----------------------------------------------------------  
cat("\n Peak calling for sample: \t", mySample, "\t over \t", myInput, "\n", sep="")

#-----------------------------------------------------------

#*----------------------------------------------------------------*
# select informative tags based on the binding characteristics
cmd <- paste("`macs14 callpeak -t ", path2sample," -c ", path2input," -f BAM -g ", genomeMACS, " -n ", path2peakcalling, "/narrowPeak/", mySample," -p ", fdr, "`", sep="")
cat("\n cmd line: \t ", cmd, "\n\n", sep="")
system(cmd)




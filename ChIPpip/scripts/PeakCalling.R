
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                      PeakCalling                               *
#*                                                                *
#******************************************************************
#*                                                                *
#*----------------------------------------------------------------*
#*                                                                *
#* The below peakcalling script is based on SPP                   *
#*                                                                *
#*----------------------------------------------------------------*

cat("\n ----------------------------------------\n", sep="")
cat("\n Entered PeakCalling.R \n", sep="")
cat("\n ----------------------------------------\n\n", sep="")

#*----------------------------------------------------------------*
# Initialize some parameters

args		<- commandArgs(TRUE)
code		<- args[0]		# path to peakcalling.R
path2expFolder	<- args[1]		# path to the main experiment folder
path2sample	<- args[2]		# path to sample
path2input	<- args[3]		# path to input 
fdr		<- as.numeric(args[4]) 	# fdr
makeunique	<- args[5]		# SPP specific
removepcr	<- args[6]		# SPP specific

#-----------------------------------------------------------
# LOAD LIBRARY

library(spp)
#library(snow)
#cluster <- makeCluster(8)
cluster <- NULL

#*----------------------------------------------------------------*
# Find name of sample and input from their resp. paths

mySample	<- tail(unlist(strsplit(path2sample, split="/")),2)[1]
myInput 	<- tail(unlist(strsplit(path2input, split="/")),2)[1]

#*----------------------------------------------------------------*
# Define paths 
path2peakcalling	<- paste( path2expFolder, "/peakcalling/", sep="" )
path2currSampleFolder	<- paste( path2peakcalling, "/", mySample, "/", sep="" )
path2currInputFolder	<- paste( path2peakcalling, "/", myInput, "/", sep="" )

#*----------------------------------------------------------------*
input.data	<- read.bam.tags(path2sample)
chip.data	<- read.bam.tags(path2input)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#           SPP estimate fragment size for all inputs             *
#*                                                                *
#******************************************************************
  
cat("\n\t Estimating fragment size for ", myInput, "\n", sep="")
  
#-----------------------------------------------------------
# Get binding info from cross-correlation profile
#	srange gives the possible range for the size of the protected region;
#	srange should be higher than tag length; making the upper boundary too high will increase calculation time
#	bin - bin tags within the specified number of basepairs to speed up calculation;
#	increasing bin size decreases the accuracy of the determined parameters  
cluster <- NULL
binding.characteristics <- get.binding.characteristics(input.data, srange=c(50,500), bin=5, cluster=cluster)

# Print binding peak separation distance
cat("\n ---------------------------------------------\n Binding peak seperation distance: \t", binding.characteristics$peak$x, "\n", sep="")

#-----------------------------------------------------------
# Plot cross-correlation profile and save as .pdf
pdffile <- paste( path2currInputFolder, "/" , myInput, ".crosscorrelation.pdf",sep="" )
pdf( file=pdffile, width=5, height=5)
par( mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8);
plot( binding.characteristics$cross.correlation , type='l', xlab="strand shift" , ylab="cross-correlation")
abline( v=binding.characteristics$peak$x, lty=2, col=2)
dev.off()


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#******************************************************************
#*                                                                *
#*                    SPP run PeakCalling                         *
#*                                                                *
#******************************************************************

# Samples
#-----------------------------------------------------------  
cat("\n\n Peak calling for sample: \t", mySample, "\t over \t", myInput, "\n", sep="")

#-----------------------------------------------------------

#*----------------------------------------------------------------*
# Get binding info from cross-correlation
cluster <- NULL  
binding.characteristics <- get.binding.characteristics( chip.data, srange=c(50,500) , bin=5 , cluster=cluster)

#*----------------------------------------------------------------*
# print out binding peak separation distance
cat("\t Binding peak seperation distance: \t", binding.characteristics$peak$x , "\n" , sep="")

#*----------------------------------------------------------------*
# plot cross-correlation profile to a .pdf file
pdffile <- paste( path2currSampleFolder , "/" , mySample, ".crosscorrelation.pdf",sep="" )
pdf(file=pdffile, width=5, height=5)
par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8);
plot(binding.characteristics$cross.correlation,type='l', xlab="strand shift", ylab="cross-correlation");
abline(v=binding.characteristics$peak$x, lty=2, col=2)
dev.off();

#*----------------------------------------------------------------*
# select informative tags based on the binding characteristics
chip.data	<- select.informative.tags(chip.data, binding.characteristics);
input.data	<- select.informative.tags(input.data, binding.characteristics);

#*----------------------------------------------------------------*
# Restrict or remove singular positions with very high tag counts
chip.data	<- remove.local.tag.anomalies(chip.data);
input.data	<- remove.local.tag.anomalies(input.data);

#*----------------------------------------------------------------*
# Compute genome wide density and tag enrichment profile
tag.shift	<- round(binding.characteristics$peak$x/2)

if( makeunique == 1 ){
  
	# Opt1: Output smoothed tag density (subtracting re-scaled input) into a WIG file (for single files)
 	#	Note that the tags are shifted by half of the peak separation distance
	#	Use scale.by.dataset.size=T option to normalize the tag density by the total dataset size (to make it comparable across samples sequenced at different sequencing depth)
	#	Note this option is typically used without the background subtraction
	smoothed.density	<- get.smoothed.tag.density( chip.data , control.tags=input.data , bandwidth=200 , step=100 , tag.shift=tag.shift)
	wigfile			<- paste( path2currSampleFolder, mySample, ".density.wig",sep="" )
	writewig(smoothed.density, wigfile, "Example smoothed, background-subtracted tag density")
	rm(smoothed.density)
  
} else if( makeunique == 2 ){
  
	# Opt2: Smoothed maxlikelihood log2 enrichment (chip over input)
	smoothed.enrichment.estimate	<- get.smoothed.enrichment.mle(chip.data, input.data, bandwidth=200, step=100, tag.shift=tag.shift)
	wigfile				<- paste( path2currSampleFolder, mySample, ".enrichment.wig", sep="" )
	writewig(smoothed.enrichment.estimate, wigfile, "Example smoothed maximum likelihood log2 enrichment estimate")
  
} else if( makeunique == 3 ){
  
	# Opt3: Output conservative enrichment estimates
	#	alpha specifies significance level at which confidence intervals will be estimated
	#	default = fws=500, step=100, alpha=0.01  
	enrichment.estimates	<- get.conservative.fold.enrichment.profile(chip.data, input.data, fws=1000, step=50, alpha=cutoff)
	wigfile			<- paste( path2currSampleFolder, mySample, ".enrichment.estimates.wig", sep="" )
	writewig(enrichment.estimates, wigfile, "Example conservative fold-enrichment/depletion estimates shown on log2 scale")
	rm(enrichment.estimates)
}

broad.clusters <- get.broad.enrichment.clusters(chip.data , input.data , window.size=1e3 , z.thr=3 , tag.shift=round(binding.characteristics$peak$x/2))


#*----------------------------------------------------------------*
# Save in broadPeak format
write.broadpeak.info( broad.clusters , paste(path2currSampleFolder , mySample , ".broadPeak" , sep=""))

#*----------------------------------------------------------------*
# Compute point binding positions
#       using "wtd" method:
#       binding detection parameters
#       the binding.characteristics contains the optimized half-size for binding detection window
detection.window.halfsize <- binding.characteristics$whs

if( removepcr == 1 ){
  
	# Determine binding positions using wtd method
	bp	<- find.binding.positions(signal.data=chip.data, control.data=input.data, fdr=fdr, whs=detection.window.halfsize, cluster=cluster)
	print(paste("detected", sum(unlist(lapply(bp$npl,function(d) length(d$x)))),"peaks"))
  
} else if( removepcr == 2 ){
  
	# By default, the methods exclude large (>1e4bp) regions that exhibit significantly higher input tag density than expected. To prevent such exclusion use tec.filter=F option
	bp	<- find.binding.positions(signal.data=chip.data, control.data=input.data, fdr=fdr, method=tag.lwcc, whs=detection.window.halfsize, cluster=cluster)
  
}

bp <- add.broad.peak.regions(chip.data, input.data, bp, window.size=1000, z.thr=3)

#*----------------------------------------------------------------*
# Save in peaks

# Save as .narrowPeak
write.narrowpeak.binding(bp, paste( path2currSampleFolder, mySample, ".narrowPeak", sep="" ))

# Output detected binding positions
output.binding.results(bp, paste( path2currSampleFolder, mySample, ".binding.positions.txt", sep="" ))



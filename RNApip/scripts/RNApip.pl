#!/usr/bin/perl -w

#my $phred="phred64";

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#************************************************************************
#*                                                                	*
#*                RNApip PERL script					*
#*                                                                	*
#************************************************************************

#************************************************************************
#*                                                                	*
#*      Open Targets.txt and start pipeline 				*
#*                                                                	*
#*----------------------------------------------------------------------*

if( $ARGV[0] ) { $path2expFolder = $ARGV[0]; }
else{ die "\n\n----------------------------------------\n\n Provide the path where to your project: </PATH/TO/PROJECT> \n\n--------------------------------------------------------------------------------\n\n"; }

#*----------------------------------------------------------------------*
# Read Targets.txt file

my $Targets = "$path2expFolder/DataStructure/Targets.txt";
open(INPUT, $Targets) || die "Error opening $Targets : $!\n\n\n";

my ($expFolder, $genome, $userFolder, $path2RNAseqScripts, $path2RNAseq, $path2fastqgz, $chrlens, $path2gtfFile)	 			= ("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA");
my ($unzip, $qc, $map, $filter)															= ("FALSE", "FALSE", "FALSE", "FALSE", "FALSE");
my (@sc, @lines2remove)																= ();
# Find paths to different folders in the Targets.txt file
while(<INPUT>) {
	if (/# My_email/) {
                $_ =~ m/"(.+?)"/;
                $email = "$1";	
	}
	if (/# My_project_title/) {
		$_ =~ m/"(.+?)"/;
		$expFolder = "$1";
	}
	if (/# Reference_genome/) {
		$_ =~ m/"(.+?)"/;
		$genome = "$1";
	}
	if (/# Path_to_proj_folder/) {	
		$_ =~ m/"(.+?)"/;
		$userFolder = "$1";
	}
	if (/# Path_to_RNApip/) {
		$_ =~ m/"(.+?)"/;
		$path2RNAseq = "$1";
		$path2RNAseqScripts = join("", $path2RNAseq, "/scripts");
	}
	if (/# Path_to_orifastq.gz/) {
		$_ =~ m/"(.+?)"/;
		$path2fastqgz = "$1";
	}
	if (/# Path_to_chrLens.dat/) {
                $_ =~ m/"(.+?)"/;
                $chrlens = "$1";
        }
	if (/# Path_to_RefGen.fa/) {
                $_ =~ m/"(.+?)"/;
                $refGenome = "$1";
        }
	if (/# Path_to_gtfFile/) {
                $_ =~ m/"(.+?)"/;
                $path2gtfFile = "$1";
        }
	if (/# Aligner_algorithm/) {
                $_ =~ m/"(.+?)"/;
                $aligner = "$1";
        }
	if (/# Paired_end_run/) {
                $_ =~ m/"(.+?)"/;
                $PE = "$1";
        }
	if (/# Steps_to_execute/) {
		$_ =~ m/"(.+?)"/;
        	@steps2execute = ();
		if (grep /\bunzip\b/i, $_ )		{ $unzip 		= "TRUE"; push @steps2execute, "Unzip";		}
		if (grep /\bqc\b/i, $_ )		{ $qc			= "TRUE"; push @steps2execute, "QC";		}
		if (grep /\bmap\b/i, $_ )		{ $map	 		= "TRUE"; push @steps2execute, "Map";		}
		if (grep /\bfilter\b/i, $_ )            { $filter               = "TRUE"; push @steps2execute, "Filter";    	}
        }

} # end of Targets.txt



my $AdvSettings = "$path2expFolder/DataStructure/AdvancedSettings.txt";
open(INPUT, $AdvSettings) || die "Error opening $AdvSettings : $!\n\n\n";

my ($removepcrdup, $makeunique, $ndiff, $aligncommand)				= ("NA", "NA", "NA", "NA", "NA", "NA", "NA");

while(<INPUT>) {

	if (/# Filter.removePCRdup/) {
		$_ =~ m/"(.+?)"/;
		$removepcr = "$1";
	}
	if (/# Filter.makeUniqueRead/) {
		$_ =~ m/"(.+?)"/;
		$makeunique = "$1";
	}
	if (/# Filter.maxEditDist/) {
		$_ =~ m/"(.+?)"/;
		$ndiff = "$1";
	}
	if (/# Align.command.opt/) {
                $_ =~ m/"(.+?)"/;
                $aligncommand = "$1";
        }

} # end of AdvancedSettings.txt




#*----------------------------------------------------------------------*
# Define paths

my $path2expFolder = "$userFolder/$expFolder";
$Targets = "$path2expFolder/DataStructure/Targets.txt";

#*----------------------------------------------------------------------*

chdir "$path2expFolder";

print "\n##################################################################################################";
print "\n# ";
print "\n#	The pipeline will run the following tasks:\t\t";
print join("  -  ", @steps2execute);
print "\n# ";
print "\n##################################################################################################\n\n";
print "\n";
print "\n My email:\t\t $email";
print "\n";
print "\n expFolder:\t\t $expFolder";
print "\n genome:\t\t $genome";
print "\n userFolder:\t\t $userFolder";
print "\n path2RNApip:\t\t $path2RNAseq";
print "\n path2expFolder:\t $path2expFolder";
print "\n path2fastq.gz:\t\t $path2fastqgz";
print "\n Targets:\t\t $path2expFolder/DataStructure/Targets.txt";
print "\n chrlens:\t\t $chrlens";
print "\n refGenome:\t\t $refGenome";
print "\n";
print "\n Paired end sequencing:\t $PE";
print "\n Aligner algorithm:\t $aligner";
print "\n Align command: \t $aligncommand";
print "\n Remove pcr dupl:\t $removepcr";
print "\n Make unique reads:\t $makeunique";
print "\n";
print "\n Current working dir:\t $path2expFolder";
print "\n";
print "\n .........................................";
print "\n Performing following tasks:";
print "\n .........................................";
print "\n unzip:\t\t\t $unzip";
print "\n qc:\t\t\t $qc";
print "\n map:\t\t\t $map";
print "\n filter:\t\t $filter";
print "\n .........................................";
print "\n";
print "\n Samples: ";
foreach my $i (0 .. $#samples) {
        print "\n\t $samples[$i]";
}
print "\n";
#print "\n";
#print "\n----------------------------------------\n";


#*----------------------------------------------------------------------*
# Parse the Targets.txt file and find unique sample names of FileName and InpName

my @Targets1 = `cut -f1 $Targets`;
	chomp(@Targets1);
my @Targets2 = `cut -f2 $Targets`;
	chomp(@Targets2);

# Store original file names in orisamples

my @orisamples;
foreach $line (@Targets1) {
	$line =~ /^$/ and die "Targets 1: Blank line detected at $.\n\n";
	$line =~ /^[# = " OriFileName FileName OriInpName InpName]/ and next;
	push(@orisamples, $line);
}

# Store original file names in samples
my @samples;
foreach $line (@Targets2) {
	$line =~ /^$/ and die "Targets 1: Blank line detected at $.\n\n";
	$line =~ /^[# = " OriFileName FileName OriInpName InpName]/ and next;
	push(@samples, $line);
}

#*----------------------------------------------------------------------*
# Remove duplicated elements in the list @samples and @inputs
%seen		= ();
@samples	= grep { ! $seen{$_} ++ } @samples;

#*----------------------------------------------------------------------*
# Store variables into @samples
my $cutoff	= 0.0;
#my @chrs	= `cut -f1 $chrlens`;
#chomp(@chrs);

#*----------------------------------------------------------------------*
# Set different paths

my $tmpscr 			= "$path2expFolder/scripts";
my $path2fastq			= "$path2expFolder/fastq";
my $path2QC			= "$path2expFolder/QC";
my $path2HTSeq			= "$path2expFolder/HTSeq";
my $path2Tophat			= "$path2expFolder/Tophat";
my $scrhead 			= "$path2RNAseqScripts/QSUB_header.sh";
my $path2iterate		= "$tmpscr/iterate/";
my $RNAseqMainIterative		= "$path2iterate/RNAseq.sh";
my $IterateSH			= "$path2iterate/IterateSH.sh";


#************************************************************************
#									*
# 			START tasks					*
#								   	*
#*----------------------------------------------------------------------*

#*----------------------------------------------------------------------*
# Unzipping and renaming fastq.gz files

if( $unzip =~ "TRUE" ){

	print "\n\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- \n";
	print "\n Unzipping and renaming files using Targets.txt \n";

	my $iterateJobName	= "Iterate_unzip";
	my $myJobName		= "unzip";
	my $path2qsub		= "$tmpscr/$myJobName/qsub";

	# Create file to store jobs in
	unless( -d "$tmpscr/$myJobName" )	{ `mkdir $tmpscr/$myJobName`; }
	unless( -d "$path2qsub" )		{ `mkdir $path2qsub`; }
	my $QSUB	= "$tmpscr/$myJobName/$myJobName\.sh";
	open $QSUB, ">", "$QSUB" or die "Can't open '$QSUB'";
	print $QSUB "#!/bin/bash\n";
	close $QSUB;
	`chmod 777 $QSUB`;
        print "\n Store all of the following '$myJobName' jobs in $QSUB \n";
        my @myJobs;
	
	print "\n\n-------------------------------------\n\n Unzipping samples: ";
	foreach my $i (0 .. $#samples) {
		
		print "\n\t $orisamples[$i] \t $samples[$i] ";		

		# Prepare a personal qsub script
		my $QSUBint  = "$tmpscr/$myJobName/$samples[$i]\_$myJobName\.sh";
		`cp $scrhead $QSUBint`;

		my $cmd		= "gunzip -c $path2fastqgz/$orisamples[$i]\.fastq\.gz > $path2fastq/$samples[$i]\.fastq";
		`echo "$cmd" >> $QSUBint`;
		
		#---------------------------------------------
		# Keep track of the jobs in @myJobs
		my $jobName	= "Sample_$myJobName$i";
		push(@myJobs, $jobName);
		$cmd		= "$jobName=`qsub -o $path2qsub -e $path2qsub $QSUBint`";
		open $QSUB, ">>", "$QSUB" or die "Can't open '$QSUB'";
		print $QSUB "$cmd\n";
		close $QSUB;     
	}

	#*----------------------------------------------------------------------*
	# Change Targets.txt file for next iteration
	print "\n--------------------------------------------------------------------------------------------------\n";
	print "\n Changing '$myJobName' variable to FALSE and proceed";
	`/usr/bin/perl -p -i -e "s/$myJobName/$myJobName\_DONE/gi" $Targets`;

	#*----------------------------------------------------------------------*
	# Prepar file containing the jobs to run

	# Add the next job line to the $mapQSUB
	foreach( @myJobs ){ $_ = "\$".$_ ; }
	my $myJobsVec	= join(":", @myJobs);
	my $finalcmd    = "FINAL=\`qsub -N $iterateJobName -o $path2qsub -e $path2qsub -W depend=afterok\:$myJobsVec $IterateSH`";

	open $QSUB, ">>", "$QSUB" or die "Can't open '$QSUB'";
	print $QSUB "$finalcmd\n";
	close $QSUB;

	#*----------------------------------------------------------------------*
	# Submit jobs to run

	print "\n\n--------------------------------------------------------------------------------------------------\n";
	print "\n Submitting job to cluster: \t `sh $QSUB` \n";
	`sh $QSUB`;

	#*----------------------------------------------------------------------*
	# Exit script

	print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	print "\n Exiting $myJobName section with no known error \n";
	print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n";

	exit 0;
} 


#*----------------------------------------------------------------------*
# Quality Control

if( $qc =~ "TRUE" ) {

	print "\n\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- \n";
	print "\n QC fastq files \n";

	my $iterateJobName	= "Iterate_QC";
	my $myJobName		= "QC";
	my $path2qsub		= "$tmpscr/$myJobName/qsub";

	# Create file to store jobs in
	unless( -d "$tmpscr/$myJobName" )	{ `mkdir $tmpscr/$myJobName`; }
	unless( -d "$path2qsub" )		{ `mkdir $path2qsub`; }
	my $QSUB        = "$tmpscr/$myJobName/$myJobName\.sh";
	open $QSUB, ">", "$QSUB" or die "Can't open '$QSUB'";
	print $QSUB "#!/bin/bash\n";
	close $QSUB;	
	`chmod 777 $QSUB`;
	print "\n Store all of the following '$myJobName' jobs in $QSUB \n";
	my @myJobs;

	unless( -d "$path2QC" ) { `mkdir $path2QC`; }
	`cp $path2RNAseqScripts/QC.R $tmpscr`;

	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	my $code 	= "$tmpscr/QC.R" ;
	my $cmd		= "Rscript $code $path2expFolder &>> $path2qsub/QCReport.log";
	`echo "$cmd" >> $QSUB`;
	#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--


	#*----------------------------------------------------------------------*
        # Change Targets.txt file for next iteration
        print "\n--------------------------------------------------------------------------------------------------\n";
        print "\n Changing '$myJobName' variable to FALSE and proceed";
        `/usr/bin/perl -p -i -e "s/$myJobName/$myJobName\_DONE/gi" $Targets`;

	#*----------------------------------------------------------------------*
	# Prepar file containing the jobs to run

	# Add the next job iteration
	my $finalcmd    = "FINAL=\`qsub -N $iterateJobName -o $path2qsub -e $path2qsub $IterateSH`";

	open $QSUB, ">>", "$QSUB" or die "Can't open '$QSUB'";
	print $QSUB "$finalcmd\n";
	close $QSUB;

	#*----------------------------------------------------------------------*
	# Submit jobs to run

	print "\n\n--------------------------------------------------------------------------------------------------\n";
	print "\n Submitting job to cluster: \t `sh $QSUB` \n";
	`sh $QSUB`;

	#*----------------------------------------------------------------------*
	# Exit script

	print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	print "\n Exiting $myJobName section with no known error \n";
	print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n";

	exit 0;
}



#*----------------------------------------------------------------------*
# Mapping sequences with tophat

if( $map =~ "TRUE" ){	

	print "\n\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n";
	print "\n Mapping fastq files using Tophat\n";

	my $iterateJobName	= "Iterate_map";
	my $myJobName		= "map";
	my $myJobName2		= "rename";
	my $path2qsub		= "$tmpscr/$myJobName/qsub";

	# Create file to store jobs in
	unless( -d "$tmpscr/$myJobName" )	{ `mkdir $tmpscr/$myJobName`; }
	unless( -d "$path2qsub" )		{ `mkdir $path2qsub`; }
	my $QSUB        = "$tmpscr/$myJobName/$myJobName\.sh";
	open $QSUB, ">", "$QSUB" or die "Can't open '$QSUB'";
	print $QSUB "#!/bin/bash\n";
	close $QSUB;
	`chmod 777 $QSUB`;
	print "\n Store all of the following '$myJobName' jobs in $QSUB \n";
	my @myJobs;
	my @myJobs2;


	#*----------------------------------------------------------------------*
	# Create a folder named * mysample * within each bwa_sam, bwa_saf and bwa_sai folders

	print "\n\n-------------------------------------\n\n Mapping samples: ";
	foreach my $i (0 .. $#samples) {
		
		print "\n\t $samples[$i] ";

		# Prepare a personal qsub script
		my $QSUBint		= "$tmpscr/$myJobName/$samples[$i]\_$myJobName\.sh";
		`cp $scrhead $QSUBint`;
		my $QSUBintRename	= "$tmpscr/$myJobName/$samples[$i]\_$myJobName2\.sh";
		`cp $scrhead $QSUBintRename`;
		
		# Create a folder for each sample to store files
		unless (-d "$path2Tophat/$samples[$i]")		{`mkdir $path2Tophat/$samples[$i]`;}
		
		if( $PE ) {
			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			my $cmd		= "$aligncommand $path2gtfFile -o $path2Tophat/$samples[$i] $refGenome/$genome $path2fastq/$samples[$i]\_R1.fastq $path2fastq/$samples[$i]\_R2.fastq";
			`echo "$cmd" >> $QSUBint`;
			#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--			

		} else {
			
			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			my $cmd		= "$aligncommand $path2gtfFile -o $path2Tophat/$samples[$i] $refGenome/$genome $path2fastq/$samples[$i].fastq";
			`echo "$cmd" >> $QSUBint`;
			#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

		}

		#---------------------------------------------
		# Rename samples
		#
		my $cmd		= "mv $path2Tophat/$samples[$i]/accepted_hits.bam $path2Tophat/$samples[$i]/$samples[$i]\_unfiltered.bam";
		`echo "$cmd" >> $QSUBintRename`;

		#---------------------------------------------
		# Keep track of the jobs in @myJobs
		my $jobName	= "$myJobName$i";
		push(@myJobs, $jobName);
		my $jobName2	= "$myJobName2$i";
		push(@myJobs2, $jobName2);
		
		my $cmd1	= "$jobName=`qsub -o $path2qsub -e $path2qsub $QSUBint`";
		$jobName = "\$".$jobName ;
		my $cmd2	= "$jobName2=`qsub -o $path2qsub -e $path2qsub -W depend=afterok\:$jobName $QSUBintRename`";
		open $QSUB, ">>", "$QSUB" or die "Can't open '$QSUB'";
		print $QSUB "$cmd1\n";
		print $QSUB "$cmd2\n";
		close $QSUB;

	}


	#*----------------------------------------------------------------------*
	# Change Targets.txt file for next iteration
	print "\n--------------------------------------------------------------------------------------------------\n";
	print "\n Changing '$myJobName' variable to FALSE and proceed";
	`/usr/bin/perl -p -i -e "s/$myJobName/$myJobName\_DONE/gi" $Targets`;

	#*----------------------------------------------------------------------*
	# Prepar file containing the jobs to run

	# Add the next job line to the $mapQSUB
	foreach( @myJobs2 ){ $_ = "\$".$_ ; }
	my $myJobsVec	= join(":", @myJobs2);
	my $finalcmd	= "FINAL=\`qsub -N $iterateJobName -o $path2qsub -e $path2qsub -W depend=afterok\:$myJobsVec $IterateSH`";
	open $QSUB, ">>", "$QSUB" or die "Can't open '$QSUB'";
	print $QSUB "$finalcmd\n";
	close $QSUB;

	#*----------------------------------------------------------------------*
	# Submit jobs to run
	
	print "\n\n--------------------------------------------------------------------------------------------------\n";
	print "\n Submitting job to cluster: \t `sh $QSUB` \n";
	`sh $QSUB`;

	#*----------------------------------------------------------------------*
	# Exit script

	print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	print "\n Exiting $myJobName section with no known error \n";
	print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n";

	exit 0;

} 



#*----------------------------------------------------------------*
# Filtering reads

if( $filter =~ "TRUE" ){
	
	my %hChrs		= ();

 	print "\n\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n";
        print "\n Filtering reads\n";

	my $iterateJobName	= "Iterate_filter";
	my $myJobName		= "filter";
	my $path2qsub		= "$tmpscr/$myJobName/qsub";

	# Create file to store jobs in
	unless( -d "$tmpscr/$myJobName" )	{ `mkdir $tmpscr/$myJobName`; }
	unless( -d "$path2qsub" )		{ `mkdir $path2qsub`; }
	my $QSUB	= "$tmpscr/$myJobName/$myJobName\.sh";
	open $QSUB, ">", "$QSUB" or die "Can't open '$QSUB'";
	print $QSUB "#!/bin/bash\n";
	close $QSUB;
	`chmod 777 $QSUB`;
	print "\n Store all of the following '$myJobName' jobs in $QSUB \n";
	my @myJobs;


	print "\n\n-------------------------------------\n\n Filtering samples: ";
	foreach my $i (0 .. $#samples) {

		print "\n\t $samples[$i] ";

		# Prepare a personal qsub script
		my $QSUBint	= "$tmpscr/$myJobName/$samples[$i]\_$myJobName\.sh";
		`cp $scrhead $QSUBint`;

		# -----------------------------------------
                # Remove pcr duplications
		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		if( $removepcr ){
			# bam to sorted bam
			my $cmd = "samtools sort $path2Tophat/$samples[$i]/$samples[$i]\_unfiltered.bam $path2Tophat/$samples[$i]/$samples[$i]\_sortedwpcr";
			`echo "$cmd" >> $QSUBint`;
			$cmd = "samtools rmdup -s $path2Tophat/$samples[$i]/$samples[$i]\_sortedwpcr.bam $path2Tophat/$samples[$i]/$samples[$i].bam";
			`echo "$cmd" >> $QSUBint`;
		} else {
			# bam to sorted bam
			my $cmd = "samtools sort $path2Tophat/$samples[$i]/$samples[$i]\_unfiltered.bam $path2Tophat/$samples[$i]/$samples[$i].bam";
			`echo "$cmd" >> $QSUBint`;
		}

		my $cmd= "samtools index $path2Tophat/$samples[$i]/$samples[$i].bam $path2Tophat/$samples[$i]/$samples[$i].bai";
		`echo "$cmd" >> $QSUBint`;
		#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--


		print "\n\n Samtools processing for $samples[$i] done. \n\n";
		
		#---------------------------------------------
		# Keep track of the jobs in @myJobs
		my $jobName	= "$myJobName$i";
		push(@myJobs, "$jobName");
		$cmd        	= "$jobName=`qsub -o $path2qsub -e $path2qsub $QSUBint`";
		open $QSUB, ">>", "$QSUB" or die "Can't open '$QSUB'";
		print $QSUB "$cmd\n";
		close $QSUB;

	}

	#*----------------------------------------------------------------------*
	# Change Targets.txt file for next iteration
	print "\n--------------------------------------------------------------------------------------------------\n";
	print "\n Changing '$myJobName' variable to FALSE and proceed";
	`/usr/bin/perl -p -i -e "s/$myJobName/$myJobName\_DONE/gi" $Targets`;

	#*----------------------------------------------------------------------*
	# Prepar file containing the jobs to run
	
	# Add the next job line to the $QSUB
	foreach( @myJobs ){ $_ = "\$".$_ ; }
	my $myJobsVec   = join(":", @myJobs);
	my $finalcmd    = "FINAL=\`qsub -N $iterateJobName -o $path2qsub -e $path2qsub -W depend=afterok\:$myJobsVec $IterateSH`";
	open $QSUB, ">>", "$QSUB" or die "Can't open '$QSUB'";
	print $QSUB "$finalcmd\n";
	close $QSUB;

	#*----------------------------------------------------------------------*
	# Submit jobs to run

	print "\n\n--------------------------------------------------------------------------------------------------\n";
	print "\n Submitting job to cluster: \t `sh $QSUB` \n";
	`sh $QSUB`;

	#*----------------------------------------------------------------------*
	# Exit script

	print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	print "\n Exiting $myJobName section with no known error \n";
	print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n";

	exit 0;
}


if($unzip =~ "FALSE"  &&  $qc =~ "FALSE"  &&  $map =~ "FALSE"  &&  $filter =~ "FALSE" ){

	print "\n\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n";
	print "\n Exiting \n";

	my $iterateJobName	= "Iterate_exit";
	my $myJobName           = "exit";
	my $path2qsub           = "$tmpscr/$myJobName/qsub";

	# Create file to store jobs in
 	unless( -d "$tmpscr/$myJobName" )	{ `mkdir $tmpscr/$myJobName`; }
 	unless( -d "$path2qsub" )		{ `mkdir $path2qsub`; }

	my $QSUBfinal	= "$tmpscr/$myJobName/$myJobName\_final.sh";
	`cp $scrhead $QSUBfinal`;
	my $cmd		= "echo `An email has been sent to $email`";
	`echo "$cmd" >> $QSUBfinal`;
	`chmod 777 $QSUBfinal`;	

	my $QSUB	= "$tmpscr/$myJobName/$myJobName\.sh";
	$cmd		= "#!/bin/bash";
	`echo "$cmd" >> $QSUB`;
	$cmd		= "FINAL=\`qsub -m e -M $email -N $iterateJobName -o $path2qsub -e $path2qsub $QSUBfinal`";
	`echo "$cmd" >> $QSUB`;
	`chmod 777 $QSUB`;


        #*----------------------------------------------------------------------*
        # Submit jobs to run

        print "\n\n--------------------------------------------------------------------------------------------------\n";
        print "\n Submitting job to cluster: \t `sh $QSUB` \n";
        `sh $QSUB`;

	print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	print "\n Exiting $myJobName section with no known error \n";
	print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n";

	exit 0;

}


exit 0;

print "\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n\n";

#*----------------------------------------------------------------------*



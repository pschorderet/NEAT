#!/usr/bin/perl -w

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#************************************************************************
#*                                                                	*
#*                ChIPseq PERL script					*
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

my ($expFolder, $genome, $userFolder, $path2ChIPseqScripts, $path2ChIPseq, $path2fastqgz) 	= ("NA", "NA", "NA", "NA", "NA", "NA");
my ($unzip, $qc, $map, $filter, $peakcalling, $cleanbigwig, $cleanfolders)			= ("FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE");
my (@sc, @lines2remove)										= ();
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
	if (/# Path_to_ChIPpip/) {
		$_ =~ m/"(.+?)"/;
		$path2ChIPseq = "$1";
		$path2ChIPseqScripts = join("", $path2ChIPseq, "/scripts");
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
	if (/# Paired_end_run/) {
                $_ =~ m/"(.+?)"/;
                $PE = "$1";
        }
	if (/# Aligner_algorithm/) {
		$_ =~ m/"(.+?)"/;
		$aligner = "$1";
	}
	if (/# Steps_to_execute/) {
		$_ =~ m/"(.+?)"/;
        	@steps2execute = ();
		if (grep /\bunzip\b/i, $_ )		{ $unzip 		= "TRUE"; push @steps2execute, "Unzip";		}
		if (grep /\bqc\b/i, $_ )		{ $qc			= "TRUE"; push @steps2execute, "QC";		}
		if (grep /\bmap\b/i, $_ )		{ $map	 		= "TRUE"; push @steps2execute, "Map";		}
		if (grep /\bfilter\b/i, $_ )		{ $filter 		= "TRUE"; push @steps2execute, "Filter";	}
		if (grep /\bpeakcalling\b/i, $_ )	{ $peakcalling		= "TRUE"; push @steps2execute, "Peakcalling";	}
		if (grep /\bcleanbigwig\b/i, $_ )	{ $cleanbigwig		= "TRUE"; push @steps2execute, "Cleanbigwig";	}
        }
	if (/# Remove_from_bigwig/) {
		$_ =~ m/"(.+?)"/;
		my $text = "$1";
		my @var = split(",", $text);
		foreach my $line (@var) {
			$line =~ s/\s+//g;
			push(@lines2remove, $line);
		}
	}

} # end of Targets.txt



my $AdvSettings = "$path2expFolder/DataStructure/AdvancedSettings.txt";
open(INPUT, $AdvSettings) || die "Error opening $AdvSettings : $!\n\n\n";

my ($removepcrdup, $makeunique, $ndiff, $aligncommand1, $aligncommand2, $fdr, $posopt, $densityopt, $enforceisize)		= ("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA");

while(<INPUT>) {

	if (/# Bwa.maxEditDist/) {
		$_ =~ m/"(.+?)"/;
		$ndiff = "$1";
	}
	if (/# Align.command.line.1/) {
		$_ =~ m/"(.+?)"/;
		$aligncommand1 = "$1";
	}
	if (/# Align.command.line.2/) {
		$_ =~ m/"(.+?)"/;
		$aligncommand2 = "$1";
        }
 	if (/# Filter.removePCRdup/) {
		$_ =~ m/"(.+?)"/;
		$removepcr = "$1";
	}
	if (/# Filter.makeUniqueRead/) {
		$_ =~ m/"(.+?)"/;
		$makeunique = "$1";
 	}
	if (/# Filter.splitbychr/) {
		$_ =~ m/"(.+?)"/;
		$splitbychr = "$1";
	}
	if (/# Filter.enforceinssize/) {
		$_ =~ m/"(.+?)"/;
		$enforceisize = "$1";
	}
	if (/# Filter.minisize/) {
		$_ =~ m/"(.+?)"/;
		$minisize = "$1";
	}
	if (/# Filter.maxisize/) {
		$_ =~ m/"(.+?)"/;
		$maxisize = "$1";
	}
	if (/# PeakCaller.fdr/) {
		$_ =~ m/"(.+?)"/;
		$fdr = "$1";
	}
	if (/# PeakCaller.posopt/) {
		$_ =~ m/"(.+?)"/;
		$posopt = "$1";
	}
	if (/# PeakCaller.densityopt/) {
		$_ =~ m/"(.+?)"/;
		$densityopt = "$1";
	}

} # end of AdvancedSettings.txt


#*----------------------------------------------------------------------*
# Parse the Targets.txt file and find unique sample names of FileName and InpName

my @Targets1 = `cut -f1 $Targets`;
        chomp(@Targets1);
my @Targets2 = `cut -f2 $Targets`;
        chomp(@Targets2);
my @Targets3 = `cut -f3 $Targets`;
        chomp(@Targets3);
my @Targets4 = `cut -f4 $Targets`;
        chomp(@Targets4);

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
my @oriinputs;
foreach $line (@Targets3) {
        $line =~ /^$/ and die "Targets 3: Blank line detected at $.\n\n";
        $line =~ /^[# = " OriFileName FileName OriInpName InpName]/ and next;
        push(@oriinputs, $line);
}
my @inputs;
foreach $line (@Targets4) {
        $line =~ /^$/ and next;
        $line =~ /^[# = " OriFileName FileName OriInpName InpName]/ and next;
        push(@inputs, $line);
}


#*----------------------------------------------------------------------*
# Remove duplicated elements in the list @samples and @inputs
%seen           = ();
@samples        = grep { ! $seen{$_} ++ } @samples;
@allinputs	= @inputs;
%seen           = ();
@inputs         = grep { ! $seen{$_} ++ } @inputs;
@samplesInputs  = @samples;
push (@samplesInputs, @inputs);

#*----------------------------------------------------------------------*
# Store variables into @samples
my $cutoff	= 0.0;
my @chrs        = `cut -f1 $chrlens`;
chomp(@chrs);



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
print "\n path2ChIPpip:\t\t $path2ChIPseq";
print "\n path2expFolder:\t $path2expFolder";
print "\n path2fastq.gz:\t\t $path2fastqgz";
print "\n Targets:\t\t $path2expFolder/DataStructure/Targets.txt";
print "\n chrlens:\t\t $chrlens";
print "\n refGenome:\t\t $refGenome";
print "\n";
print "\n Paired end sequencing:\t $PE";
print "\n Aligner algorithm:\t $aligner";
print "\n Align command line 1:\t $aligncommand1";
print "\n Align command line 2:\t $aligncommand2";
print "\n";
print "\n Remove pcr dupl:\t $removepcr";
print "\n Make unique reads:\t $makeunique";
print "\n PeakCaller.fdr:\t $fdr";
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
print "\n peakcalling:\t\t $peakcalling";
print "\n cleanbigwig:\t\t $cleanbigwig \t (remove: @lines2remove)";
print "\n .........................................";
print "\n Samples: ";
foreach my $i (0 .. $#samples) {
        print "\n\t $samples[$i] \t - \t $inputs[$i]";
}
print "\n";
#print "\n";
#print "\n----------------------------------------\n";


#*----------------------------------------------------------------------*
# Set different paths

my $tmpscr 			= "$path2expFolder/scripts";
my $path2fastq			= "$path2expFolder/fastq";
my $path2QC			= "$path2expFolder/QC";
my $path2aligned		= "$path2expFolder/aligned";
my $path2peakcalling            = "$path2expFolder/peakcalling";
my $scrhead 			= "$path2ChIPseqScripts/QSUB_header.sh";
my $path2iterate		= "$tmpscr/iterate/";
my $ChIPseqMainIterative	= "$path2iterate/ChIPseq.sh";
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

	foreach my $i (0 .. $#samples) {

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

	foreach my $i (0 .. $#inputs) {

		# Prepare a personal qsub script
		my $QSUBint  = "$tmpscr/$myJobName/$inputs[$i]\_$myJobName\.sh";
		`cp $scrhead $QSUBint`;

		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		my $cmd         = "gunzip -c $path2fastqgz/$oriinputs[$i]\.fastq\.gz > $path2fastq/$inputs[$i]\.fastq";
		`echo "$cmd" >> $QSUBint`;
		#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

		#---------------------------------------------
		# Keep track of the jobs in @myJobs
		my $jobName     = "Input_$myJobName$i";
		push(@myJobs, $jobName);
		$cmd            = "$jobName=`qsub -o $path2qsub -e $path2qsub $QSUBint`";
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
	`cp $path2ChIPseqScripts/QC.R $tmpscr`;

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
# Mapping sequences with bwa

if( $map =~ "TRUE" ){	

	print "\n\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n";
	print "\n Mapping fastq files\n";

	my $iterateJobName	= "Iterate_map";
	my $myJobName		= "map";
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


	foreach my $i (0 .. $#samplesInputs) {
		
		# Prepare a personal qsub script
		my $QSUBint	= "$tmpscr/$myJobName/$samplesInputs[$i]\_$myJobName\.sh";
		`cp $scrhead $QSUBint`;
		
		# Create a directory
		my $path2currentSampleDir	= "$path2aligned/$samplesInputs[$i]";
		unless( -d "$path2currentSampleDir" )	{ `mkdir $path2currentSampleDir`; }


		if( $PE ) {
			
			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			my $cmd		= "$aligncommand1 -n $ndiff $refGenome $path2fastq/$samplesInputs[$i]\_1.fastq > $path2currentSampleDir/$samplesInputs[$i]\_1.sai";
			`echo "$cmd" >> $QSUBint`;
			$cmd		= "$aligncommand1 -n $ndiff $refGenome $path2fastq/$samplesInputs[$i]\_2.fastq > $path2currentSampleDir/$samplesInputs[$i]\_2.sai";
			`echo "$cmd" >> $QSUBint`;
			$cmd		= "$aligncommand2 $refGenome $path2currentSampleDir/$samplesInputs[$i]\_1.sai $path2currentSampleDir/$samplesInputs[$i]\_2.sai $path2fastq/$samplesInputs\_1.fastq $path2fastq/$samplesInputs[$i]\_2.fastq > $path2currentSampleDir/$samplesInputs[$i]\.sam";
			`echo "$cmd" >> $QSUBint`;
			#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

		} else {

			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			my $cmd		= "$aligncommand1 -n $ndiff $refGenome $path2fastq/$samplesInputs[$i]\.fastq > $path2currentSampleDir/$samplesInputs[$i]\.sai";
			`echo "$cmd" >> $QSUBint`;
			my $cmd2	= "$aligncommand2 $refGenome $path2currentSampleDir/$samplesInputs[$i]\.sai $path2fastq/$samplesInputs[$i]\.fastq > $path2currentSampleDir/$samplesInputs[$i]\.sam";
			`echo "$cmd2" >> $QSUBint`;
			#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

		}

		#---------------------------------------------
		# Keep track of the jobs in @myJobs
		my $jobName     = "$myJobName$i";
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

	foreach my $i (0 .. $#samplesInputs) {

		# Prepare a personal qsub script
		my $QSUBint	= "$tmpscr/$myJobName/$samplesInputs[$i]\_$myJobName\.sh";
		`cp $scrhead $QSUBint`;
				
		my $path2currentSampleDir	= "$path2aligned/$samplesInputs[$i]";
		my $j=0;
		
		# -----------------------------------------		
		# Get unique matches only
		my $samplep = $samplesInputs[$i];
		if( $makeunique && ($enforceisize==0) ){

			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			# print "\n\n Getting uniquely mapped reads for $path2aligned/$samplesInputs[$i]\.sam";
			my $cmd		= "grep -E '\\sX0:i:1\\s' $path2currentSampleDir/$samplesInputs[$i]\.sam > $path2currentSampleDir/$samplesInputs[$i]\.u.sam";
                        `echo "$cmd" >> $QSUBint`;
			#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

		}			


		#------------------------------------------------------------   Start Ayla's stuff    ---------------------------------------------------------------
		# Not sure exactly what this is...
		
		# enforce insert size rules 
		if( $enforceisize ) {
			print "\n Enforcing size rules for $path2currentSampleDir/$samplesInputs[$i]\.sam \n";
			open( F1,"< $path2currentSampleDir/$samplesInputs[$i]\.sam" );
			open( F2,"> $path2currentSampleDir/$samplesInputs[$i]\.i\.sam" );
			my $lastname = "";
			my $lastline = "";
			my %seen=();
			while( my $line=<F1> ) {
				if( ($makeunique and $line=~/\sX0:i:1\s/) or $makeunique==0 ){
					my @a = split( /\t/,$line );
					my $isize = abs($a[8]);
			
					if( $isize > $maxisize or $isize < $minisize ){ next; }
					if( ($a[8]>0 and $a[3]>$a[7]) or ($a[8]<0 and $a[3]<$a[7]) ){ next; } 
					# I think I am trying to throw away non-paired reads. 
					# Because I dumped non-uniquely mapping reads, there could be unpaired guys here...
					# Use a hash to store, hope the pair is close enough that hash doens't grow to big.
					if( exists $seen{$a[0]} ) {
						print F2 "$seen{$a[0]}";
						print F2 "$line";
						delete $seen{$a[0]};
					} else {
						$seen{$a[0]} = $line;
					}
			
				}
			}
			$prefixp = "$prefixp\.i";
		}

		# remove pcr stacks
		if( 0 and $removepcr ){

			# first sort to remove stacks
			my $cmd = "/usr/local/bin/IGVTools/igvtools sort $path2currentSampleDir/$samplesInputs[$i]\.sam $path2currentSampleDir/$samplesInputs[$i]\.sorted.sam";
			print "\n $cmd \n";
			`$cmd`;

			open( F1,"< $path2currentSampleDir/$samplesInputs[$i]\.sorted.sam" );
			open( F2,"> $path2currentSampleDir/$samplesInputs[$i]\.p.sam" );
			my $lastn = 'chrname';
			my $lastp = -1;
			my $lastr = 'readseq';
			while( my $line=<F1> ) {
				if( $line=~/^@/ ){ print F2 $line; }
				else {
					my @a= split( /\t/,$line );
					if( $a[2]=~/^$lastn$/ and $a[3]==$lastp and $a[9]=~/^$lastr$/ ){ next; }
					else {
						$hChrs{$a[2]}=1;
						print F2 $line;
						$lastn = $a[2];
						$lastp = $a[3];
						$lastr = $a[9];
					}
				}
			}
			close F1;
			close F2;
			$prefixp = "$prefixp\.p";
		}
		#------------------------------------------------------------	End Ayla's stuff    ---------------------------------------------------------------

		
		# .sam to .bam		
		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		# sam to bam
		$cmd		= "samtools view -b $path2currentSampleDir/$samplesInputs[$i]\.u.sam -T $refGenome -o $path2currentSampleDir/$samplesInputs[$i]\.u.unsorted\.bam";
		`echo "$cmd" >> $QSUBint`;
		#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--
	
		# -----------------------------------------
		if( $removepcr ) {
			
			# .bam to sorted .bam (with removed pcr dup)
			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			my $cmd		= "samtools sort  $path2currentSampleDir/$samplesInputs[$i]\.u.unsorted\.bam $path2currentSampleDir/$samplesInputs[$i]\.u.sortedwpcr";
			`echo "$cmd" >> $QSUBint`;
			$cmd		= "samtools rmdup -s $path2currentSampleDir/$samplesInputs[$i]\.u.sortedwpcr\.bam $path2currentSampleDir/$samplesInputs[$i]\.bam";
			`echo "$cmd" >> $QSUBint`;
			#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

		} else {

			# .bam to sorted .bam
			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			# bam to sorted bam
			$cmd		= "samtools sort  $path2currentSampleDir/$samplesInputs[$i]\.u.unsorted\.bam $path2currentSampleDir/$samplesInputs[$i]\.bam";
			`echo "$cmd" >> $QSUBint`;
			#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--	
	
		}

		# -----------------------------------------
		# Sorted .bam to .bai
		my $cmdf		= "samtools index $path2currentSampleDir/$samplesInputs[$i]\.bam $path2currentSampleDir/$samplesInputs[$i]\.bai";
		`echo "$cmdf" >> $QSUBint`;

		# -----------------------------------------
		# Split indexed .bam by chr
		if( $splitbychr ) {
			unless( -d "$path2currentSampleDir/splitbychr" )	{ `mkdir $path2currentSampleDir/splitbychr`; }
			foreach my $chr( keys %hChrs ) {
				print "\n\n Entered split by chromosome \n\n";
				#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
				#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+	
				# .bam to sorted .bam by chr
				$cmd		= "samtools view -b $path2currentSampleDir/$samplesInputs[$i]\.bam $chr > $path2currentSampleDir/splitbychr/$samplesInputs[$i]\.$chr\.bam";
				`echo "$cmd" >> $QSUBint`;
				$cmd		= "samtools index $path2currentSampleDir/splitbychr/$samplesInputs[$i]\.$chr\.bam $path2currentSampleDir/splitbychr/$samplesInputs[$i]\.$chr\.bai";
				`echo "$cmd" >> $QSUBint`;
				#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

			}
		}

		print "\n\n Samtools processing for $samplesInputs[$i] done. \n\n";
		
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


#*----------------------------------------------------------------------*
# Running PeakCaller

if( $peakcalling =~ "TRUE" ){

	print "\n\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n";
	print "\n Run PeakCaller (SPP)\n";

	my $iterateJobName	= "Iterate_peakcalling";
	my $myJobName		= "peakcalling";
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
	my @myJobsInputs;
	my @myJobsSamples;

	# Copy script and rreate folder
	`cp $path2ChIPseqScripts/PeakCalling.R $tmpscr`;

	foreach my $i (0 .. $#samples) {

		my $path2currentSampleDir	= "$path2aligned/$samples[$i]";
		my $path2currentInputDir	= "$path2aligned/$allinputs[$i]";	

		unless( -d "$path2peakcalling/$samples[$i]" )	{ `mkdir $path2peakcalling/$samples[$i]`; }
		unless( -d "$path2peakcalling/$allinputs[$i]" )	{ `mkdir $path2peakcalling/$allinputs[$i]`; }		

		#-----------------------------------------------------------
		# Prepare a personal qsub script
		my $QSUBint  = "$tmpscr/$myJobName/$samples[$i]\_$myJobName\.sh";
		`cp $scrhead $QSUBint`;
		
		#-----------------------------------------------------------
                # Parameters
		my $path2currentSample	= "$path2currentSampleDir/$samples[$i]\.bam";
		my $path2currentInput	= "$path2currentInputDir/$allinputs[$i]\.bam";
		my $code		= "$tmpscr/PeakCalling.R";

		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+		
		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		my $cmd		= "Rscript $code $path2expFolder $path2currentSample $path2currentInput $fdr $posopt $densityopt &>> $path2qsub/$samples[$i]\_peakcalling.log";
		`echo "$cmd" >> $QSUBint`;
                #--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

		#---------------------------------------------
		# Keep track of the jobs in @myJobs
		my $jobName	= "$myJobName$i";
		push(@myJobsSamples, "$jobName");
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

        # Add the next job line to the $filterQSUB
        foreach( @myJobsSamples ){ $_ = "\$".$_ ; }
        my $myJobsVec   = join(":", @myJobsSamples);
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
# Clean up bigwig files

if( $cleanbigwig =~ "TRUE" ){

	print "\n\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n";
	print "\n Bigwig file cleaning\n";

	my $iterateJobName	= "Iterate_cleanbigwig";
	my $myJobName		= "cleanbigwig";
	my $myJobName2		= "wigToBigwig";
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


	my $QSUBintWig  = "$tmpscr/$myJobName/$myJobName2\.sh";
	open $QSUBintWig, ">", "$QSUBintWig" or die "Can't open '$QSUBintWig'";
	print $QSUBintWig "#!/bin/bash\n";
	close $QSUBintWig;
	`chmod 777 $QSUBintWig`;
	
	#*----------------------------------------------------------------------*
        # Open special character file and store in @specialchar
	
	print "\n These chromosomes will be deleted: @sc \n";
	
	# Create bigwig folder
	print "\n Create folders\n";
	print "\n\t bigwig folder: \t $path2expFolder/bigwig/";    
	unless( -d "$path2peakcalling/bigwig" )		{ `mkdir $path2peakcalling/bigwig`; }

	# Create narrowPeak folder
	print "\n\t narrowPeak folder: \t $path2peakcalling/narrowPeak/";    
	unless( -d "$path2peakcalling/narrowPeak" )	{ `mkdir $path2peakcalling/narrowPeak`; }

	# Create broadPeak folder
	print "\n\t broadPeak folder: \t $path2peakcalling/broadPeak/ \n";    
	unless( -d "$path2peakcalling/broadPeak" )	{ `mkdir $path2peakcalling/broadPeak`; }
    
	foreach my $i (0 .. $#samples) {

		# Prepare a personal qsub script
		my $QSUBint 	= "$tmpscr/$myJobName/$samples[$i]\_$myJobName\.sh";
		`cp $scrhead $QSUBint`;		

		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		# Move narrow and broadPeaks to their resp. folder
		my $cmd         = "mv $path2peakcalling/$samples[$i]/$samples[$i]\.broadPeak $path2peakcalling/broadPeak/";
		`echo "$cmd" >> $QSUBint`;
		$cmd		= "mv $path2peakcalling/$samples[$i]/$samples[$i]\.narrowPeak $path2peakcalling/narrowPeak/";
		`echo "$cmd" >> $QSUBint`;

		# Remove all lines starting with the character in @line2remove
		# Rename the original .wig file
		$cmd		= "cp $path2peakcalling/$samples[$i]/$samples[$i]\.density.wig $path2peakcalling/$samples[$i]/$samples[$i]\.density0.wig";
		`echo "$cmd"  >> $QSUBint`;
		foreach my $j (0 .. $#lines2remove) {
			# Remove lines from bigwig that do not correspond to specified chromosomes
			my $k = $j+1;
			my $cmd		= "sed '/$lines2remove[$j]/d' $path2peakcalling/$samples[$i]/$samples[$i]\.density$j.wig > $path2peakcalling/$samples[$i]/$samples[$i]\.density$k.wig";
			`echo "$cmd"  >> $QSUBint`;
			# Delete intermediary files
			$cmd		= "rm $path2peakcalling/$samples[$i]/$samples[$i]\.density$j.wig";
			`echo "$cmd"  >> $QSUBint`;
		}

		my $finalcount	= $#lines2remove+1;
		$cmd		= "$QSUBintWig";
		`echo "$cmd"  >> $QSUBint`;

		$cmd		= "wigToBigWig $path2peakcalling/$samples[$i]/$samples[$i]\.density$finalcount.wig $chrlens $path2peakcalling/bigwig/$samples[$i]\.bw";
		`echo "$cmd"  >> $QSUBintWig`;

		#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

		# Keep track of the jobs in @myJobs
		#---------------------------------------------
		my $jobName	= "$myJobName$i";
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

	# Add the next job line to the $filterQSUB
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

	exit 1;

}



if($unzip =~ "FALSE"  &&  $qc =~ "FALSE"  &&  $map =~ "FALSE"  &&  $filter =~ "FALSE"  &&  $peakcalling =~ "FALSE"  &&  $cleanbigwig =~ "FALSE"){

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
	my $cmd		= "An email has been sent to $email";
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



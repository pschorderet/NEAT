#!/usr/bin/perl -w

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

my ($expFolder, $genome, $userFolder, $path2RNAseqScripts, $path2RNAseq, $path2fastqgz, $chrlens, $path2gtfFile)                                = ("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA");
my ($unzip, $qc, $map, $filter, $cleanfiles, $granges)                                                                                          = ("FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE");
my (@sc, @lines2remove)                                                                                                                         = ();

# Find paths to different folders in the Targets.txt file
while(<INPUT>) {
        if (/# My_personal_email\b/) {
                $_ =~ m/"(.+?)"/;
                $email = "$1";
        }
	elsif (/# My_project_title\b/) {
                $_ =~ m/"(.+?)"/;
                $expFolder = "$1";
        }
	elsif (/# Reference_genome\b/) {
                $_ =~ m/"(.+?)"/;
                $genome = "$1";
        }
	elsif (/# Remote_path_to_proj\b/) {
                $_ =~ m/"(.+?)"/;
                $userFolder = "$1";
        }
	elsif (/# Remote_path_to_NEAT\b/) {
                $_ =~ m/"(.+?)"/;
            	$path2NEAT = "$1";
		$path2RNAseq = "$1\/RNApip";
                $path2RNAseqScripts = join("", $path2RNAseq, "/scripts");
        }
	elsif (/# Remote_path_to_orifastq_gz\b/) {
                $_ =~ m/"(.+?)"/;
                $path2fastqgz = "$1";
        }
	elsif (/# Remote_path_to_chrLens_dat\b/) {
                $_ =~ m/"(.+?)"/;
                $chrlens = "$1";
        }
	elsif (/# Remote_path_to_RefGen_fasta\b/) {
                $_ =~ m/"(.+?)"/;
                $refGenome = "$1";
        }
	elsif (/# Remote_path_to_ann_gtf_file\b/) {
                $_ =~ m/"(.+?)"/;
                $path2gtfFile = "$1";
        }
	elsif (/# Aligner_algo_short\b/) {
                $_ =~ m/"(.+?)"/;
                $aligner = "$1";
        }
	elsif (/# Paired_end_seq_run\b/) {
                $_ =~ m/"(.+?)"/;
                $PE = "$1";
        }
	elsif (/# Steps_to_execute_pipe\b/) {
                $_ =~ m/"(.+?)"/;
                @steps2execute = ();
                if (grep /\bunzip\b/i, $_ )             { $unzip                = "TRUE"; push @steps2execute, "Unzip";         }
                if (grep /\bqc\b/i, $_ )                { $qc                   = "TRUE"; push @steps2execute, "QC";            }
                if (grep /\bmap\b/i, $_ )               { $map                  = "TRUE"; push @steps2execute, "Map";           }
                if (grep /\bfilter\b/i, $_ )            { $filter               = "TRUE"; push @steps2execute, "Filter";        }
                if (grep /\bcleanfiles\b/i, $_ )        { $cleanfiles           = "TRUE"; push @steps2execute, "Cleanfiles";    }
                if (grep /\bgranges\b/i, $_ )           { $granges              = "TRUE"; push @steps2execute, "GRanges";	}

        }

} # end of Targets.txt



my $AdvSettings = "$path2expFolder/DataStructure/AdvancedSettings.txt";
open(INPUT, $AdvSettings) || die "Error opening $AdvSettings : $!\n\n\n";

my ($removepcr, $makeunique, $aligncommand)                             = ("NA", "NA", "NA");
my ($SUBkey, $SUBheader, $SUBdependCondition, $SUBcommand)		= ("NA", "NA", "NA", "NA");

while(<INPUT>) {
	if (/# Ressource_manager\b/) {
                $_ =~ m/"(.+?)"/;
		if (grep /\bqsub/i, $_ )        { $SUBkey = "qsub"; $SUBheader  = "QSUB_header.sh";     $QSUBdependCondition = "afterok";    }
                if (grep /\bbsub/i, $_ )        { $SUBkey = "bsub"; $SUBheader  = "BSUB_header.sh";     $BSUBdependCondition = "ended";      }
                $SUBcommand = "$1";
        }
	if (/# Unzip_comand\b/) {
                $_ =~ m/"(.+?)"/;
                $unzipCommand = "$1";
        }
	elsif (/# Zip_file_extension\b/) {
                $_ =~ m/"(.+?)"/;
                $zipExtension = "$1";
        }
        elsif (/# Filter_removePCRdup\b/) {
                $_ =~ m/"(.+?)"/;
                $removepcr = "$1";
        }
	elsif (/# Filter_makeUniqueRead\b/) {
                $_ =~ m/"(.+?)"/;
                $makeunique = "$1";
        }
        elsif (/# Align_command_line_1\b/) {
                $_ =~ m/"(.+?)"/;
                $aligncommand = "$1";
        }
	elsif (/# Wigfile_binSize\b/) {
                $_ =~ m/"(.+?)"/;
                $wigBinSize = "$1";
        }	

} # end of AdvancedSettings.txt


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
my @samples;
foreach $line (@Targets2) {
        $line =~ /^$/ and die "Targets 1: Blank line detected at $.\n\n";
        $line =~ /^[# = " OriFileName FileName OriInpName InpName]/ and next;
        push(@samples, $line);
}
#*----------------------------------------------------------------------*
# Remove duplicated elements in the list @samples and @inputs
%seen           = ();
@samples        = grep { ! $seen{$_} ++ } @samples;

# Remove samples that have "_R2" as these are the paired lanes of "_R1"
my @samples2unzip	= @samples;
my @samplesPE;
my @samplesNoPE;
if( $PE ){

	print "\nPE experiment. \n";
        foreach my $i (0 .. $#samples) {
                if ( grep /\_R2$/, $samples[$i] ){
                        # print "\t '_R2' sample found. \t ($samples[$i]) \n";
                        push(@samplesPE, $samples[$i]);
                }
                else{
                     	# print "\t Main sample found.  \t ($samples[$i]) \n";
                        push(@samplesNoPE, $samples[$i]);
                }
        }
	@samples = @samplesNoPE;
}

#print "\n\n\norisamples:   @orisamples\n";
#print "samples2unzip:   @samples2unzip\n";
#print "samples: \t @samples\n";
#print "samplesNoPE: \t @samplesNoPE\n";
#print "samplesPE: \t @samplesPE\n";

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
print "\n Ressource manager:\t $SUBcommand ($SUBkey)";
print "\n";
print "\n expFolder:\t\t $expFolder";
print "\n genome:\t\t $genome";
print "\n userFolder:\t\t $userFolder";
print "\n path2RNApip:\t\t $path2RNAseq";
print "\n path2expFolder:\t $path2expFolder";
print "\n path2fastq_gz:\t\t $path2fastqgz";
print "\n Targets:\t\t $path2expFolder/DataStructure/Targets.txt";
print "\n chrlens:\t\t $chrlens";
print "\n refGenome:\t\t $refGenome";
print "\n path2gtf:\t\t $path2gtfFile";
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
print "\n Performing following modules:";
print "\n .........................................";
print "\n unzip:\t\t\t $unzip \t ($unzipCommand *filename*$zipExtension)";
print "\n qc:\t\t\t $qc";
print "\n map:\t\t\t $map";
print "\n filter:\t\t $filter";
print "\n cleanfiles:\t\t $cleanfiles";
print "\n granges:\t\t $granges";
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
my $scrhead 			= "$path2RNAseqScripts/$SUBheader";
my $path2iterate		= "$tmpscr/iterate";
my $RNAseqMainIterative		= "$path2iterate/RNAseq.sh";
my $IterateSH			= "$path2iterate/Iterate\_$expFolder.sh";
my $path2CustFct                = "$path2NEAT/CustomFunctions";
my $path2aligned  	        = "$path2expFolder/aligned";
my $path2bam                    = "$path2expFolder/bam";
#my $path2bamRX                  = "$path2expFolder/bam_RX";
my $path2GRanges                = "$path2expFolder/GRangesRData";
#my $path2GRangesRX              = "$path2expFolder/GRangesRData_RX";
my $firstcmd                    = NULL;
my $cmd                         = NULL;
my $cmd2			= NULL;
my $finalcmd                    = NULL;
my $finalcmd2                   = NULL;

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


	#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	my $myJobName		= "unzip";
	my $iterateJobName	= "Iterate_$myJobName\_$expFolder";
	my $path2qsub		= "$tmpscr/$myJobName/$SUBkey";

	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        # Create folders
	unless( -d "$path2fastq" )	{ `mkdir $path2fastq`; }

	#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
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
	
	#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
	print "\n\n-------------------------------------\n\n Unzipping samples: ";
	foreach my $i (0 .. $#samples2unzip) {
		
		print "\n\t $orisamples[$i] \t $samples2unzip[$i] ";		

		# Prepare a personal qsub script
		my $QSUBint  = "$tmpscr/$myJobName/$samples2unzip[$i]\_$myJobName\.sh";
		`cp $scrhead $QSUBint`;

		my $cmd         = "$unzipCommand $path2fastqgz/$orisamples[$i]$zipExtension > $path2fastq/$samples2unzip[$i]\.fastq";
		`echo "$cmd" >> $QSUBint`;
		
		#---------------------------------------------
		# Keep track of the jobs in @myJobs
#		my $jobName	= "Sample_$myJobName$i";
		my $jobName     = "$samples2unzip[$i]_$myJobName$i";

		push(@myJobs, $jobName);
		if($SUBkey =~ "qsub"){  $cmd	= "$jobName=`$SUBcommand -o $path2qsub -e $path2qsub $QSUBint`";}
		if($SUBkey =~ "bsub"){  $cmd	= "$SUBcommand -J $jobName -o $path2qsub\/$jobName.out -e $path2qsub\/$jobName.err $QSUBint";} 
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
	if($SUBkey =~ "qsub"){
                foreach( @myJobs ){ $_ = "\$".$_ ; }
                my $myJobsVec   = join(":", @myJobs);
                $finalcmd	= "FINAL=\`$SUBcommand -N $iterateJobName -o $path2qsub -e $path2qsub -W depend=$QSUBdependCondition\:$myJobsVec $IterateSH`";
        }
	if($SUBkey =~ "bsub"){
                foreach( @myJobs ){ $_ = "$BSUBdependCondition\(\"".$_."\")" ; }
                my $myJobsVec   = join(" && ", @myJobs);
                $finalcmd       = "$SUBcommand -J $iterateJobName -o $path2qsub\/$iterateJobName.out -e $path2qsub\/$iterateJobName.err -w \'$myJobsVec\' $IterateSH";
        }	
	open $QSUB, ">>", "$QSUB" or die "Can't open '$QSUB'";
	print $QSUB "$finalcmd\n";
	close $QSUB;

	#*----------------------------------------------------------------------*
	# Submit jobs to run

	print "\n\n--------------------------------------------------------------------------------------------------\n";
	print "\n Submitting job to $SUBkey cluster: \t `sh $QSUB` \n";
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

	#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	my $myJobName		= "QC";
	my $iterateJobName	= "Iterate_$myJobName\_$expFolder";
	my $path2qsub		= "$tmpscr/$myJobName/$SUBkey";

	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        # Create folders
        unless( -d "$path2QC" ) { `mkdir $path2QC`; }
	`cp $path2RNAseqScripts/QC.R $tmpscr`;
	
	#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
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

	#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
	foreach my $i (0) {
	
		# Prepare a personal qsub script
                my $QSUBint     = "$tmpscr/$myJobName/$myJobName\_qsub.sh";
               `cp $scrhead $QSUBint`;
	
		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		`cp $path2CustFct/QC.R $tmpscr/`;
                my $code	= "$tmpscr/QC.R";
		$cmd		= "Rscript $code $path2expFolder &>> $path2qsub/QCReport.log";
		`echo "$cmd" >> $QSUBint`;
		#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

		#---------------------------------------------
                # Keep track of the jobs in @myJobs
                my $jobName     = "$myJobName$i";

                push(@myJobs, $jobName);
		if($SUBkey =~ "qsub"){  $cmd            = "$jobName=`$SUBcommand -o $path2qsub -e $path2qsub $QSUBint`";}
		if($SUBkey =~ "bsub"){  $cmd		= "$SUBcommand -J $jobName -o $path2qsub\/$jobName.out -e $path2qsub\/$jobName.err $QSUBint";}
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
	# Other jobs do not depend on the completion of this, so no dependencies are needed
	if($SUBkey =~ "qsub"){
                foreach( @myJobs ){ $_ = "\$".$_ ; }
                my $myJobsVec   = join(":", @myJobs);
                $finalcmd	= "FINAL=\`$SUBcommand -N $iterateJobName -o $path2qsub -e $path2qsub $IterateSH`";
        }
	if($SUBkey =~ "bsub"){
                foreach( @myJobs ){ $_ = "$BSUBdependCondition\(\"".$_."\")" ; }
                my $myJobsVec   = join(" && ", @myJobs);
                $finalcmd       = "$SUBcommand -J $iterateJobName -o $path2qsub\/$iterateJobName.out -e $path2qsub\/$iterateJobName.err $IterateSH";
        }

        open $QSUB, ">>", "$QSUB" or die "Can't open '$QSUB'";
        print $QSUB "$finalcmd\n";
        close $QSUB;

	#*----------------------------------------------------------------------*
	# Submit jobs to run

	print "\n\n--------------------------------------------------------------------------------------------------\n";
	print "\n Submitting job to $SUBkey cluster: \t `sh $QSUB` \n";
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

	#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	my $myJobName		= "map";
	my $iterateJobName	= "Iterate_$myJobName\_$expFolder";
	my $myJobName2		= "rename";
	my $path2qsub		= "$tmpscr/$myJobName/$SUBkey";

	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        # Create folders
        unless( -d "$path2Tophat" ) { `mkdir $path2Tophat`; }

	#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
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


	#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
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

			# Lane 2 of sample i should be named the same except having a '_R2' instead od a '_R1'
			$cmd		= "$aligncommand $path2gtfFile -o $path2Tophat/$samples[$i] $path2Tophat/$samplesPE[$i] $refGenome/$genome $path2fastq/$samples[$i].fastq $path2fastq/$samplesPE[$i].fastq";
			`echo "$cmd" >> $QSUBint`;
			#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--			

		} else {
			
			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			$cmd		= "$aligncommand $path2gtfFile -o $path2Tophat/$samples[$i] $refGenome/$genome $path2fastq/$samples[$i].fastq";
			`echo "$cmd" >> $QSUBint`;
			#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

		}

		#---------------------------------------------
		# Rename samples
		#
		$cmd		= "mv $path2Tophat/$samples[$i]/accepted_hits.bam $path2Tophat/$samples[$i]/$samples[$i]\.bam";
		`echo "$cmd" >> $QSUBintRename`;

		#---------------------------------------------
		# Keep track of the jobs in @myJobs
#		my $jobName	= "$myJobName$i";
		my $jobName     = "$samples[$i]_$myJobName$i";
		push(@myJobs, $jobName);
#		my $jobName2	= "$myJobName2$i";
		my $jobName2     = "$samples[$i]_$myJobName2$i";
		push(@myJobs2, $jobName2);
		
		if($SUBkey =~ "qsub"){  $cmd1           = "$jobName=`$SUBcommand -o $path2qsub -e $path2qsub $QSUBint`";}
		if($SUBkey =~ "bsub"){  $cmd1		= "$SUBcommand -J $jobName -o $path2qsub\/$jobName.out -e $path2qsub\/$jobName.err $QSUBint";}

		if($SUBkey =~ "qsub"){
#			foreach( @myJobs2 ){ $_ = "\$".$_ ; }
#			my $myJobsVec2  = join(":", @myJobs2);
                	$cmd2		= "$jobName2=\`$SUBcommand -o $path2qsub -e $path2qsub -W depend=$QSUBdependCondition\:$jobName $QSUBintRename`";
        	}
		if($SUBkey =~ "bsub"){
#			foreach( @myJobs2 ){ $_ = "$BSUBdependCondition\(\"".$_."\")" ; }
#			my $myJobsVec2   = join(" && ", @myJobs2);
			$cmd2		= "$SUBcommand -J $jobName2 -o $path2qsub\/$jobName2.out -e $path2qsub\/$jobName2.err -w $BSUBdependCondition\(\"$jobName\") $QSUBintRename";
        	}

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

	if($SUBkey =~ "qsub"){
                foreach( @myJobs2 ){ $_ = "\$".$_ ; }
                my $myJobsVec   = join(":", @myJobs2);
                $finalcmd	= "FINAL=\`$SUBcommand -N $iterateJobName -o $path2qsub -e $path2qsub -W depend=$QSUBdependCondition\:$myJobsVec $IterateSH`";
        }
	if($SUBkey =~ "bsub"){
                foreach( @myJobs ){ $_ = "$BSUBdependCondition\(\"".$_."\")" ; }
                my $myJobsVec   = join(" && ", @myJobs);
                $finalcmd       = "$SUBcommand -J $iterateJobName -o $path2qsub\/$iterateJobName.out -e $path2qsub\/$iterateJobName.err -w \'$myJobsVec\' $IterateSH";
        }
	open $QSUB, ">>", "$QSUB" or die "Can't open '$QSUB'";
	print $QSUB "$finalcmd\n";
	close $QSUB;

	#*----------------------------------------------------------------------*
	# Submit jobs to run
	
	print "\n\n--------------------------------------------------------------------------------------------------\n";
	print "\n Submitting job to $SUBkey cluster: \t `sh $QSUB` \n";
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

	#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	my $myJobName		= "filter";
	my $iterateJobName	= "Iterate_$myJobName\_$expFolder";
	my $path2qsub		= "$tmpscr/$myJobName/$SUBkey";

	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        # Create folders

	#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
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

	#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
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
		$cmd		= "mv $path2Tophat/$samples[$i]/$samples[$i]\.bam $path2Tophat/$samples[$i]/$samples[$i]\_unsorted.bam";
                `echo "$cmd" >> $QSUBint`;
		if( $removepcr ){
			# bam to sorted bam
			$cmd		= "samtools sort $path2Tophat/$samples[$i]/$samples[$i]\_unsorted.bam $path2Tophat/$samples[$i]/$samples[$i]\_sortedwpcr";
			`echo "$cmd" >> $QSUBint`;
			$cmd		= "samtools rmdup -s $path2Tophat/$samples[$i]/$samples[$i]\_sortedwpcr.bam $path2Tophat/$samples[$i]/$samples[$i].bam";
			`echo "$cmd" >> $QSUBint`;
			
		} else {
			# bam to sorted bam
			$cmd = "samtools sort $path2Tophat/$samples[$i]/$samples[$i]\_unfiltered.bam $path2Tophat/$samples[$i]/$samples[$i].bam";
			`echo "$cmd" >> $QSUBint`;
		}

		$cmd= "samtools index $path2Tophat/$samples[$i]/$samples[$i].bam $path2Tophat/$samples[$i]/$samples[$i].bai";
		`echo "$cmd" >> $QSUBint`;
		#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--


		print "\n\n Samtools processing for $samples[$i] done. \n\n";
		
		#---------------------------------------------
		# Keep track of the jobs in @myJobs
#		my $jobName	= "$myJobName$i";
		my $jobName     = "$samples[$i]_$myJobName$i";

		push(@myJobs, "$jobName");
		#$cmd        	= "$jobName=`qsub -o $path2qsub -e $path2qsub $QSUBint`";
		if($SUBkey =~ "qsub"){  $cmd           = "$jobName=`$SUBcommand -o $path2qsub -e $path2qsub $QSUBint`";}
		if($SUBkey =~ "bsub"){  $cmd           = "$SUBcommand -J $jobName -o $path2qsub\/$jobName.out -e $path2qsub\/$jobName.err $QSUBint";}
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
	if($SUBkey =~ "qsub"){
                foreach( @myJobs ){ $_ = "\$".$_ ; }
                my $myJobsVec   = join(":", @myJobs);
                $finalcmd	= "FINAL=\`$SUBcommand -N $iterateJobName -o $path2qsub -e $path2qsub -W depend=$QSUBdependCondition\:$myJobsVec $IterateSH`";
        }
	if($SUBkey =~ "bsub"){
                foreach( @myJobs ){ $_ = "$BSUBdependCondition\(\"".$_."\")" ; }
                my $myJobsVec   = join(" && ", @myJobs);
                $finalcmd       = "$SUBcommand -J $iterateJobName -o $path2qsub\/$iterateJobName.out -e $path2qsub\/$iterateJobName.err -w \'$myJobsVec\' $IterateSH";
        }

	open $QSUB, ">>", "$QSUB" or die "Can't open '$QSUB'";
	print $QSUB "$finalcmd\n";
	close $QSUB;

	#*----------------------------------------------------------------------*
	# Submit jobs to run

	print "\n\n--------------------------------------------------------------------------------------------------\n";
	print "\n Submitting job to $SUBkey cluster: \t `sh $QSUB` \n";
	`sh $QSUB`;

	#*----------------------------------------------------------------------*
	# Exit script

	print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	print "\n Exiting $myJobName section with no known error \n";
	print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n";

	exit 0;
}


if($cleanfiles =~ "TRUE"){

	#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        my $myJobName           = "cleanfiles";
	my $iterateJobName	= "Iterate_$myJobName\_$expFolder";
        my $path2qsub           = "$tmpscr/$myJobName/$SUBkey";
	my $myJobName2          = "rename_Tophat";
	my $iterateJobName2	= "Iterate_$myJobName2\_$expFolder";

	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        # Create folders
#	unless( -d "$path2aligned" )    { `mkdir $path2aligned`;        }
        unless( -d "$path2bam" )        { `mkdir $path2bam`;            }

	#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        # Create file to store jobs in
        unless( -d "$tmpscr/$myJobName" )	{ `mkdir $tmpscr/$myJobName`; }
        unless( -d "$path2qsub" )               { `mkdir $path2qsub`; }
        my $QSUB        = "$tmpscr/$myJobName/$myJobName\.sh";
        open $QSUB, ">", "$QSUB" or die "Can't open '$QSUB'";
        print $QSUB "#!/bin/bash\n";
        close $QSUB;
        `chmod 777 $QSUB`;
        print "\n Store all of the following '$myJobName' jobs in $QSUB \n";
        my @myJobs;

	 #*----------------------------------------------------------------------*
        # Rename the 'Tophat' folder to 'aligned'
        my $QSUBintRename	= "$tmpscr/$myJobName/$myJobName2\.sh";
        open $QSUBintRename, ">", "$QSUBintRename" or die "Can't open '$QSUBintRename'";
        print $QSUBintRename "#!/bin/bash\n";
        close $QSUBintRename;
        `chmod 777 $QSUBintRename`;


        print "\n\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n";
        print "\n Moving .bam and .bai files from the Tophat to the aligned folder \n";

	#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        foreach my $i (0 .. $#samples) {

                # Prepare a personal qsub script
                my $QSUBint     = "$tmpscr/$myJobName/$samples[$i]\_$myJobName\.sh";
               `cp $scrhead $QSUBint`;

		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		my $cmd = "`cp $path2Tophat/$samples[$i]/$samples[$i].bam $path2bam/`";
		`echo "$cmd" >> $QSUBint`;
                $cmd = "`cp $path2Tophat/$samples[$i]/$samples[$i].bai $path2bam/`";
		`echo "$cmd" >> $QSUBint`;
                #--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

                #---------------------------------------------
                # Keep track of the jobs in @myJobs
#               my $jobName     = "$myJobName$i";
		my $jobName     = "$samples[$i]_$myJobName$i";

                push(@myJobs, "$jobName");
		if($SUBkey =~ "qsub"){  $cmd           = "$jobName=`$SUBcommand -o $path2qsub -e $path2qsub $QSUBint`";}
		if($SUBkey =~ "bsub"){  $cmd           = "$SUBcommand -J $jobName -o $path2qsub\/$jobName.out -e $path2qsub\/$jobName.err $QSUBint";}
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
        # Prepar renaming job
	my $cmd = "mv $path2Tophat/ $path2aligned";
	`echo "$cmd" >> $QSUBintRename`;
	
	#*----------------------------------------------------------------------*
        # Prepar file containing the jobs to run

	# Add the next job line to the $QSUB
	if($SUBkey =~ "qsub"){
                foreach( @myJobs ){ $_ = "\$".$_ ; }
                my $myJobsVec   = join(":", @myJobs);
                $finalcmd	= "RENAME=\`$SUBcommand -N $iterateJobName2 -o $path2qsub -e $path2qsub -W depend=$QSUBdependCondition\:$myJobsVec $QSUBintRename`";
        }
	if($SUBkey =~ "bsub"){
                foreach( @myJobs ){ $_ = "$BSUBdependCondition\(\"".$_."\")" ; }
                my $myJobsVec   = join(" && ", @myJobs);
                $finalcmd       = "$SUBcommand -J $iterateJobName2 -o $path2qsub\/$iterateJobName2.out -e $path2qsub\/$iterateJobName2.err -w \'$myJobsVec\' $QSUBintRename";
        }

	if($SUBkey =~ "qsub"){
                $finalcmd2	= "FINAL=\`$SUBcommand -N $iterateJobName -o $path2qsub -e $path2qsub -W depend=$QSUBdependCondition\:\$RENAME $IterateSH`";
        }	
	if($SUBkey =~ "bsub"){
                $finalcmd2       = "$SUBcommand -J $iterateJobName -o $path2qsub\/$iterateJobName.out -e $path2qsub\/$iterateJobName.err -w $BSUBdependCondition\(\"$iterateJobName2\") $IterateSH";
        }

        open $QSUB, ">>", "$QSUB" or die "Can't open '$QSUB'";
        print $QSUB "$finalcmd\n";
        print $QSUB "$finalcmd2\n";
	close $QSUB;

        #*----------------------------------------------------------------------*
        # Submit jobs to run

        print "\n\n--------------------------------------------------------------------------------------------------\n";
        print "\n Submitting job to $SUBkey cluster: \t `sh $QSUB` \n";
        `sh $QSUB`;

        #*----------------------------------------------------------------------*
        # Exit script

        print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
        print "\n Exiting $myJobName section with no known error \n";
        print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n";

        exit 0;

}


#*----------------------------------------------------------------------*
# Create GRanges objects

if( $granges =~ "TRUE" ){

        print "\n\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n";
        print "\n Creating GRanges\n";

	#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        my $myJobName           = "granges";
	my $iterateJobName	= "Iterate_$myJobName\_$expFolder";
        my $path2qsub           = "$tmpscr/$myJobName/$SUBkey";

	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        # Create folders
	unless( -d "$path2GRanges" )	{ `mkdir $path2GRanges`; }

	#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
        # Create file to store jobs in
        unless( -d "$tmpscr/$myJobName" )	{ `mkdir $tmpscr/$myJobName`; }
        unless( -d "$path2qsub" )               { `mkdir $path2qsub`; }
        my $QSUB        = "$tmpscr/$myJobName/$myJobName\.sh";
        open $QSUB, ">", "$QSUB" or die "Can't open '$QSUB'";
        print $QSUB "#!/bin/bash\n";
        close $QSUB;
        `chmod 777 $QSUB`;
        print "\n Store all of the following '$myJobName' jobs in $QSUB \n";
        my @myJobs;

	#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
	foreach my $i (0) {
                #-----------------------------------------------------------
                # Prepare a personal qsub script
                my $QSUBint  = "$tmpscr/$myJobName/$myJobName\_qsub.sh";
                `cp $scrhead $QSUBint`;

                #-----------------------------------------------------------
                # Parameters
		`cp $path2CustFct/Bam2GRangesRemote.R $tmpscr/`;
                my $code	= "$tmpscr/Bam2GRangesRemote.R";

                #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                my $cmd         = "Rscript $code $path2expFolder $path2bam $path2GRanges $path2CustFct $wigBinSize &>> $path2qsub/GRanges.log";
                `echo "$cmd" >> $QSUBint`;
                #--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

		#---------------------------------------------
                # Keep track of the jobs in @myJobs
                my $jobName     = "$myJobName$i";
                push(@myJobs, "$jobName");
		if($SUBkey =~ "qsub"){  $cmd           = "$jobName=`$SUBcommand -o $path2qsub -e $path2qsub $QSUBint`";}
		if($SUBkey =~ "bsub"){  $cmd           = "$SUBcommand -J $jobName -o $path2qsub\/$jobName.out -e $path2qsub\/$jobName.err $QSUBint";}
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
	if($SUBkey =~ "qsub"){
                foreach( @myJobs ){ $_ = "\$".$_ ; }
                my $myJobsVec   = join(":", @myJobs);
                $finalcmd	= "FINAL=\`$SUBcommand -N $iterateJobName -o $path2qsub -e $path2qsub -W depend=$QSUBdependCondition\:$myJobsVec $IterateSH`";
        }
	if($SUBkey =~ "bsub"){
                foreach( @myJobs ){ $_ = "$BSUBdependCondition\(\"".$_."\")" ; }
                my $myJobsVec   = join(" && ", @myJobs);
                $finalcmd       = "$SUBcommand -J $iterateJobName -o $path2qsub\/$iterateJobName.out -e $path2qsub\/$iterateJobName.err -w \'$myJobsVec\' $IterateSH";
        }

	open $QSUB, ">>", "$QSUB" or die "Can't open '$QSUB'";
        print $QSUB "$finalcmd\n";
        close $QSUB;

        #*----------------------------------------------------------------------*
        # Submit jobs to run

        print "\n\n--------------------------------------------------------------------------------------------------\n";
        print "\n Submitting job to $SUBkey cluster: \t `sh $QSUB` \n";
        `sh $QSUB`;

        #*----------------------------------------------------------------------*
        # Exit script

        print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
        print "\n Exiting $myJobName section with no known error \n";
        print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n";

        exit 0;

}

exit 0;


print "\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n\n";

#*----------------------------------------------------------------------*



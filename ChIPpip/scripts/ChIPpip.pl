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

my ($expFolder, $genome, $genomeRX, $userFolder, $path2ChIPseqScripts, $path2ChIPseq, $path2fastqgz)	= ("NA", "NA", "NA", "NA", "NA", "NA", "NA");
my ($unzip, $qc, $chiprx, $map, $filter, $peakcalling, $cleanbigwig, $cleanfiles, $granges)		= ("FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE");
my (@sc, @lines2remove)											= ();

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
	elsif (/# Reference_genome_rx\b/) {
                $_ =~ m/"(.+?)"/;
                $genomeRX = "$1";
        }
	elsif (/# Remote_path_to_proj\b/) {
                $_ =~ m/"(.+?)"/;
                $userFolder = "$1";
        }
	elsif (/# Remote_path_to_NEAT\b/) {
                $_ =~ m/"(.+?)"/;
		$path2NEAT = "$1";
                $path2ChIPseq = "$1\/ChIPpip";
                $path2ChIPseqScripts = join("", $path2ChIPseq, "/scripts");
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
	elsif (/# Remote_path_to_chrLens_dat_ChIP_rx\b/) {
                $_ =~ m/"(.+?)"/;
                $chrlensRX = "$1";
        }
	elsif (/# Remote_path_to_RefGen_fasta_ChIP_rx\b/) {
                $_ =~ m/"(.+?)"/;
                $refGenomeRX = "$1";
        }
	elsif (/# Aligner_algo_short\b/) {
                $_ =~ m/"(.+?)"/;
		if (grep /\bbwa\b/i, $_ )	{ $aligner	= "BWA"; }
		if (grep /\bbowtie\b/i, $_ )	{ $aligner      = "BOWTIE"; }
		if (grep /\bbowtie2\b/i, $_ )	{ $aligner      = "BOWTIE"; }
        }
	elsif (/# Paired_end_seq_run\b/) {
                $_ =~ m/"(.+?)"/;
                $PE = "$1";
        }
	elsif (/# Remove_from_bigwig\b/) {
                $_ =~ m/"(.+?)"/;
                my $text = "$1";
                my @var = split(",", $text);
                foreach my $line (@var) {
                        $line =~ s/\s+//g;
                        push(@lines2remove, $line);
                }
        }
	elsif (/# PeakCaller_R_script\b/) {
                $_ =~ m/"(.+?)"/;
                $peakcaller = "$1";
        }

	elsif (/# Steps_to_execute_pipe\b/) {
                $_ =~ m/"(.+?)"/;
                @steps2execute = ();
                if (grep /\bunzip\b/i, $_ )             { $unzip                = "TRUE"; push @steps2execute, "Unzip";         }
                if (grep /\bqc\b/i, $_ )                { $qc                   = "TRUE"; push @steps2execute, "QC";            }
                if (grep /\bchiprx\b/i, $_ )            { $chiprx               = "TRUE"; push @steps2execute, "ChIPrx";        }
                if (grep /\bmap\b/i, $_ )               { $map                  = "TRUE"; push @steps2execute, "Map";           }
                if (grep /\bfilter\b/i, $_ )            { $filter               = "TRUE"; push @steps2execute, "Filter";        }
                if (grep /\bpeakcalling\b/i, $_ )	{ $peakcalling          = "TRUE"; push @steps2execute, "Peakcalling";   }
                if (grep /\bcleanbigwig\b/i, $_ )	{ $cleanbigwig          = "TRUE"; push @steps2execute, "Cleanbigwig";   }
		if (grep /\bcleanfiles\b/i, $_ )	{ $cleanfiles		= "TRUE"; push @steps2execute, "Cleanfiles";	}
		if (grep /\bgranges\b/i, $_ )		{ $granges		= "TRUE"; push @steps2execute, "GRanges";	}
        }

} # end of Targets.txt



my $AdvSettings = "$path2expFolder/DataStructure/AdvancedSettings.txt";
open(INPUT, $AdvSettings) || die "Error opening $AdvSettings : $!\n\n\n";

my ($removepcr, $makeunique, $aligncommand1, $fdr, $posopt, $densityopt, $enforceisize)		= ("NA", "NA", "NA", "NA", "NA", "NA", "NA");
my ($SUBkey, $SUBheader, $SUBdependCondition, $SUBcommand)					= ("NA", "NA", "NA", "NA");

while(<INPUT>) {
	if (/# Ressource_manager\b/) {
                $_ =~ m/"(.+?)"/;
                if (grep /\bqsub/i, $_ )        { $SUBkey = "qsub"; $SUBheader	= "QSUB_header.sh";	$SUBdependCondition = "afterok";	}
                if (grep /\bbsub/i, $_ )        { $SUBkey = "bsub"; $SUBheader	= "BSUB_header.sh";	$SUBdependCondition = "done";		}
                $SUBcommand = "$1";
        }
	elsif (/# Unzip_comand/) {
                $_ =~ m/"(.+?)"/;
                $unzipCommand = "$1";
        }
        elsif (/# Zip_file_extension\b/) {
                $_ =~ m/"(.+?)"/;
                $zipExtension = "$1";
	}
        elsif (/# Align_command_line_1\b/) {
                $_ =~ m/"(.+?)"/;
                $aligncommand1 = "$1";
	}
        elsif (/# Filter_removePCRdup\b/) {
                $_ =~ m/"(.+?)"/;
                $removepcr = "$1";
	}
        elsif (/# Filter_makeUniqueRead\b/) {
                $_ =~ m/"(.+?)"/;
                $makeunique = "$1";
        }
        elsif (/# PeakCaller_fdr\b/) {
                $_ =~ m/"(.+?)"/;
                $fdr = "$1";
        }
        elsif (/# PeakCaller_posopt\b/) {
                $_ =~ m/"(.+?)"/;
                $posopt = "$1";
        }
        elsif (/# PeakCaller_densityopt\b/) {
                $_ =~ m/"(.+?)"/;
                $densityopt = "$1";
        }
        elsif (/# PeakCaller_enfSize\b/) {
                $_ =~ m/"(.+?)"/;
                $enforceisize = "$1";
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
my @Targets3 = `cut -f3 $Targets`;
        chomp(@Targets3);
my @Targets4 = `cut -f4 $Targets`;
        chomp(@Targets4);

# Store original file names in orisamples
my @orisamples;
foreach $line (@Targets1) {
        $line =~ /^$/ and die "Targets 1: Blank line detected at $.\n\n";
        $line =~ /^[# : = " OriFileName FileName OriInpName InpName]/ and next;
        push(@orisamples, $line);
}
my @samples;
foreach $line (@Targets2) {
        $line =~ /^$/ and die "Targets 1: Blank line detected at $.\n\n";
        $line =~ /^[# : = " OriFileName FileName OriInpName InpName]/ and next;
        push(@samples, $line);
}
my @oriinputs;
foreach $line (@Targets3) {
        $line =~ /^$/ and die "Targets 3: Blank line detected at $.\n\n";
        $line =~ /^[# : = " OriFileName FileName OriInpName InpName]/ and next;
        push(@oriinputs, $line);
}
my @inputs;
foreach $line (@Targets4) {
        $line =~ /^$/ and next;
        $line =~ /^[# : = " OriFileName FileName OriInpName InpName]/ and next;
        push(@inputs, $line);
}

#*----------------------------------------------------------------------*
# Remove duplicated elements in the list @samples and @inputs
%seen           = ();
@samples        = grep { ! $seen{$_} ++ } @samples;
@allinputs	= @inputs;
%seen           = ();
@inputs         = grep { ! $seen{$_} ++ } @inputs;

# Remove samples that have "_R2" as these are the paired lanes of "_R1"
my @samplesPE;
my @samplesNoPE;
my @inputsPE;
my @inputsNoPE;
my @samplesInputsPE;
my @samples2unzip = @samples;

if( $PE ){
	print "\nPE experiment. \n";
        foreach my $i (0 .. $#samples) {
                if ( grep /\_R2$/, $samples[$i] ){
                        #print "\t '_R2' sample found. \t ($samples[$i]) \n";
                        push(@samplesPE, $samples[$i]);
                }
                else{
                       #print "\t Main sample found.  \t ($samples[$i]) \n";
                        push(@samplesNoPE, $samples[$i]);
                }
        }
	@samples = @samplesNoPE;
	@samplesInputsPE = @samplesPE;

	foreach my $i (0 .. $#inputs) {
                if ( grep /\_R2$/, $inputs[$i] ){
                        push(@inputsPE, $inputs[$i]);
                }
                else{
                        push(@inputsNoPE, $inputs[$i]);
                }
        }
	@inputs = @inputsNoPE;
	push(@samplesInputsPE, @inputsPE)
}

#print "\n\n\norisamples:   @orisamples\n";
#print "samples2unzip:   @samples2unzip\n";
#print "samples: \t @samples\n";
#print "samplesNoPE: \t @samplesNoPE\n";
#print "samplesPE: \t @samplesPE\n";

@samplesInputs  = @samples;
push (@samplesInputs, @inputs);

#*----------------------------------------------------------------------*
# Store variables into @samples
my $cutoff	= 0.0;
my @chrs        = `cut -f1 $chrlens`;
chomp(@chrs);

#*----------------------------------------------------------------------*
# Define paths
my $path2expFolder 	= "$userFolder/$expFolder";
$Targets 		= "$path2expFolder/DataStructure/Targets.txt";
#$AdvancedSettings	= "$path2expFolder/DataStructure/AdvancedSettings.txt";

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
print "\n path2ChIPpip:\t\t $path2ChIPseq";
print "\n path2expFolder:\t $path2expFolder";
print "\n path2fastq.gz:\t\t $path2fastqgz";
print "\n Targets:\t\t $path2expFolder/DataStructure/Targets.txt";
print "\n chrlens:\t\t $chrlens";
print "\n refGenome:\t\t $refGenome";
if($chiprx =~ "TRUE"){
        print "\n\t - chrlensRX:\t $chrlensRX";
        print "\n\t - refGenomeRX:\t $refGenomeRX";
}
print "\n";
print "\n Unzip command:\t\t $unzipCommand";
print "\n Paired end sequencing:\t $PE";
print "\n Aligner algorithm:\t $aligner";
print "\n Remove pcr dupl:\t $removepcr";
print "\n Make unique reads:\t $makeunique";
print "\n Peak caller algo:\t $peakcaller";
print "\n";
#print "\n Current working dir:\t $path2expFolder";
#print "\n";
print "\n .........................................";
print "\n Performing following modules:";
print "\n .........................................";
print "\n unzip:\t\t\t $unzip \t ($unzipCommand *filename*$zipExtension)";
print "\n qc:\t\t\t $qc";
print "\n chiprx:\t\t $chiprx \t ($genomeRX)";
print "\n map:\t\t\t $map \t ($genome) ($aligncommand1)";
print "\n filter:\t\t $filter";
print "\n peakcalling:\t\t $peakcalling \t ($peakcaller)";
print "\n cleanbigwig:\t\t $cleanbigwig \t (remove: @lines2remove)";
print "\n cleanfiles:\t\t $cleanfiles";
print "\n granges:\t\t $granges";
print "\n .........................................";
print "\n";
print "\n Samples: ";
foreach my $i (0 .. $#samples) {
        print "\n\t $samples[$i] \t - \t $inputs[$i]";
}
print "\n";
#print "\n----------------------------------------\n";



#*----------------------------------------------------------------------*
# Set different paths

my $tmpscr			= "$path2expFolder/scripts";
my $path2fastq			= "$path2expFolder/fastq";
my $path2QC			= "$path2expFolder/QC";
my $path2aligned		= "$path2expFolder/aligned";
my $scrhead			= "$path2ChIPseqScripts/$SUBheader";
my $path2iterate		= "$tmpscr/iterate";
my $ChIPseqMainIterative	= "$path2iterate/RNAseq.sh";
my $IterateSH			= "$path2iterate/Iterate\_$expFolder.sh";
my $path2CustFct		= "$path2NEAT/CustomFunctions";
my $path2DataStructure		= "$path2expFolder/DataStructure";
my $path2peakcalling		= "$path2expFolder/peakcalling";
my $path2bam			= "$path2expFolder/bam";
my $path2bamRX			= "$path2expFolder/bam_RX";
my $path2GRanges		= "$path2expFolder/GRangesRData";
my $path2GRangesRX		= "$path2expFolder/GRangesRData_RX";
my $firstcmd    		= NULL;
my $cmd				= NULL;
my $finalcmd			= NULL; 
my $shellHeader			= "#!/bin/bash\n";

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
	foreach my $i (0 .. $#samples2unzip) {

		# Prepare a personal qsub script
		my $QSUBint  = "$tmpscr/$myJobName/$samples2unzip[$i]\_$myJobName\.sh";
		`cp $scrhead $QSUBint`;
		`chmod 777 $QSUBint`;

		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		$cmd		= "gunzip -c $path2fastqgz/$orisamples[$i]$zipExtension > $path2fastq/$samples2unzip[$i]\.fastq";
		`echo "$cmd" >> $QSUBint`;
		#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--
		
		#---------------------------------------------
		# Keep track of the jobs in @myJobs
		my $jobName	= "Sample_$myJobName$i";
#		my $jobName     = "$samples2unzip[$i]_$myJobName$i";
		push(@myJobs, $jobName);
		if($SUBkey =~ "qsub"){	$cmd	= "$jobName=`$SUBcommand -o $path2qsub -e $path2qsub $QSUBint`";}
		if($SUBkey =~ "bsub"){  $cmd	= "$SUBcommand -J $jobName -o $path2qsub\/$jobName.out -e $path2qsub\/$jobName.err $QSUBint";}
		open $QSUB, ">>", "$QSUB" or die "Can't open '$QSUB'";
		print $QSUB "$cmd\n";
		close $QSUB;     
	}

	#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
	foreach my $i (0 .. $#inputs) {

		# Prepare a personal qsub script
		my $QSUBint  = "$tmpscr/$myJobName/$inputs[$i]\_$myJobName\.sh";
		`cp $scrhead $QSUBint`;
		`chmod 777 $QSUBint`;

		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		my $cmd         = "$unzipCommand $path2fastqgz/$oriinputs[$i]$zipExtension > $path2fastq/$inputs[$i]\.fastq";
		`echo "$cmd" >> $QSUBint`;
		#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

		#---------------------------------------------
		# Keep track of the jobs in @myJobs
		my $jobName     = "Inp_$myJobName$i";
#		my $jobName     = "$inputs[$i]_$myJobName$i";
		push(@myJobs, $jobName);
		if($SUBkey =~ "qsub"){	$cmd	= "$jobName=`$SUBcommand -o $path2qsub -e $path2qsub $QSUBint`";}
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
		my $myJobsVec	= join(":", @myJobs);
		$finalcmd	= "FINAL=\`$SUBcommand -o $path2qsub -e $path2qsub -W depend=$SUBdependCondition\:$myJobsVec $IterateSH`";
	}
	if($SUBkey =~ "bsub"){
		foreach( @myJobs ){ $_ = "$SUBdependCondition\(\"".$_."\")" ; }
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
        `cp $path2ChIPseqScripts/QC.R $tmpscr`;

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
		`chmod 777 $QSUBint`;

		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		`cp $path2CustFct/QC.R $tmpscr/`;
		my $code	= "$tmpscr/QC.R";
		my $cmd		= "Rscript $code $path2expFolder &>> $path2qsub/QCReport.log";
		`echo "$cmd" >> $QSUBint`;
		#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

		#---------------------------------------------
                # Keep track of the jobs in @myJobs
                my $jobName     = "$myJobName$i";
                push(@myJobs, $jobName);
		if($SUBkey =~ "qsub"){  $cmd	= "$jobName=`$SUBcommand -o $path2qsub -e $path2qsub $QSUBint`";}
		if($SUBkey =~ "bsub"){  $cmd    = "$SUBcommand -J $jobName -o $path2qsub\/$jobName.out -e $path2qsub\/$jobName.err $QSUBint";}
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

	# Add the next job iteration
	# Other	jobs do	not depend on the completion of	this, so no dependencies are needed
	# Add the next job line to the $mapQSUB
        if($SUBkey =~ "qsub"){
                foreach( @myJobs ){ $_ = "\$".$_ ; }
                my $myJobsVec   = join(":", @myJobs);
                $finalcmd	= "FINAL=\`$SUBcommand -N $iterateJobName -o $path2qsub -e $path2qsub $IterateSH`";
        }
	 if($SUBkey =~ "bsub"){
                foreach( @myJobs ){ $_ = "$SUBdependCondition\(\"".$_."\")" ; }
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
# Mapping sequences with bwa onto chiprx genome

if( $chiprx =~ "TRUE" ){

	print "\n\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n";
        print "\n Mapping fastq files to the ChIP-RX genome\n";

	#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        my $myJobName           = "chiprx";
	my $iterateJobName	= "Iterate_$myJobName\_$expFolder";
	my $myJobName2		= "chiprx_clean";
        my $path2qsub           = "$tmpscr/$myJobName/$SUBkey";

	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        # Create folders
	unless( -d "$path2aligned" )	{ `mkdir $path2aligned`; }

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
        foreach my $i (0 .. $#samplesInputs) {

                # Prepare a personal qsub script
                my $QSUBint     = "$tmpscr/$myJobName/$samplesInputs[$i]\_$myJobName\.sh";
                `cp $scrhead $QSUBint`;
		`chmod 777 $QSUBint`;

                # Create a directory
                my $path2currentSampleDir	= "$path2aligned/$samplesInputs[$i]\_RX";
                unless( -d "$path2currentSampleDir" )   { `mkdir $path2currentSampleDir`; }


                if( $PE ) {
			
			if($aligner =~ "BWA"){
                        	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                        	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                        	my $cmd         = "$aligncommand1 $refGenomeRX $path2fastq/$samplesInputs[$i]\.fastq > $path2currentSampleDir/$samplesInputs[$i]\.sai";
                        	`echo "$cmd" >> $QSUBint`;
                	        $cmd            = "$aligncommand1 $refGenomeRX $path2fastq/$samplesInputsPE[$i]\.fastq > $path2currentSampleDir/$samplesInputsPE[$i]\.sai";
        	                `echo "$cmd" >> $QSUBint`;
	                        $cmd            = "bwa sampe $refGenomeRX $path2currentSampleDir/$samplesInputs[$i]\.sai $path2currentSampleDir/$samplesInputsPE[$i]\.sai $path2fastq/$samplesInputs[$i]\_1.fastq $path2fastq/$samplesInputsPE[$i]\.fastq > $path2currentSampleDir/$samplesInputs[$i]\.sam";
	                        `echo "$cmd" >> $QSUBint`;
				# Uniquely mapped reads
                	        $cmd         = "grep -E '\\sX0:i:1\\s' $path2currentSampleDir/$samplesInputs[$i]\.sam > $path2currentSampleDir/$samplesInputs[$i]\.u.sam";
                        	`echo "$cmd" >> $QSUBint`;
                        	# sam to bam
                	        $cmd            = "samtools view -b $path2currentSampleDir/$samplesInputs[$i]\.u.sam -T $refGenomeRX -o $path2currentSampleDir/$samplesInputs[$i]\.bam";
        	                `echo "$cmd" >> $QSUBint`;
	                        #--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--
			}

                } else {
			
			if($aligner =~ "BWA"){			
                        	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                        	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                        	my $cmd         = "$aligncommand1 $refGenomeRX $path2fastq/$samplesInputs[$i]\.fastq > $path2currentSampleDir/$samplesInputs[$i]\.sai";
                        	`echo "$cmd" >> $QSUBint`;
                        	my $cmd2        = "`bwa samse $refGenomeRX $path2currentSampleDir/$samplesInputs[$i]\.sai $path2fastq/$samplesInputs[$i]\.fastq > $path2currentSampleDir/$samplesInputs[$i]\.sam`";
                        	`echo "$cmd2" >> $QSUBint`;
				# Uniquely mapped reads
				$cmd         = "`grep -E '\\sX0:i:1\\s' $path2currentSampleDir/$samplesInputs[$i]\.sam > $path2currentSampleDir/$samplesInputs[$i]\.u.sam`";
                        	`echo "$cmd" >> $QSUBint`;
				# sam to bam
                		$cmd		= "`samtools view -b $path2currentSampleDir/$samplesInputs[$i]\.u.sam -T $refGenomeRX -o $path2currentSampleDir/$samplesInputs[$i]\.bam`";
                		`echo "$cmd" >> $QSUBint`;
                        	#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--
			}
                }

                #---------------------------------------------
                # Keep track of the jobs in @myJobs
#		my $jobName     = "$myJobName$i";
		my $jobName     = "$samplesInputs[$i]_$myJobName$i";
                push(@myJobs, $jobName);
		if($SUBkey =~ "qsub"){  $cmd	= "$jobName=`$SUBcommand -o $path2qsub -e $path2qsub $QSUBint`";}
		if($SUBkey =~ "bsub"){  $cmd    = "$SUBcommand -J $jobName -o $path2qsub\/$jobName.out -e $path2qsub\/$jobName.err $QSUBint";}
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
                $finalcmd	= "FINAL=\`$SUBcommand -N $iterateJobName -o $path2qsub -e $path2qsub -W depend=$SUBdependCondition\:$myJobsVec $IterateSH`";
        }
	 if($SUBkey =~ "bsub"){
                foreach( @myJobs ){ $_ = "$SUBdependCondition\(\"".$_."\")" ; }
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
# Mapping sequences with bwa

if( $map =~ "TRUE" ){	

	print "\n\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n";
	print "\n Mapping fastq files\n";

	#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	my $myJobName		= "map";
	my $iterateJobName	= "Iterate_$myJobName\_$expFolder";
	my $path2qsub		= "$tmpscr/$myJobName/$SUBkey";

	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        # Create folders
        unless( -d "$path2aligned" )    { `mkdir $path2aligned`; }

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
	foreach my $i (0 .. $#samplesInputs) {
		
		# Prepare a personal qsub script
		my $QSUBint	= "$tmpscr/$myJobName/$samplesInputs[$i]\_$myJobName\.sh";
		`cp $scrhead $QSUBint`;
		`chmod 777 $QSUBint`;
		
		# Create a directory
		my $path2currentSampleDir	= "$path2aligned/$samplesInputs[$i]";
		unless( -d "$path2currentSampleDir" )	{ `mkdir $path2currentSampleDir`; }


		if( $PE ) {

			if($aligner =~ "BWA"){

				#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
				#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
				my $cmd		= "$aligncommand1 $refGenome $path2fastq/$samplesInputs[$i]\.fastq > $path2currentSampleDir/$samplesInputs[$i]\.sai";
				`echo "$cmd" >> $QSUBint`;
				$cmd		= "$aligncommand1 $refGenome $path2fastq/$samplesInputsPE[$i]\.fastq > $path2currentSampleDir/$samplesInputsPE[$i]\.sai";
				`echo "$cmd" >> $QSUBint`;
				$cmd		= "bwa sampe $refGenome $path2currentSampleDir/$samplesInputs[$i]\.sai $path2currentSampleDir/$samplesInputsPE[$i]\.sai $path2fastq/$samplesInputs[$i]\.fastq $path2fastq/$samplesInputsPE[$i]\.fastq > $path2currentSampleDir/$samplesInputs[$i]\.sam";
				`echo "$cmd" >> $QSUBint`;
				#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--
			}

			if($aligner =~ "BOWTIE"){

				#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                                #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
				my $cmd         = "$aligncommand1 $refGenome -1 $path2fastq/$samplesInputs[$i]\.fastq -2 $path2fastq/$samplesInputsPE[$i]\.fastq -S $path2currentSampleDir/$samplesInputs[$i]\.sam";
                                `echo "$cmd" >> $QSUBint`;
                                #--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

			}

		} else {

			if($aligner =~ "BWA"){
				#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
				#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
				my $cmd		= "$aligncommand1 $refGenome $path2fastq/$samplesInputs[$i]\.fastq > $path2currentSampleDir/$samplesInputs[$i]\.sai";
				`echo "$cmd" >> $QSUBint`;
				my $cmd2	= "bwa samse $refGenome $path2currentSampleDir/$samplesInputs[$i]\.sai $path2fastq/$samplesInputs[$i]\.fastq > $path2currentSampleDir/$samplesInputs[$i]\.sam";
				`echo "$cmd2" >> $QSUBint`;
				#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--
			}
			
			if($aligner =~ "BOWTIE"){

                                #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                                #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                                my $cmd         = "$aligncommand1 $refGenome $path2fastq/$samplesInputs[$i]\.fastq -S $path2currentSampleDir/$samplesInputs[$i]\.sam";
                                `echo "$cmd" >> $QSUBint`;
                                #--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

                        }

		}

		#---------------------------------------------
		# Keep track of the jobs in @myJobs
#		my $jobName     = "$myJobName$i";
		my $jobName     = "$samplesInputs[$i]_$myJobName$i";
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
                $finalcmd	= "FINAL=\`$SUBcommand -N $iterateJobName -o $path2qsub -e $path2qsub -W depend=$SUBdependCondition\:$myJobsVec $IterateSH`";
        }
	 if($SUBkey =~ "bsub"){
		foreach( @myJobs ){ $_ = "$SUBdependCondition\(\"".$_."\")" ; }
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
	foreach my $i (0 .. $#samplesInputs) {

		# Prepare a personal qsub script
		my $QSUBint	= "$tmpscr/$myJobName/$samplesInputs[$i]\_$myJobName\.sh";
		`cp $scrhead $QSUBint`;
                `chmod 777 $QSUBint`;
				
		my $path2currentSampleDir	= "$path2aligned/$samplesInputs[$i]";
		my $j=0;
		
		# -----------------------------------------		
		# Get unique matches only
		#my $samplep = $samplesInputs[$i];
#		if( $makeunique && ($enforceisize == "0") ){

			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
			# print "\n\n Getting uniquely mapped reads for $path2aligned/$samplesInputs[$i]\.sam";
			my $cmd		= "grep -E '\\sX0:i:1\\s' $path2currentSampleDir/$samplesInputs[$i]\.sam > $path2currentSampleDir/$samplesInputs[$i]\.u.sam";
                        `echo "$cmd" >> $QSUBint`;
			#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

#		}			
		
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

		print "\n\n Samtools processing for $samplesInputs[$i] done. \n\n";
		
		#---------------------------------------------
		# Keep track of the jobs in @myJobs
#		my $jobName	= "$myJobName$i";
		my $jobName     = "$samplesInputs[$i]_$myJobName$i";
		push(@myJobs, "$jobName");
		if($SUBkey =~ "qsub"){  $cmd	= "$jobName=`$SUBcommand -o $path2qsub -e $path2qsub $QSUBint`";}
		if($SUBkey =~ "bsub"){  $cmd    = "$SUBcommand -J $jobName -o $path2qsub\/$jobName.out -e $path2qsub\/$jobName.err $QSUBint";}
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
                $finalcmd	= "FINAL=\`$SUBcommand -N $iterateJobName -o $path2qsub -e $path2qsub -W depend=$SUBdependCondition\:$myJobsVec $IterateSH`";
        }
	if($SUBkey =~ "bsub"){
                foreach( @myJobs ){ $_ = "$SUBdependCondition\(\"".$_."\")" ; }
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
# Running PeakCaller

if( $peakcalling =~ "TRUE" ){

	print "\n\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n";
	print "\n Run PeakCaller\n";
	
	#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	my $myJobName		= "peakcalling";
	my $iterateJobName	= "Iterate_$myJobName\_$expFolder";
	my $path2qsub		= "$tmpscr/$myJobName/$SUBkey";

	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        # Create folders
        unless( -d "$path2peakcalling" )    { `mkdir $path2peakcalling`; }
	
        unless( -d "$path2peakcalling/bigwig" )         { `mkdir $path2peakcalling/bigwig`; }
        unless( -d "$path2peakcalling/narrowPeak" )     { `mkdir $path2peakcalling/narrowPeak`; }
        unless( -d "$path2peakcalling/broadPeak" )	{ `mkdir $path2peakcalling/broadPeak`; }

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
	my @myJobsInputs;
	my @myJobsSamples;

	# Copy script and create folder
	`cp $path2CustFct/$peakcaller $tmpscr/`;

	#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
	foreach my $i (0 .. $#samples) {

		my $path2currentSampleDir	= "$path2aligned/$samples[$i]";
		my $path2currentInputDir	= "$path2aligned/$allinputs[$i]";	

		unless( -d "$path2peakcalling/$samples[$i]" )	{ `mkdir $path2peakcalling/$samples[$i]`; }
		unless( -d "$path2peakcalling/$allinputs[$i]" )	{ `mkdir $path2peakcalling/$allinputs[$i]`; }		

		#-----------------------------------------------------------
		# Prepare a personal qsub script
		my $QSUBint  = "$tmpscr/$myJobName/$samples[$i]\_$myJobName\.sh";
		`cp $scrhead $QSUBint`;
                `chmod 777 $QSUBint`;
		
		#-----------------------------------------------------------
                # Parameters
		my $path2currentSample	= "$path2currentSampleDir/$samples[$i]\.bam";
		my $path2currentInput	= "$path2currentInputDir/$allinputs[$i]\.bam";
		my $code		= "$tmpscr/$peakcaller";

		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+		
		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		my $cmd		= "Rscript $code $path2expFolder $path2currentSample $path2currentInput $fdr $posopt $densityopt &>> $path2qsub/$samples[$i]\_peakcalling.log";
		`echo "$cmd" >> $QSUBint`;
                #--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

		#---------------------------------------------
		# Keep track of the jobs in @myJobs
#		my $jobName	= "$myJobName$i";
		my $jobName     = "$samples[$i]_$myJobName$i";
		push(@myJobsSamples, "$jobName");
		if($SUBkey =~ "qsub"){  $cmd	= "$jobName=`$SUBcommand -o $path2qsub -e $path2qsub $QSUBint`";}
		if($SUBkey =~ "bsub"){  $cmd    = "$SUBcommand -J $jobName -o $path2qsub\/$jobName.out -e $path2qsub\/$jobName.err $QSUBint";}
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
	if($SUBkey =~ "qsub"){
                foreach( @myJobs ){ $_ = "\$".$_ ; }
                my $myJobsVec   = join(":", @myJobs);
                $finalcmd	= "FINAL=\`$SUBcommand -N $iterateJobName -o $path2qsub -e $path2qsub -W depend=$SUBdependCondition\:$myJobsVec $IterateSH`";
        }
	if($SUBkey =~ "bsub"){
                foreach( @myJobs ){ $_ = "$SUBdependCondition\(\"".$_."\")" ; }
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
# Clean up bigwig files

if( $cleanbigwig =~ "TRUE" ){

	print "\n\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n";
	print "\n Bigwig file cleaning\n";

	#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	my $myJobName		= "cleanbigwig";
	my $iterateJobName	= "Iterate_$myJobName\_$expFolder";
	my $myJobName2		= "wigToBigwig";
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

	my $QSUBintWig  = "$tmpscr/$myJobName/$myJobName2\.sh";
	open $QSUBintWig, ">", "$QSUBintWig" or die "Can't open '$QSUBintWig'";
	print $QSUBintWig "#!/bin/bash\n";
	close $QSUBintWig;
	`chmod 777 $QSUBintWig`;
	
	#*----------------------------------------------------------------------*
        # Open special character file and store in @specialchar
	
	print "\n These chromosomes will be deleted: @sc \n";

	#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>	
	foreach my $i (0 .. $#samples) {

		# Prepare a personal qsub script
		my $QSUBint 	= "$tmpscr/$myJobName/$samples[$i]\_$myJobName\.sh";
		`cp $scrhead $QSUBint`;
                `chmod 777 $QSUBint`;

		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		# Move narrow and broadPeaks to their resp. folder
		my $cmd         = "mv $path2peakcalling/$samples[$i]/$samples[$i]\.broadPeak $path2peakcalling/broadPeak/";
		`echo "$cmd" >> $QSUBint`;
		`chmod 777 $QSUBint`;
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
		`cp $scrhead $QSUBintWig`;
                `chmod 777 $QSUBintWig`;

		$cmd		= "wigToBigWig $path2peakcalling/$samples[$i]/$samples[$i]\.density$finalcount.wig $chrlens $path2peakcalling/bigwig/$samples[$i]\.bw";
		`echo "$cmd"  >> $QSUBintWig`;

		#--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--

		# Keep track of the jobs in @myJobs
		#---------------------------------------------
#		my $jobName	= "$myJobName$i";
		my $jobName     = "$samples[$i]_$myJobName$i";
		push(@myJobs, $jobName);
		if($SUBkey =~ "qsub"){  $cmd	= "$jobName=`$SUBcommand -o $path2qsub -e $path2qsub $QSUBint`";}
		if($SUBkey =~ "bsub"){  $cmd    = "$SUBcommand -J $jobName -o $path2qsub\/$jobName.out -e $path2qsub\/$jobName.err $QSUBint";}
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
	if($SUBkey =~ "qsub"){
                foreach( @myJobs ){ $_ = "\$".$_ ; }
                my $myJobsVec   = join(":", @myJobs);
                $finalcmd	= "FINAL=\`$SUBcommand -N $iterateJobName -o $path2qsub -e $path2qsub -W depend=$SUBdependCondition\:$myJobsVec $IterateSH`";
        }
	if($SUBkey =~ "bsub"){
                foreach( @myJobs ){ $_ = "$SUBdependCondition\(\"".$_."\")" ; }
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

	exit 1;

}



if($cleanfiles =~ "TRUE"){

	print "\n\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n";
        print "\n Cleaning files and folders\n";

	#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        my $myJobName           = "cleanfiles";
	my $iterateJobName	= "Iterate_$myJobName\_$expFolder";
        my $path2qsub           = "$tmpscr/$myJobName/$SUBkey";
	
	#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        # Create folders
        unless( -d "$path2aligned" )    { `mkdir $path2aligned`;        }
        unless( -d "$path2bam" )        { `mkdir $path2bam`;            }
        unless( -d "$path2bamRX" )	{ `mkdir $path2bamRX`;          }


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

	print "\n\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n";
        print "\n Moving .bam and .bai files from the Tophat to the aligned folder \n";

	#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
	foreach my $i (0 .. $#samplesInputs) {
	
		# Count reads in .bam and store in LibrarySize.txt
		my $path2currentSampleDir	= "$path2aligned/$samplesInputs[$i]";

		# Prepare a personal qsub script
		my $QSUBint	= "$tmpscr/$myJobName/$samplesInputs[$i]\_$myJobName\.sh";
		`cp $scrhead $QSUBint`;
                `chmod 777 $QSUBint`;

		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

		my $cmd		= "cp $path2aligned/$samplesInputs[$i]/$samplesInputs[$i].bam $path2bam/";
		`echo "$cmd" >> $QSUBint`;
		$cmd		= "cp $path2aligned/$samplesInputs[$i]/$samplesInputs[$i].bai $path2bam/";
		`echo "$cmd" >> $QSUBint`;
		# Check existence of _RX files. If they exists, move them
		if(-e "$path2aligned/$samplesInputs[$i]\_RX/$samplesInputs[$i].bam"){
			$cmd		= "cp $path2aligned/$samplesInputs[$i]\_RX/$samplesInputs[$i].bam $path2bamRX/";
			`echo "$cmd" >> $QSUBint`;
			$cmd		= "cp $path2aligned/$samplesInputs[$i]\_RX/$samplesInputs[$i].bai $path2bamRX/";        
			`echo "$cmd" >> $QSUBint`;
		}
		#---------------------------------------------
		# Keep track of the jobs in @myJobs
#		my $jobName     = "$myJobName$i";
		my $jobName     = "$samplesInputs[$i]_$myJobName$i";
		push(@myJobs, "$jobName");
		if($SUBkey =~ "qsub"){  $cmd	= "$jobName=`$SUBcommand -o $path2qsub -e $path2qsub $QSUBint`";}
		if($SUBkey =~ "bsub"){  $cmd    = "$SUBcommand -J $jobName -o $path2qsub\/$jobName.out -e $path2qsub\/$jobName.err $QSUBint";}
		open $QSUB, ">>", "$QSUB" or die "Can't open '$QSUB'";
		print $QSUB "$cmd\n";
		close $QSUB;

	}
	
	#-----------------------------------------------------
	# Copy .bam and .bai files to bam folder
	
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
                $finalcmd	= "FINAL=\`$SUBcommand -N $iterateJobName -o $path2qsub -e $path2qsub -W depend=$SUBdependCondition\:$myJobsVec $IterateSH`";
        }
	if($SUBkey =~ "bsub"){
                foreach( @myJobs ){ $_ = "$SUBdependCondition\(\"".$_."\")" ; }
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
# Normalize bam files using the chiprx genome reads

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
        unless( -d "$path2GRangesRX" )	{ `mkdir $path2GRangesRX`; }	

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
                my $QSUBint     = "$tmpscr/$myJobName/$myJobName\_qsub.sh";
                `cp $scrhead $QSUBint`;
                `chmod 777 $QSUBint`;

                #-----------------------------------------------------------
                # Parameters
                `cp $path2CustFct/Bam2GRangesRemote.R $tmpscr/`;
                my $code        = "$tmpscr/Bam2GRangesRemote.R";

                #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+          IMPORTANT CODE HERE         -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

                my $cmd         = "Rscript $code $path2expFolder $path2bam $path2GRanges $path2CustFct $wigBinSize &>> $path2qsub/GRanges.log";
                `echo "$cmd" >> $QSUBint`;
                #--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--


                #---------------------------------------------
                # Keep track of the jobs in @myJobs
                my $jobName     = "$myJobName$i";
                push(@myJobs, "$jobName");
                if($SUBkey =~ "qsub"){  $cmd    = "$jobName=`$SUBcommand -o $path2qsub -e $path2qsub $QSUBint`";}
                if($SUBkey =~ "bsub"){  $cmd    = "$SUBcommand -J $jobName -o $path2qsub\/$jobName.out -e $path2qsub\/$jobName.err $QSUBint";}
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
                $finalcmd	= "FINAL=\`$SUBcommand -N $iterateJobName -o $path2qsub -e $path2qsub -W depend=$SUBdependCondition\:$myJobsVec $IterateSH`";
        }
	if($SUBkey =~ "bsub"){
                foreach( @myJobs ){ $_ = "$SUBdependCondition\(\"".$_."\")" ; }
                my $myJobsVec   = join(" && ", @myJobs);
                $finalcmd	= "$SUBcommand -J $iterateJobName -o $path2qsub\/$iterateJobName.out -e $path2qsub\/$iterateJobName.err -w \'$myJobsVec\' $IterateSH";
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

           
                
                
                
                
                
                
                
                
                
                
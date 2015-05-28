#!/usr/bin/perl -w

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#************************************************************************
#*                                                                      *
#*                ChIPseq PERL script                                   *
#*                                                                      *
#************************************************************************

#************************************************************************
#*                                                                      *
#*	Open Targets.txt and start pipeline                             *
#*                                                                      *
#*----------------------------------------------------------------------*

if( $ARGV[0] ) { $path2expFolder = $ARGV[0]; }
else{ die "\n\n----------------------------------------\n\n Provide the path where to your project: </PATH/TO/PROJECT> \n\n--------------------------------------------------------------------------------\n\n"; }
        
#*----------------------------------------------------------------------*
# Read Targets.txt file

my $Targets = "$path2expFolder/DataStructure/Targets.txt";
open(INPUT, $Targets) || die "Error opening $Targets : $!\n\n\n";

my ($expFolder, $genome, $genomeRX, $userFolder, $path2ChIPseqScripts, $path2ChIPseq, $path2fastqgz)    = ("NA", "NA", "NA", "NA", "NA", "NA", "NA");
my ($unzip, $qc, $chiprx, $map, $filter, $peakcalling, $cleanbigwig, $cleanfiles, $granges)		= ("FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE");
my (@sc, @lines2remove)                                                                                 = ();
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
#               $path2NEAT = "$1";
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
                $aligner = "$1";
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
                if (grep /\bcleanfiles\b/i, $_ )        { $cleanfiles           = "TRUE"; push @steps2execute, "Cleanfiles";    }
                if (grep /\bgranges\b/i, $_ )           { $granges              = "TRUE"; push @steps2execute, "GRanges";	}
        }

} # end of Targets.txt



my $AdvSettings = "$path2expFolder/DataStructure/AdvancedSettings.txt";
open(INPUT, $AdvSettings) || die "Error opening $AdvSettings : $!\n\n\n";

my ($removepcr, $makeunique, $ndiff, $aligncommand1, $fdr, $posopt, $densityopt, $enforceisize)		= ("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA");

while(<INPUT>) {
	if (/# Ressource_manager/) {
                $_ =~ m/"(.+?)"/;
		if (grep /\bqsub/i, $_ )	{ $SUBkey	= "qsub"; $SUBheader	= "QSUB_header.sh";	}
		if (grep /\bbsub/i, $_ )      	{ $SUBkey    	= "bsub"; $SUBheader    = "BSUB_header.sh";	}
		$SUBcommand = "$1";
        }
	elsif (/# Unzip_comand/) {
                $_ =~ m/"(.+?)"/;
                $unzipCommand = "$1";
        }
	elsif (/# Zip_file_extension/) {
                $_ =~ m/"(.+?)"/;
                $zipExtension = "$1";
        }
	if (/# Bwa_maxEditDist/) {
                $_ =~ m/"(.+?)"/;
                $ndiff = "$1";
        }
	if (/# Align_command_line_1/) {
                $_ =~ m/"(.+?)"/;
                $aligncommand1 = "$1";
	}
	if (/# Filter_removePCRdup/) {
                $_ =~ m/"(.+?)"/;
                $removepcr = "$1";
        }
	if (/# Filter_makeUniqueRead/) {
                $_ =~ m/"(.+?)"/;
                $makeunique = "$1";
        }
	if (/# PeakCaller_fdr/) {
                $_ =~ m/"(.+?)"/;
                $fdr = "$1";
        }
	if (/# PeakCaller.posopt/) {
                $_ =~ m/"(.+?)"/;
                $posopt = "$1";
        }
	if (/# PeakCaller_densityopt/) {
                $_ =~ m/"(.+?)"/;
                $densityopt = "$1";
        }
	if (/# PeakCaller_enfSize/) {
                $_ =~ m/"(.+?)"/;
                $enforceisize = "$1";
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

# Remove samples that have "_R2" as these are the paired lanes of "_R1"
my @samplesPE;
my @samplesNoPE;
my @samples2unzip = @samples;
if( $PE ){

	print "\nPE experiment. \n";
        foreach my $i (0 .. $#samples) {
                if ( grep /\_R2$/, $samples[$i] ){
                        # print "\t $i\. '_R2' sample found. \t ($samples[$i]) \n";
                        push(@samplesPE, $samples[$i]);
                }
                else{
                     	# print "\t $i\. Main sample found.  \t ($samples[$i]) \n";
                        push(@samplesNoPE, $samples[$i]);
                }
        }

	@samples = @samplesNoPE;

}

#print "\n\n\norisamples:   @orisamples\n";
#print "samples2unzip:   @samples2unzip\n";
#print "samples: \t @samples\n";
#print "inputs: \t @inputs\n";
#print "samplesNoPE: \t @samplesNoPE\n";
#print "samplesPE: \t @samplesPE\n";

#*----------------------------------------------------------------------*
# Define paths
my $path2expFolder	= "$userFolder/$expFolder";
$Targets                = "$path2expFolder/DataStructure/Targets.txt";

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
print "\n path2fastq_gz:\t\t $path2fastqgz";
print "\n Targets:\t\t $path2expFolder/DataStructure/Targets.txt";
print "\n chrlens:\t\t $chrlens";
print "\n refGenome:\t\t $refGenome";
if($chiprx =~ "TRUE"){
	print "\n\t - chrlensRX:\t $chrlensRX";
	print "\n\t - refGenomeRX:\t $refGenomeRX";
}
print "\n";
print "\n Paired end sequencing:\t $PE";
print "\n Aligner algorithm:\t $aligner";
print "\n Remove pcr dupl:\t $removepcr";
print "\n Make unique reads:\t $makeunique";
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
# Create different folders if they do not exist
# subdir names

my $tmpscr		= "$path2expFolder/scripts";
my $scrhead		= "$path2ChIPseqScripts/$SUBheader";
my $path2iterate	= "$tmpscr/iterate";
my $path2qsub		= "$path2iterate/$SUBkey";
my $path2DataStructure	= "$path2expFolder/DataStructure";
my $path2chrlens        = "$path2DataStructure/$genome";

unless( -d "$tmpscr" )				{ `mkdir $tmpscr`;			}
unless( -d "$path2iterate" )			{ `mkdir $path2iterate`;		}
unless( -d "$path2qsub" )			{ `mkdir $path2qsub`;			}
unless( -d "$path2chrlens" )                    { `mkdir $path2chrlens`;                }

#------------------------------------------------------------------------
# Copy and create various files with execution permissions


# ------ Copy temp.sh file to $tmpscr

`cp $scrhead $tmpscr`;
my $nameOfChIPseqFile	= "ChIPpip";


# ------ ChIPseqMainCopy.pl to iterate later

`cp $path2ChIPseqScripts/$nameOfChIPseqFile\.pl $path2iterate/`;
`mv $path2iterate/$nameOfChIPseqFile\.pl $path2iterate/$nameOfChIPseqFile\_$expFolder\.pl`;


# ------ Copy chr_lens.dat file to $path2DataStructure

`cp $chrlens $path2chrlens`;
if($chiprx =~ "TRUE"){
	my $path2chrlensRX		= "$path2DataStructure/$genomeRX";
	unless( -d "$path2chrlensRX" )	{ `mkdir $path2chrlensRX`; }
	`cp $chrlensRX $path2chrlensRX`;
}


# ------ ChIPseqMainIterative.sh to iterate later

my $ChIPseqMainIterative = "$path2iterate/$nameOfChIPseqFile\_$expFolder\.sh";
`cp $scrhead $ChIPseqMainIterative`;
open $ChIPseqMainIterative, ">>", "$ChIPseqMainIterative" or die "Can't open '$ChIPseqMainIterative'\n";
#print $ChIPseqMainIterative "`echo \"perl $path2iterate\/$nameOfChIPseqFile\_$expFolder\.pl $path2expFolder\"`";
print $ChIPseqMainIterative "perl $path2iterate\/$nameOfChIPseqFile\_$expFolder\.pl $path2expFolder";
close $ChIPseqMainIterative;
# Change permissions of file so that it can be executed later on
`chmod 777 $ChIPseqMainIterative`;


#*----------------------------------------------------------------------*
# Prepar file containing the jobs to run

# Add the first iteration of the script to $SubmitJobsToCluster
my $firstcmd    = NULL;
if($SUBkey =~ "qsub"){	$firstcmd	= "FIRST=`$SUBcommand -o $path2qsub -e $path2qsub $ChIPseqMainIterative`";}
if($SUBkey =~ "bsub"){  $firstcmd	= "$SUBcommand -J Iterate_$expFolder -o $path2qsub\/Iterate_$expFolder.out -e $path2qsub\/Iterate_$expFolder.err $ChIPseqMainIterative";}


my $IterateSH	= "$path2iterate/Iterate_$expFolder.sh";
open $IterateSH, ">", "$IterateSH" or die "Can't open '$IterateSH'";
print $IterateSH "#!/bin/bash\n";
print $IterateSH "$firstcmd\n";
close $IterateSH;
`chmod 777 $IterateSH`;


#*----------------------------------------------------------------------*
# Submit jobs to run 

#print "\n\n--------------------------------------------------------------------------------------------------\n";
#print "\n  Submitting job to $SUBkey cluster: \t `sh $IterateSH` \n";
`sh $IterateSH`;

#*----------------------------------------------------------------------*
# Exit script

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "\n Exiting INITIAL section with no known error \n";
print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n";

exit 0;


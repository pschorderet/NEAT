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

my ($expFolder, $genome, $userFolder, $path2ChIPseqScripts, $path2ChIPseq, $path2fastqgz)	= ("NA", "NA", "NA", "NA", "NA", "NA");
my ($unzip, $qc, $map, $filter, $peakcalling, $cleanbigwig, $cleanfolders)                      = ("FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE");
my (@sc, @lines2remove)                                                                         = ();
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
                if (grep /\bunzip\b/i, $_ )             { $unzip                = "TRUE"; push @steps2execute, "Unzip";         }
                if (grep /\bqc\b/i, $_ )                { $qc                   = "TRUE"; push @steps2execute, "QC";            }
                if (grep /\bmap\b/i, $_ )               { $map                  = "TRUE"; push @steps2execute, "Map";           }
                if (grep /\bfilter\b/i, $_ )            { $filter               = "TRUE"; push @steps2execute, "Filter";        }
                if (grep /\bpeakcalling\b/i, $_ )	{ $peakcalling          = "TRUE"; push @steps2execute, "Peakcalling";   }
                if (grep /\bcleanbigwig\b/i, $_ )	{ $cleanbigwig          = "TRUE"; push @steps2execute, "Cleanbigwig";   }
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
	if (/# PeakCaller.enfSize/) {
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
        $line =~ /^[# = " OriFileName FileName OriInpName InpName]/ and next;
        push(@orisamples, $line);
}
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
# Define paths

my $path2expFolder 	= "$userFolder/$expFolder";
$Targets 		= "$path2expFolder/DataStructure/Targets.txt";

#*----------------------------------------------------------------------*

chdir "$path2expFolder";
`clear`;
print "\n\n\n#################################################################################################";
print "\n#												#";
print "\n#				ChIPseq pipeline v1.0.1 (Jan 2015)				#";
print "\n#												#";
print "\n#################################################################################################\n";
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
#print "\n Current working dir:\t $path2expFolder";
#print "\n";
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
my $scrhead		= "$path2ChIPseqScripts/QSUB_header.sh";
my $path2iterate	= "$tmpscr/iterate";
my $path2qsub		= "$path2iterate/qsub";
my $path2DataStructure	= "$path2expFolder/DataStructure";
my $path2aligned	= "$path2expFolder/aligned";
my $path2fastq		= "$path2expFolder/fastq";

unless( -d "$path2fastq" )			{ `mkdir $path2fastq`;			}
unless( -d "$path2aligned" )			{ `mkdir $path2aligned`;		}
unless( -d "$path2expFolder/peakcalling" )	{ `mkdir $path2expFolder/peakcalling`;	}
unless( -d "$tmpscr" )				{ `mkdir $tmpscr`;			}
unless( -d "$path2iterate" )			{ `mkdir $path2iterate`;		}
unless( -d "$path2qsub" )			{ `mkdir $path2qsub`;			}


#------------------------------------------------------------------------
# Copy and create various files with execution permissions


# ------ Copy temp.sh file to $tmpscr
`cp $scrhead $tmpscr`;
my $nameOfChIPseqFile	= "ChIPpip";

# ------ ChIPseqMainCopy.pl to iterate later
`cp $path2ChIPseqScripts/$nameOfChIPseqFile\.pl $path2iterate/`;

# ------ Copy chr_lens.dat file to $path2DataStructure
`cp $chrlens $path2DataStructure`;

# ------ ChIPseqMainIterative.sh to iterate later
my $ChIPseqMainIterative = "$path2iterate/$nameOfChIPseqFile\.sh";
`cp $scrhead $ChIPseqMainIterative`;
open $ChIPseqMainIterative, ">>", "$ChIPseqMainIterative" or die "Can't open '$ChIPseqMainIterative'\n";
print $ChIPseqMainIterative "`echo \"perl $path2iterate\/$nameOfChIPseqFile\.pl $path2expFolder\"`";
close $ChIPseqMainIterative;
# Change permissions of file so that it can be executed later on
`chmod 777 $ChIPseqMainIterative`;


#*----------------------------------------------------------------------*
# Prepar file containing the jobs to run

# Add the first iteration of the script to $SubmitJobsToCluster
my $firstcmd    = "FIRST=`qsub -N Iterate -o $path2qsub -e $path2qsub $ChIPseqMainIterative`";
my $IterateSH	= "$path2iterate/IterateSH.sh";
open $IterateSH, ">", "$IterateSH" or die "Can't open '$IterateSH'";
print $IterateSH "#!/bin/bash\n";
print $IterateSH "$firstcmd\n";
close $IterateSH;
`chmod 777 $IterateSH`;


#*----------------------------------------------------------------------*
# Submit jobs to run 

#print "\n\n--------------------------------------------------------------------------------------------------\n";
#print "\n  Submitting job to cluster: \t `sh $IterateSH` \n";
`sh $IterateSH`;

#*----------------------------------------------------------------------*
# Exit script

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "\n Exiting INITIAL section with no known error \n";
print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n";

exit 1;


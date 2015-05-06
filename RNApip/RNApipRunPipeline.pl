#!/usr/bin/perl -w

#my $phred="phred64";

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#************************************************************************
#*                                                                      *
#*                RNApip PERL script                                    *
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

my ($expFolder, $genome, $userFolder, $path2RNAseqScripts, $path2RNAseq, $path2fastqgz, $chrlens, $path2gtfFile)				= ("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA");
my ($unzip, $qc, $map, $filter, $cleanfiles, $granges)												= ("FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE");

my (@sc, @lines2remove)                                                                                                                         = ();
# Find paths to different folders in the Targets.txt file
while(<INPUT>) {
        if (/# My_personal_email/) {
                $_ =~ m/"(.+?)"/;
                $email = "$1";
        }
	elsif (/# My_project_title/) {
                $_ =~ m/"(.+?)"/;
                $expFolder = "$1";
        }
	elsif (/# Reference_genome/) {
                $_ =~ m/"(.+?)"/;
                $genome = "$1";
        }
	elsif (/# Remote_path_to_proj/) {
                $_ =~ m/"(.+?)"/;
                $userFolder = "$1";
        }
	elsif (/# Remote_path_to_NEAT/) {
                $_ =~ m/"(.+?)"/;
                $path2RNAseq = "$1\/RNApip";
                $path2RNAseqScripts = join("", $path2RNAseq, "/scripts");
        }
	elsif (/# Remote_path_to_orifastq_gz/) {
                $_ =~ m/"(.+?)"/;
                $path2fastqgz = "$1";
        }
	elsif (/# Remote_path_to_chrLens_dat/) {
                $_ =~ m/"(.+?)"/;
                $chrlens = "$1";
        }
	elsif (/# Remote_path_to_RefGen_fasta/) {
                $_ =~ m/"(.+?)"/;
                $refGenome = "$1";
        }
	elsif (/# Remote_path_to_ann_gtf_file/) {
                $_ =~ m/"(.+?)"/;
                $path2gtfFile = "$1";
        }
	elsif (/# Aligner_algo_short/) {
                $_ =~ m/"(.+?)"/;
                $aligner = "$1";
        }
	elsif (/# Paired_end_seq_run/) {
                $_ =~ m/"(.+?)"/;
                $PE = "$1";
        }
	elsif (/# Steps_to_execute_pipe/) {
                $_ =~ m/"(.+?)"/;
                @steps2execute = ();
                if (grep /\bunzip\b/i, $_ )             { $unzip                = "TRUE"; push @steps2execute, "Unzip";         }
                if (grep /\bqc\b/i, $_ )                { $qc                   = "TRUE"; push @steps2execute, "QC";            }
                if (grep /\bmap\b/i, $_ )               { $map                  = "TRUE"; push @steps2execute, "Map";           }
                if (grep /\bfilter\b/i, $_ )            { $filter               = "TRUE"; push @steps2execute, "Filter";        }
		if (grep /\bcleanfiles\b/i, $_ )        { $cleanfiles		= "TRUE"; push @steps2execute, "Cleanfiles";    }
                if (grep /\bgranges\b/i, $_ )           { $granges		= "TRUE"; push @steps2execute, "GRanges";	}

        }

} # end of Targets.txt



my $AdvSettings = "$path2expFolder/DataStructure/AdvancedSettings.txt";
open(INPUT, $AdvSettings) || die "Error opening $AdvSettings : $!\n\n\n";

my ($removepcr, $makeunique, $ndiff, $aligncommand)				= ("NA", "NA", "NA", "NA");

while(<INPUT>) {
	if (/# Ressource_manager/) {
                $_ =~ m/"(.+?)"/;
                if (grep /\bqsub/i, $_ )        { $SUBkey	= "qsub";}
                if (grep /\bbsub/i, $_ )        { $SUBkey	= "bsub";}
                $SUBcommand = "$1";
        }
	if (/# Unzip_comand/) {
                $_ =~ m/"(.+?)"/;
                $unzipCommand = "$1";
        }
	elsif (/# Zip_file_extension/) {
                $_ =~ m/"(.+?)"/;
                $zipExtension = "$1";
        }
        elsif (/# Filter_removePCRdup/) {
                $_ =~ m/"(.+?)"/;
                $removepcr = "$1";
        }
	elsif (/# Filter_makeUniqueRead/) {
                $_ =~ m/"(.+?)"/;
                $makeunique = "$1";
        }
	elsif (/# Filter_maxEditDist/) {
                $_ =~ m/"(.+?)"/;
                $ndiff = "$1";
        }
	elsif (/# Align_command_opt/) {
                $_ =~ m/"(.+?)"/;
                $aligncommand = "$1";
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
%seen		= ();
@samples	= grep { ! $seen{$_} ++ } @samples;

# Remove samples that have "_R2" as these are the paired lanes of "_R1"
my @samples2unzip	= @samples;
my @samplesPE;
my @samplesNoPE;
if( $PE ){

	print "\nPE experiment. \n";
	foreach my $i (0 .. $#samples) {
        	if ( grep /\_R2$/, $samples[$i] ){ 	
			print "\t $i\. '_R2' sample found. \t ($samples[$i]) \n";
			push(@samplesPE, $samples[$i]);
		}
		else{	
			print "\t $i\. Main sample found.  \t ($samples[$i]) \n";
			push(@samplesNoPE, $samples[$i]);
		}
	}
	@samples = @samplesNoPE;
}

#print "\n\n\norisamples:   @orisamples\n";
print "samples2unzip:   @samples2unzip\n";
print "samples: \t @samples\n";
#print "samplesNoPE: \t @samplesNoPE\n";
print "samplesPE: \t @samplesPE\n";

#*----------------------------------------------------------------------*
# Define paths

my $path2expFolder	= "$userFolder/$expFolder";
$Targets		= "$path2expFolder/DataStructure/Targets.txt";

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
#print "\n Current working dir:\t $path2expFolder";
#print "\n";
print "\n .........................................";
print "\n Performing following modules:";
print "\n .........................................";
print "\n unzip:\t\t\t $unzip \t ($unzipCommand filename.fastq$zipExtension)";
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
#print "\n----------------------------------------\n";



#*----------------------------------------------------------------------*
# Create different folders if they do not exist
# subdir names

my $tmpscr		= "$path2expFolder/scripts";
my $scrhead		= "$path2RNAseqScripts/QSUB_header.sh";
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
my $nameOfRNAseqFile	= "RNApip";


# ------ RNAseqMainCopy.pl to iterate later

`cp $path2RNAseqScripts/$nameOfRNAseqFile\.pl $path2iterate/`;
`mv $path2iterate/$nameOfRNAseqFile\.pl $path2iterate/$nameOfRNAseqFile\_$expFolder\.pl`;


# ------ Copy chr_lens.dat file to $path2DataStructure

`cp $chrlens $path2chrlens`;


# ------ RNAseqMainIterative.sh to iterate later

my $RNAseqMainIterative = "$path2iterate/$nameOfRNAseqFile\_$expFolder\.sh";
`cp $scrhead $RNAseqMainIterative`;
open $RNAseqMainIterative, ">>", "$RNAseqMainIterative" or die "Can't open '$RNAseqMainIterative'\n";
print $RNAseqMainIterative "`echo \"perl $path2iterate\/$nameOfRNAseqFile\_$expFolder\.pl $path2expFolder\"`";
close $RNAseqMainIterative;
# Change permissions of file so that it can be executed later on
`chmod 777 $RNAseqMainIterative`;


#*----------------------------------------------------------------------*
# Prepar file containing the jobs to run

# Add the first iteration of the script to $SubmitJobsToCluster
my $firstcmd    = NULL;
if($SUBkey =~ "qsub"){	$firstcmd	= "FIRST=`$SUBcommand -N Iterate_$expFolder -o $path2qsub -e $path2qsub $RNAseqMainIterative`";}
if($SUBkey =~ "bsub"){	$firstcmd	= "FIRST=`$SUBcommand -J Iterate_$expFolder -o $path2qsub -e $path2qsub $RNAseqMainIterative`";}
my $IterateSH	= "$path2iterate/Iterate\_$expFolder\.sh";
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

exit 1;


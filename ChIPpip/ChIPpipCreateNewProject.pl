#!/usr/bin/perl -w
use Cwd;
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#************************************************************************
#*                                                                      *
#*			ChIPseqCreateNewProject				*
#*                                                                      *
#************************************************************************

#*----------------------------------------------------------------------#
#                                                                       #
#               Create new arborescence for a ChIPseq project           #
#                                                                       #
#*----------------------------------------------------------------------*

if( $ARGV[0] ){ $path2NewMainFolder = $ARGV[0]; }
else{ die "\n\n--------------------------------------------------------------------------------\n \n Provide the path where you want to create your new project: </PATH/TO/NEWPROJECT> \n\n-------------------------------------------------------------------------------- \n\n\n"; }

print "\n\n *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- \n";
print "\n Creating new ChIPpip project \n";
print "\n -------------------------------------------------------------------------------------------------- \n";

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
`cp -r ./scripts/NewChIPpipProject/ $path2NewMainFolder`;
#
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

#------------------------------------------------------------
# Start filling in the Targets.txt file
#
my $Targets		= "$path2NewMainFolder/DataStructure/Targets.txt";

#-----------------------------------------
#------------------Title------------------
#
@path			= split("/", $path2NewMainFolder);
my $projtitle		= pop @path;
my $toreplace		= "<MAIN_TITLE>";
`/usr/bin/perl -p -i -e "s/$toreplace/$projtitle/" $Targets`;

#-----------------------------------------
#-------------path2expFolder--------------
#
$path2expFolder		= join("\\/",@path);
$toreplace		= "<PATH_TO_PROJECTFOLDER>";
`/usr/bin/perl -p -i -e "s/$toreplace/$path2expFolder/" $Targets`;

#-----------------------------------------
#-----------path2ChIPseqscripts-----------
#
$replace		= cwd;
@path			= split("/", $replace);
$replace		= join("\\/",@path);
$toreplace		= "<PATH_TO_CHIPPIP>";
`/usr/bin/perl -p -i -e "s/$toreplace/$replace/gi" $Targets`;

#------------------------------------------------------------
# Print for user to see

print "\n New ChIPpip project has been built as follows: \n";
print "\n\t $path2NewMainFolder";
print "\n\t\t └─ DataStructure";
print "\n\t\t\t └─ Targets.txt";
print "\n\t\t └─ fastq";
print "\n\t\t └─ scripts";
print "\n\n -------------------------------------------------------------------------------------------------- \n";
print "\n\n IMPORTANT: \t Fill in the Targets.txt file before running ChIPpipRunPipeline.pl ";
print "\n\t\t To modify Targets.txt, copy paste the following command in your terminal:";
print "\n\n\t\t nano $path2NewMainFolder/DataStructure/Targets.txt ";
print "\n";
`ls $path2NewMainFolder`;
print "\n\n\n *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- \n\n\n";

exit 0;




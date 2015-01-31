## RNApip: RNAseq pipelines made easy

RNApip is a perl/R package that supports the analysis of next generation sequencing (NGS) data as part of NEAT (NGS easy analysis toolkit). It is a versatile and easily configurable tool that allows users to go from compressed fastq files (.fastq.gz) to filter .bam files using a single command line. One central feature of RNApip is its ability to perform various tasks on complex sample setups while managing batch submissions and cluster queuing. RNApip can easily be implemented in any institution with limited to no programming knowledge. The workflow has been designed to efficiently run on a computer cluster using a distributed resource manager such as TORQUE. RNApip has been developed by and for wet-lab scientists as well as bioinformaticiens to ensure user-friendliness, management of complicated experimental setups and reproducibility in the big data era. To start using RNApip, please follow the tutorial.



### What does RNApip do
RNApip can run the following tasks using a single command line:

[ 1 ]       Unzip and rename fastq.gz files

[ 2 ]       Quality control of sequencing reads

[ 3 ]       Map reads (Tophat)

[ 4 ]       Filter reads

[ 5 ]       Email notification when pipeline has finished



### After RNApip
Once RNApip has been run, users are literally two clicks away from differnetial gene expression calling, using RNAmE. RNAmE is part of the NEAT toolkit and has been developed as a downstream module for RNApip. RNAmE supports various steps in the process between generating .bam files to obtaining readable data for wet-lab scientists including .pdf graphs (smear plots, etc), venn diagrams (gene expression overlap), count tables and consolidated RPKM tables. 


### Tutorial
A tutorial can be found in the Vignette folder. This accompanies new users through the entire analysis of a test dataset provided as part of the package.
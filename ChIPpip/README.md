## ChIPpip: ChIPseq pipelines made easy

ChIPpip is a perl/R package that supports the analysis of next generation sequencing (NGS) data as part of NEAT (NGS easy analysis toolkit). It is a versatile and easily configurable tool that allows users to go from compressed fastq files (.fastq.gz) to bigwigs using a single command line. One central feature of ChIPpip is its ability to perform various tasks on complex sample setups while managing batch submissions and cluster queuing. ChIPpip can easily be implemented in any institution with limited to no programming knowledge. The workflow has been designed to efficiently run on a computer cluster using a distributed resource manager such as TORQUE. ChIPpip has been developed by and for wet-lab scientists as well as bioinformaticiens to ensure user-friendliness, management of complicated experimental setups and reproducibility in the big data era. To start using ChIPpip, please read the README file.



### What does ChIPpip do
ChIPpip can run the following tasks using a single command line:

[ 1 ]       Unzip and rename fastq.gz files

[ 2 ]       Quality control of sequencing reads

[ 3 ]       Map reads (bwa)

[ 4 ]       Filter reads

[ 5 ]       Peakcalling (SPP)

[ 6 ]       Clean bigwig files

[ 7 ]       Email notification when pipeline has finished



### After ChIPpip
Once ChIPpip has been run, users are literally two clicks away from metagene analysis using ChIPmE. ChIPmE is part of the NEAT toolkit and has been developed as a downstream module for ChIPpip. ChIpmE supports various steps in the process between generating .bam files to obtaining readable data for wet-lab scientists including .pdf graphs (enrichments over features), venn diagrams (overlap of peaks) and count tables. 


### Tutorial
A tutorial can be found in the Vignette folder. This accompanies new users through the entire analysis of a test dataset provided as part of the package.
## NEAT: NGS pipelines made easy

NEAT is a next generation analysis toolkit that supports the analysis of large data including metagene analysis (ChIPseq) and differential gene expression analysis (RNAseq). 
NEAT can be run either on a cluster via the command line or directly through the respective applications. This allows users to generate metagene analysis (for ChIPseq data) as well as differentially expressed gene analysis (for RNAseq data) using a simple and reproducable pipeline. All files including count tables, RPKM values, DEG, venn diagrams, feature-centered enrichment plots, smear plots, etc are automatically saved and archived for a userfriendly and reproducable flow. 
The NEAT toolkit can easily be implemented in any institution with limited to no programming knowledge.
NEAT has been developped in collaboration with wet-lab scientists as well as bioinformaticiens to insure user-friendliness, management of complicated experimental setups and reproducibility in the big data era.
To start using NEAT, please read the README files corresponding to your experiemental setup found in the RNApip and ChIPpip folders.


### What NEAT does
NEAT has been developed as two main modules each conatining two sub-modules:

[ 1 ]       Build a pipeline (ChIPpip or RNApip)

    (i)     Unzip and rename

    (ii)    QC

    (iii)   Map

    (iv)    FIlter

    (v)     Peakcalling

    (vi)    Bigwig

[ 2 ]       Analyze a project (ChIPmE or RNAmE)

    (i)     Create human-readable data (pdf, xls sheets, etc)

### Tutorial
Novice users can follow through the entire analysis of a test dataset provided as part of the package. Please follow the tutorials found in the 'vignette' folder


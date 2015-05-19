## NEAT: NGS pipelines made easy

NEAT is a next generation analysis toolkit that supports the analysis of large data including metagene analysis (ChIPseq) and differential gene expression analysis (RNAseq). 

NEAT can be run either on a cluster (qsub or bsub) via the command line or directly through the respective applications. This allows users to generate metagene analysis (for ChIPseq data) as well as differentially expressed gene analysis (for RNAseq data) using a simple and reproducable pipeline without even needing to connect to a remote cluster via the terminal. All files including count tables, RPKM values, DEG, venn diagrams, feature-centered enrichment plots, smear plots, etc are automatically saved and archived for a user-friendly and reproducible flow. 

The NEAT toolkit can easily be implemented in any institution with limited to no programming knowledge.
NEAT has been developed in collaboration with wet-lab scientists as well as bioinformaticiens to insure user-friendliness, management of complicated experimental setups and reproducibility in the big data era.

To start using NEAT, please follow the tutorials found in the 'Vignette' folder


### What NEAT can do
NEAT has been developed as four main modules:

[ 1 ]       Create a new project

[ 2 ]       Build a pipeline

    (i)     Unzip and rename

    (ii)    QC
    
    (iii)	ChIP-rx normalization

    (iv)	Map

    (v)		Filter

    (vi)	Peakcalling

    (vii)	Bigwig
    
    (viii)	Wig
    
    (viii)	GRanges

[ 3 ]       Transfer project to local computer

[ 4 ]       Analyze project

    (i)     Create human-readable data (pdf graphes, xls sheets, etc)


### Tutorial
Novice users can follow through the entire analysis of a test dataset provided as part of the package. Please follow the tutorials found in the 'vignette' folder



### Devel mode
NEAT is still under developent. If you encounter any bug, please contact us.


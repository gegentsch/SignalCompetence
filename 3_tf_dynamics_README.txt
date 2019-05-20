System requirements
-------------------

Running of script requires R. It can be downloaded from one of the many Comprehensive R Archive Network (CRAN) mirror websites of the R Project for Statistical Computing. 

https://cran.r-project.org/mirrors.html

The script was tested with R version 3.4.1 on a Darwin operating system supporting 64-bit x86-64 variant of the Intel x86 processors.


Running R script
----------------

Before running this script, align ChIP-Seq reads (FASTQ) to the genome assembly 7.1 of X. tropicalis. 
FASTQ files can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113186 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48560.

Add path/BAM files to the 'DiffBind' input (.csv) files provided under ./xenTro71/supplementary_files
myod.csvpol2.csvsox3_ap.csvtf_pol2.csv

This script can be run from the command line:

R CMD batch 3_tf_dynamics.R


or with R studio or the R GUI:

source("3_tf_dynamics.R")


Expected output
---------------

This script visualises the dynamic recruitment of various transcription factors (TFs) to the genome of early X. tropicalis embryos. 

The following files will be produced upon run execution.

:: tf_motifs_enrichment.pdf, tf_motifs_peakVsBgd.pdf, tf_motifs_pval.pdf, myod_motifs_enrichment.pdf, myod_motifs_peakVsBgd.pdf, myod_motifs_pval.pdf, sox3_ap_motifs_enrichment.pdf, sox3_ap_motifs_peakVsBgd.pdf and sox3_ap_motifs_pval.pdf (Figures 3g, 4e and 5c and Supplementary Figures 5b and 7c)
Heatmap (enrichment, significance) and bubble plots (peak versus background) of DNA motif occurences for
(1) tf: TFs and signal mediators
(2) myod: Sox3 and RNAPII in the presence or absence of MyoD-HA
(3) sox3_ap: Sox3 in the head, trunk and tailbud

:: tf_pca.pdf, pol2_tf_pca.pdf, myod_pca.pdf and sox3_ap_pca.pdf (Figures 3d,e and 4b)
PCA plot of chromatin binding
(1) tf: TFs and signal mediators
(2) pol2_tf: TFs, signal mediators and RNAPII
(3) myod: Sox3 and RNAPII in the presence or absence of MyoD-HA
(4) sox3_ap: Sox3 in the head, trunk and tailbud

:: tf_spearman.pdf, tf_spearman_noCellNote.pdf and myod_spearman.pdf, myod_spearman_noCellNote.pdf, sox3_ap_spearman.pdf, sox3_ap_spearman_noCellNote.pdf (Supplementary Figures 5a and 7a,b)
Spearman correlation heatmaps of chromatin engagement (with or without cell note)
(1) tf: TFs and signal mediators
(2) myod: Sox3 and RNAPII in the presence or absence of MyoD-HA
(3) sox3_ap: Sox3 in the head, trunk and tailbud

:: myoD_heatmap.png and sox3_ap_heatmap.png (Figure 4c and Supplementary Figure 7d)
Ectopic expression of MyoD: Heatmap of DNA occupancy, significance of differential DNA occupancy (p-value) and DNA motif occurence

:: ColourBar_motif_pval.pdf, ColourBar_DNAocc_scaled.pdf and ColourBar_diffPval.pdf
Colour bar for DNA motif p-values, scaled DNA occupancy and significance of differential DNA occupancy 


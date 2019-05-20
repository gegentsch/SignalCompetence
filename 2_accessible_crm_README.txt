System requirements
-------------------

Running of script requires R. It can be downloaded from one of the many Comprehensive R Archive Network (CRAN) mirror websites of the R Project for Statistical Computing. 

https://cran.r-project.org/mirrors.html

The script was tested with R version 3.4.1 on a Darwin operating system supporting 64-bit x86-64 variant of the Intel x86 processors.

It took ~70 sec to complete execution of the R script from the command line (see below).


Running R script
----------------

This script can be run from the command line:
R CMD batch 2_accessible_crm.R


or with R studio or the R GUI:
source("2_accessible_crm.R")


Expected output
---------------

This script quantifies the DNA motif enrichment at accessible (conserved versus non-conserved) cis-regulatory modules (CRMs) at the mid-blastula transition (MBT). 

The following files will be produced upon run execution.

:: dhs_h3k4me1_pol2_conservation.png (Figure 1f)
Heatmap of chromatin accessibility, H3K4me1 and RNAPII occupancy and PhastCons score separated by conservation threshold of PhastCons >= 0.4 (conserved, top; non-conserved, bottom)

:: MBT_consVsNonCons_accessible_crm_motifs.pdf (Figure 2b)
Bubble plot shows the occurence of the DNA motifs in accessible CRMs (versus background) in conserved (cons) and non-conserved (ncons).

:: legend_dhs_scaled.pdf, legend_h3k_scaled.pdf and legend_pol2_scaled.pdf
Colour bars for chromatin accessiblity, H3K4me1 and RNAPII occupancy, respectively.


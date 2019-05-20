System requirements
-------------------

Running of script requires R. It can be downloaded from one of the many Comprehensive R Archive Network (CRAN) mirror websites of the R Project for Statistical Computing. 

https://cran.r-project.org/mirrors.html

The script was tested with R version 3.4.1 on a Darwin operating system supporting 64-bit x86-64 variant of the Intel x86 processors. 

It took ~20 min to complete execution of the R the script from the command line (see below).


Running R script
----------------

This script can be run from the command line:
R CMD batch 4_zga_lof.R


or with R studio or the R GUI:
source("4_zga_lof.R")


Expected output
---------------

This script covers the transcriptional analysis of zygotic genome activation (ZGA) under various loss-of-functions (LOFs). 
The following files will be produced upon run execution.

:: amanitin_mz_tpm.pdf (Figure 6e)
Dotplot: reduced ZGA under the influence of alfa-amanitin (% relative to control vs. total RNA level (TPM) around MBT/early gastula stage
Piechart of downregulated zygotic only and maternal-zygotic genes
Fold change comparison between mPouV LOF vs control and mPouV/Sox3 LOF vs control (with and without gene name labels)

:: mPSLOF_vs_BmpLOF_spatioAnVgZoom_noLab.pdf, mPSLOF_vs_BmpLOF_spatioAnVgZoom.pdf
:: mPSLOF_vs_BmpLOF_spatioDVZoom_noLab.pdf, mPSLOF_vs_BmpLOF_spatioDVZoom.pdf
:: mVegTLOF_vs_BmpLOF_spatioDVZoom_noLab.pdf, mVegTLOF_vs_BmpLOF_spatioDVZoom.pdf
Comparison of relative expression levels (percentage of control levels) between indicated LOF. Gene dots are coloured according to its expression along the animal-vegetal (AnVg) or dorso-ventral (DV) axis (with and without gene name labels).

Bargraph showing the level of gene misregulation under indicated LOFs

:: GO_term_LOF_ZGA.pdf (Supplementary Figure 11b)
Table of LOF effects on ZGA (restricted to alfa-amanitin downregulated (<= 50%, FDR <= 10%) genes (Supplementary Table 9)

:: spatial_down_ALL_LOF.txt, spatial_down_amanitin.txt, spatial_down_Bmp_LOF.txt, spatial_down_cMO.txt, spatial_down_mPouV_LOF.txt, spatial_down_mPS_LOF.txt, spatial_down_mTF_LOF.txt, spatial_down_mVegT_LOF.txt, spatial_down_Nodal_LOF.txt, spatial_down_Signals_LOF.txt and spatial_down_Wnt_LOF.txt
Text files the numbers/percentage of genes downregulated by indicated LOFs

:: AnVg_log2_colorBar.pdf, DV_log2_colorBar.pdf and FC_ratios_bar.pdf

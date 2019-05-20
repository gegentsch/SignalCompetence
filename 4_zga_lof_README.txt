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
Dotplot: reduced ZGA under the influence of alfa-amanitin (% relative to control vs. total RNA level (TPM) around MBT/early gastula stage:: zga_ama50.pdf (Figure 6e)
Piechart of downregulated zygotic only and maternal-zygotic genes:: mPouVkd_vs_mPSkd_FCLog.pdf and mPouVLOF_vs_mPSLOF_FCLog_noLab.pdf (Figure 6f)
Fold change comparison between mPouV LOF vs control and mPouV/Sox3 LOF vs control (with and without gene name labels)

:: mPSLOF_vs_BmpLOF_spatioAnVgZoom_noLab.pdf, mPSLOF_vs_BmpLOF_spatioAnVgZoom.pdf
:: mPSLOF_vs_BmpLOF_spatioDVZoom_noLab.pdf, mPSLOF_vs_BmpLOF_spatioDVZoom.pdf:: mPSLOF_vs_NodalLOF_spatioAnVgZoom_noLab.pdf, mPSLOF_vs_NodalLOF_spatioAnVgZoom.pdf:: mPSLOF_vs_NodalLOF_spatioDVZoom_noLab.pdf, mPSLOF_vs_NodalLOF_spatioDVZoom.pdf:: mPSLOF_vs_WntLOF_spatioAnVgZoom_noLab.pdf, mPSLOF_vs_WntLOF_spatioAnVgZoom.pdf:: mPSLOF_vs_WntLOF_spatioDVZoom_noLab.pdf, mPSLOF_vs_WntLOF_spatioDVZoom.pdf:: mVegTLOF_vs_BmpLOF_spatioAnVgZoom_noLab.pdf, mVegTLOF_vs_BmpLOF_spatioAnVgZoom.pdf
:: mVegTLOF_vs_BmpLOF_spatioDVZoom_noLab.pdf, mVegTLOF_vs_BmpLOF_spatioDVZoom.pdf:: mVegTLOF_vs_NodalLOF_spatioAnVgZoom_noLab.pdf, mVegTLOF_vs_NodalLOF_spatioAnVgZoom.pdf:: mVegTLOF_vs_NodalLOF_spatioDVZoom_noLab.pdf, mVegTLOF_vs_NodalLOF_spatioDVZoom.pdf:: mVegTLOF_vs_WntLOF_spatioAnVgZoom_noLab.pdf, mVegTLOF_vs_WntLOF_spatioAnVgZoom.pdf:: mVegTLOF_vs_WntLOF_spatioDVZoom_noLab.pdf and mVegTLOF_vs_WntLOF_spatioDVZoom.pdf (Figure 7a and Supplementary Figure 9)
Comparison of relative expression levels (percentage of control levels) between indicated LOF. Gene dots are coloured according to its expression along the animal-vegetal (AnVg) or dorso-ventral (DV) axis (with and without gene name labels).
:: venn_mPS_mVegT_minWNB_LOF.pdf, venn_WNB_mPS_LOF.pdf and venn_WNB_mVegT_LOF.pdf (Figure 7b and Supplementary Figure 10a,b)Venn diagrams of downregulated genes under indicated LOFs:: zga_lof_spatial_AnVg_perc.pdf and zga_lof_spatial_DV_perc.pdf (Figure 7e)Grouped bargraph showing the percentage (and # genes) of regional (across the animal-vegetal or dorso-ventral axis) ZGA downregulated by indicated LOFs:: zga_ama50_misregulation.pdf (Supplementary Figure 11a)
Bargraph showing the level of gene misregulation under indicated LOFs

:: GO_term_LOF_ZGA.pdf (Supplementary Figure 11b)GO term analysis of down- and upregulated genes during ZGA under indicated LOFs:: zga_lof.csv
Table of LOF effects on ZGA (restricted to alfa-amanitin downregulated (<= 50%, FDR <= 10%) genes (Supplementary Table 9)

:: spatial_down_ALL_LOF.txt, spatial_down_amanitin.txt, spatial_down_Bmp_LOF.txt, spatial_down_cMO.txt, spatial_down_mPouV_LOF.txt, spatial_down_mPS_LOF.txt, spatial_down_mTF_LOF.txt, spatial_down_mVegT_LOF.txt, spatial_down_Nodal_LOF.txt, spatial_down_Signals_LOF.txt and spatial_down_Wnt_LOF.txt
Text files the numbers/percentage of genes downregulated by indicated LOFs

:: AnVg_log2_colorBar.pdf, DV_log2_colorBar.pdf and FC_ratios_bar.pdfColour bars for spatial expression along the animal-vegetal axis and the dorso-ventral axis, and fold change ratios


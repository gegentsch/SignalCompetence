System requirements
-------------------

Running of script requires R. It can be downloaded from one of the many Comprehensive R Archive Network (CRAN) mirror websites of the R Project for Statistical Computing. 

https://cran.r-project.org/mirrors.html

The script was tested with R version 3.4.1 on a Darwin operating system supporting 64-bit x86-64 variant of the Intel x86 processors.


Running R script
----------------

This script can be run from the command line:
R CMD batch 6_extras.R

or with R studio or the R GUI:
source("6_extras.R")


Expected output
---------------

The following files will be produced upon run execution. Input files were generated as indicated here and in more details in the Online Methods.


# Supplementary Figure 1a: RNAPII+ pCRMs

RNAPII_pCRMs_dynamics.txt
Temporal dynamics of RNAPII-engaged (RNAPII+) pCRMs (≥1 ChIP tag/million) from the 32-cell to the late gastrula stage. 

The input file (./xenTro71/pol2_tagCount_ismara.txt.gz) was generated from RNAPII ChIP-Seq data using HOMER (tags/10 million reads).


# Supplementary Figure 1d: RNAPII/H3K4me1 levels at accessible pCRMs

dhs_pol2_h3k4me1_st8p.pdf
Biplot shows the DNA occupancy levels of RNAPII and H3K4me1 at accessible pCRMs. pCRMs (dots) are colored according to their chromatin accessibility except for the pCRMs that showed RNAPII and/or H3K4me1 signals piled up over a distance of 1 kb below background.

chromatin_access_colorbar.pdf
Color scale for chromatin accessibility (log scale)

The input file (./xenTro71/supplementary_files/dhs_pol2_h3k4me1_st8p_input_tags_1Kbp.csv) was generated using HOMER (tags/10 million reads).

# Figure 1f: Chromatin accessibility versus H3K4me1 and RNAPII

metaDHSst8ptop_DNase_logstd.pdf
metaDHSst8ptop_RNAPII_logstd.pdf
metaDHSst8ptop_H3K4me1_logstd.pdf

Metaplots of chromatin accessibility, RNAPII and H3K4me1

The input file (./xenTro71/supplementary_files/dhs_st8p_tagDensity1K25bp_DNase_pol2_h3k4me1.txt.gz) was generated using HOMER (tags/10 million reads).


# Figure 3c: Maternal TFs/signal mediators versus RNAPII

metabcat_st8_top2K_tagDensity2K25bp_bcat_pol2_logstd.pdfmetabcat_st10_top2K_tagDensity2K25bp_bcat_pol2_logstd.pdfmetafoxh1_st8_top2K_tagDensity2K25bp_foxh1_pol2_logstd.pdfmetafoxh1_st10_top2K_tagDensity2K25bp_foxh1_pol2_logstd.pdfmetasmad1_st8_top2K_tagDensity2K25bp_smad1_pol2_logstd.pdfmetasmad1_st10_top2K_tagDensity2K25bp_smad1_pol2_logstd.pdfmetasmad2_st8_top2K_tagDensity2K25bp_smad2_pol2_logstd.pdfmetasmad2_st10_top2K_tagDensity2K25bp_smad2_pol2_logstd.pdfmetasox3_st8_top2K_tagDensity2K25bp_sox3_pol2_logstd.pdfmetasox3_st10_top2K_tagDensity2K25bp_sox3_pol2_logstd.pdfmetavegT_st8_top2K_tagDensity2K25bp_vegT_pol2_logstd.pdfmetavegT_st10_top2K_tagDensity2K25bp_vegT_pol2_logstd.pdf
Meta-plots summarize the level of RNAPII engagement across 2,000 pCRMs most frequently occupied by the indicated TFs or signal mediators at the 1,024-cell and early gastrula stage, respectively.

The input files (./xenTro71/supplementary_files/_pol2/*.txt.gz) were generated using HOMER (tags/10 million reads).


# Figure 3f: TFs versus signal mediators (at TF+ pCRMs)

metaFoxH1st8top2K_foxH1_bcat_logstd.pdfmetaFoxH1st8top2K_foxH1_smad1_logstd.pdfmetaFoxH1st8top2K_foxH1_smad2_logstd.pdfmetaFoxH1st10top2K_foxH1_bcat_logstd.pdfmetaFoxH1st10top2K_foxH1_smad1_logstd.pdfmetaFoxH1st10top2K_foxH1_smad2_logstd.pdfmetaSox3st8top2K_sox3_bcat_logstd.pdfmetaSox3st8top2K_sox3_smad1_logstd.pdfmetaSox3st8top2K_sox3_smad2_logstd.pdfmetaSox3st10top2K_sox3_bcat_logstd.pdfmetaSox3st10top2K_sox3_smad1_logstd.pdfmetaSox3st10top2K_sox3_smad2_logstd.pdfmetaVegTst8top2K_vegT_bcat_logstd.pdfmetaVegTst8top2K_vegT_smad1_logstd.pdfmetaVegTst8top2K_vegT_smad2_logstd.pdfmetaVegTst10top2K_vegT_bcat_logstd.pdfmetaVegTst10top2K_vegT_smad1_logstd.pdfmetaVegTst10top2K_vegT_smad2_logstd.pdf
Meta-plots summarizes the level of signal mediator engagement across 2,000 pCRMs most frequently occupied by the indicated TF at the 1,024-cell and early gastrula stage, respectively.

The input files (./xenTro71/supplementary_files/_signal_mediators/*.txt.gz) were generated using HOMER (tags/10 million reads).


# Figure 8g: Effect of mPouV/Sox3 LOF on chromatin accessibility and RNAPII-mediated gene expression

AllZygoticGenes_DHS_mPSLOF.pdf
AllZygoticGenes_mPSLOF_log2.pdf
FDR_colorbar.pdf
Localization of accessible pCRMs (affected, dot colored in orange to red with FDR decreasing from 10%; and unaffected, grey dot) relative to the zygotic TSSs that are active by the MBT and produce enough RNA transcripts to show significant >=two-fold reductions upon alfa-amanitin injection). Gene loci are sorted by mPouV/Sox3 LOF-induced transcript fold changes as shown in the log-scaled bar graph.


# Supplementary Figure 15: mPouV/Sox3-induced chromatin accessibility is required for the expression of Nodal-responsive genes

NodalLOF_down_DHS_mPSLOF.pdf
NodalLOF_down_DHS_mPSLOF_RNA.png


# Supplementary Figure 16a: mPouV/Sox3-induced chromatin accessibility is required for the expression of Wnt-responsive genes

WntLOF_down_DHS_mPSLOF.pdf
WntLOF_down_DHS_mPSLOF_RNA.png


# Supplementary Figure 16b: mPouV/Sox3-induced chromatin accessibility is required for the expression of Bmp-responsive genes

BmpLOF_down_DHS_mPSLOF.pdf
BmpLOF_down_DHS_mPSLOF_RNA.png


# Supplementary Figure 17: mPouV/Sox3-induced chromatin accessibility is required for the expression of signal non-responsive genes

SignalLOF_NOTdown_DHS_mPSLOF.pdf
SignalLOF_NOTdown_mPSLOF.png
SignalLOF_NOTdown_DHS_mPSLOF_highRes.pdf (high resolution of SignalLOF_NOTdown_DHS_mPSLOF.pdf, gene names are readable)


All genes listed in Supplementary Figure 15-17 are active by the MBT and their transcript levels are significantly reduced (≥two-fold; FDR ≤10%) upon alfa-amanitin injection. The plot aligned shows the localization and DNase sensitivity (bubble size) of accessible pCRMs (affected, dot colored in orange to red with FDR decreasing from 10%; and unaffected, grey dot) relative to zygotic TSSs. Gene loci are sorted by mPouV/Sox3 LOF-induced transcript fold changes as shown in the heat map.

The input files for Figure 8g and Supplementary Figure 15-17 (./xenTro71/supplementary_files/... dhs_fdr_genes_50kb_distanceTSS.bed.gz, zga_lof.csv, pol2_st8to8p_rzRNA_av3tpm.bed, dhs_motif_matrix.txt.gz) were generated using previous R scripts and HOMER (see Online Methods).


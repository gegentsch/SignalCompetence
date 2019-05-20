System requirements
-------------------

Running of script requires R. It can be downloaded from one of the many Comprehensive R Archive Network (CRAN) mirror websites of the R Project for Statistical Computing. 

https://cran.r-project.org/mirrors.html

The script was tested with R version 3.4.1 on a Darwin operating system supporting 64-bit x86-64 variant of the Intel x86 processors.


Running R script
----------------

Before running this script, align ChIP-Seq reads (FASTQ) to the genome assembly 7.1 of X. tropicalis. 
FASTQ files can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113186

To make experiment 2 (paired-end reads) comparable to experiment 1 (single reads), we did not use read 2 of experiment 2 in this analysis.

Add path/BAM files to the 'DiffBind' input (.csv) files provided under ./xenTro71/supplementary_files
dnase.csv
dnase_mPS_LOF_affected_capC.csvdnase_mPS_LOF_NOTaffected_capC.csv


This script can be run from the command line:
R CMD batch 5_mPS_LOF.R

or with R studio or the R GUI:
source("5_mPS_LOF.R")


Expected output
---------------

This script quantifies the changes to the chromatin landscape that were detectable at the mid-blastula transition (MBT) in the absence of maternal Pou5f3 and Sox3 (mPS LOF). 

The following files will be produced upon run execution.

:: dhs_mPS_LOF_heatmap.png (Figure 8d,e)
Heatmap organised in columns (#) showing the following from left to right. The CRMs are sorted by the FDR reflecting the significant loss of chromatin accessibility caused by mPS LOF.
(1) 	DNA occupancy of Sox3 (ChIP-Seq)
(2-3)	Chromatin accessibility (DNase) for control (uninjected) and mPS LOF embryos
(4)	p-value of differential accessibility between mPS LOF and control embryos
(5-6)	H3K4me1 deposition (ChIP-Seq)	for control (uninjected) and mPS LOF embryos
(7)	p-value of differential H3K4me1 between mPS LOF and control embryos
(8-9)	DNA occupancy of Smad2 (ChIP-Seq)
(10)	p-value of differential Smad2 between mPS LOF and control embryos
(11-12)	DNA occupancy of beta-catenin (ChIP-Seq) for control (uninjected) and mPS LOF embryos
(13)	p-value of differential beta-catenin between mPS LOF and control embryos
(14-15)	DNA occupancy of RNAPII (ChIP-Seq) for control (uninjected) and mPS LOF embryos
(16)	p-value of differential RNAPII between mPS LOF and control embryos
(17-18)	(Non-coding) RNA transcript levels for control (uninjected) and mPS LOF embryos
(19-26)	DNA motif occurence for POU.pou5f3, SOX.sox3, POU-SOX.pou5f3-sox3, FOXH.foxH1, T.box, ZF.3xC2H2.Sp, bZIP.max, CCAAT.NFY

:: metaSox3_ctrl_logstd.pdf, metaDNase_ctrl_logstd.pdf, metaDNase_mPSkd_logstd.pdf, metaH3K4me1_ctrl_logstd.pdf, metaH3K4me1_mPSkd_logstd.pdf, metaPol2_ctrl_logstd.pdf, metaPol2_mPSkd_logstd.pdf, metaSmad2_ctrl_logstd.pdf, metaSmad2_mPSkd_logstd.pdf, metaBcat_ctrl_logstd.pdf, metaBcat_mPSkd_logstd.pdf, metaRNA_ctrl_logstd.pdf and metaRNA_mPSkd_logstd.pdf (Figure 8d)
Meta-profile (density) of various chromatin features across accessible CRMs (mean +/- standard deviation on log scale)

:: metaPOU.pdf, metaSOX.pdf, metaPOU-SOX.pdf, metaFOX.pdf, metaT.pdf, metaZnF.pdf, metaBZIP.pdf and metaNFY.pdf (Figure 8e)
Meta-profile (density) of various DNA motifs across accessible CRMs (mean +/- standard error)

:: mPS_LOF_tss_dist_pie_chart.pdf (Figure 8f)
Pie chart of TSS proximity of affected and unaffected CRMs

:: dnase_mPS_LOF_spearman.pdf (Supplementary Figure 12a)
Heatmap of Spearman correlations between indicated DNase-Seq samples

:: DBA_dnase_mPS_LOF_counts.csv (Supplementary Table 11)
Spreadsheet containing normalised DNase cleavage numbers for all detected DNase-hypersensitive cis-regulatory modules (CRMs)

:: dnase_mPS_LOF_volcano.pdf (Figure 8b)
Volcano plot of DNase cleavage fold changes between mPS LOF and control (uninjected) embryos at MBT

:: dnase_mPS_LOF_vioplot.pdf (Figure 8c)
Violin plot showing the number of DNase cleavages at all and affected (FDR <= 10%) CRMs in mPS LOF and control embryos at MBT

:: DBA_dhs_mPS_LOF_capC_affected.csv and DBA_dhs_mPS_LOF_capC_NOTaffected.csv (Supplementary Table 12)
Spreadsheet containing normalised chromatin contact numbers at affected (FDR <= 10%) and unaffected CRMs

:: dhs_capC_mPS_LOF_vioplot.pdf (Figure 9c)
Violin plot showing the number of chromatin contacts at all and affected (FDR <= 10%, reduced accessibility) CRMs 

:: DHS_capC_statistics.txt (Figure 9c)
Statistical significance of chromatin contact losses at CRMs in mPS LOF embryos

:: ColourBar_mPS_LOF_DNAocc_scaled.pdf, ColourBar_mPS_LOF_TSSproximity.pdf and ColourBar_mPS_LOF_diffPval.pdf
Colour bars for scaled DNA occupancy, CRM-TSS proximity and differential p-value
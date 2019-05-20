# CHROMATIN LANDSCAPE CHANGES INDUCED BY MATERNAL POU5F3 and SOX3
# Author: G.E. Gentsch
#
#
#
# ChIP-Seq and DNase-Seq reads were aligned and peaks detected as described in the Methods section of the paper.
#
# The input templates (.csv) for the DiffBind R package and the genome coordinates in BED format are provided under directory ./xenTro71/supplementary_files. The BAM files need to be generated de novo from the raw FASTQ files deposited to GEO archive: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113186


# Installing (if required) and uploading several R packages
required.pkg <- c( "seriation", "DESeq2", "extrafont", "gplots", "DiffBind", "vioplot", "scales" )

for ( pkg in required.pkg ) {
    if ( pkg %in% rownames( installed.packages() ) == FALSE) {
        install.packages( pkg ) }
    if ( pkg %in% rownames(.packages()) == FALSE ) {
        library( pkg, character.only = TRUE ) }
}

source( "./R_utils/colour_schemes.R" )
source( "./R_utils/plot_corr.R" )
source( "./R_utils/plot_motif.R" )
source( "./R_utils/p_value_matrix.R" )
source( "./R_utils/meta_density.R" )




# To make experiment 2 (paired-end reads) comparable to experiment 1 (single reads), we only used read 1 of experiment 2.
d           <- dba( sampleSheet = "./xenTro71/supplementary_files/dnase.csv" )

d           <- dba.count( d, minOverlap=0, score=DBA_SCORE_READS, bRemoveDuplicates=TRUE, mapQCth=10, fragmentSize=1, summits=250, bParallel=TRUE )

# set contrast (none versus mPSkd; DBA_TREATMENT)
d           <- dba.contrast( d, categories=DBA_TREATMENT, minMembers=2 )

d           <- dba.analyze( d, method=DBA_DESEQ2, bFullLibrarySize=TRUE, bTagwise=FALSE, bSubControl=TRUE )

d.spearman  <- dba.plotHeatmap(d, distMethod="spearman")

brk         <- c( seq( min(d.spearman), 1, length.out=80 ) )
pdf( "dnase_mPS_LOF_spearman.pdf", width=6, height=6)
heatmap.2(
 d.spearman, Rowv=FALSE, Colv=FALSE, dendrogram=c("none"), trace="none", symm=TRUE, scale="none", breaks=brk, col=col.corr2, key=TRUE, symkey=FALSE, density.info="none", cellnote=round(d.spearman, 2), notecol="white", notecex=1.2, cexRow=1.2, cexCol=1.2, margin=c(13,13) )
dev.off()

d.counts    <- dba.report( d, th=1, bCounts=TRUE, file="dnase_mPS_LOF_counts" )

dhs         <- read.csv( "DBA_dnase_mPS_LOF_counts.csv", header=TRUE )
dhs$dnase_st8p_mPS_LOF  <- log2( ( dhs$dnase_st8p_mPS_LOF_1 + dhs$dnase_st8p_mPS_LOF_2 ) / 2 )
dhs$dnase_st8p_ctrl     <- log2( ( dhs$dnase_st8p_ctrl_1 + dhs$dnase_st8p_ctrl_2 ) / 2 )

pdf( "dnase_mPS_LOF_volcano.pdf", width=3, height=3, useDingbats=FALSE, family="Arial")
    # Make a basic volcano plot
    with(subset(dhs, FDR>=.1), plot(Fold, -log10(FDR), pch=20, main="Volcano plot", xlim=c(-7.5,7.5), ylim=c(0,15), col=alpha("black",.4),cex=.6,lwd=0))
    # Add colored points: red, if FDR<10%
    with(subset(dhs, FDR<.1), points(Fold, -log10(FDR), pch=20, col=alpha("#EE3B38",.4),cex=.6,lwd=0))
dev.off()

v1          <- dhs$dnase_st8p_ctrl
v2          <- dhs$dnase_st8p_mPS_LOF
v3          <- dhs$dnase_st8p_ctrl[ dhs$FDR<=0.1 ]
v4          <- dhs$dnase_st8p_mPS_LOF[ dhs$FDR<=0.1 ]

pdf( "dnase_mPS_LOF_vioplot.pdf", width=2.5, height=3.5, useDingbats=FALSE, family="Arial" )
    vioplot(v1, v2, v3, v4, names=c( "ALL.ctrl", "ALL.mPS.LOF", "affected.ctrl", "affected.mPS.LOF" ), col="#C9992B" )
    title("Chromatin accessibility")
dev.off()

wt1         <- wilcox.test(v2, v1, alternative="less")
wt2         <- wilcox.test(v4, v3, alternative="less")

win.dnase   <- read.table( "./xenTro71/supplementary_files/dnase_mPS_LOF_matrix.txt.gz", header=TRUE )
win.h3k4me1 <- read.table( "./xenTro71/supplementary_files/h3k4me1_mPS_LOF_matrix.txt.gz", header=TRUE )
win.pol2    <- read.table( "./xenTro71/supplementary_files/pol2_mPS_LOF_matrix.txt.gz", header=TRUE )
win.smad2   <- read.table( "./xenTro71/supplementary_files/smad2_mPS_LOF_matrix.txt.gz", header=TRUE )
win.bcat    <- read.table( "./xenTro71/supplementary_files/bcat_mPS_LOF_matrix.txt.gz", header=TRUE )

dnase.pval  <- pValueMatrix(win.dnase)
h3k4me1.pval <- pValueMatrix(win.h3k4me1)
pol2.pval   <- pValueMatrix(win.pol2)
smad2.pval  <- pValueMatrix(win.smad2)
bcat.pval   <- pValueMatrix(win.bcat)

mPS          <- read.table( "./xenTro71/supplementary_files/dhs_chromatin_RNA_mPS_LOF_matrix.txt.gz", header=T, row.names="Gene" )

# Calculate average of top 1000 peaks (choose 25bp-bin with highest read count) to scale profiles
n.profiles   <- as.numeric( unlist( strsplit( colnames( mPS )[ ncol(mPS) ], ".", fixed=TRUE ))[2] ) + 1
n.bins       <- ncol(mPS)/n.profiles
mPS.scale    <- mPS

# Scaling
for (i in seq(1,ncol(mPS),n.bins)) {
    mPS.scale[,i:(i+n.bins-1)] <- 100 / mean(sort(apply(mPS[,i:(i+n.bins-1)],1,max),decreasing=T)[1:1000]) * mPS[,i:(i+n.bins-1)]}

# Change scaling for RNAPII (mean of top 1000 divided by 5)
for (i in seq((5*n.bins+1),(5*n.bins+2*n.bins),n.bins)) {
    mPS.scale[,i:(i+n.bins-1)] <- 100 / mean(sort(apply(mPS[,i:(i+n.bins-1)],1,max),decreasing=T)[1:1000]/5) * mPS[,i:(i+n.bins-1)]}

# Change scaling for RNA (mean of top 1000 divided by 500)
for (i in seq((11*n.bins+1),(11*n.bins+2*n.bins),n.bins)) {
    mPS.scale[,i:(i+n.bins-1)] <- 100 / mean(sort(apply(mPS[,i:(i+n.bins-1)],1,max),decreasing=T)[1:1000]/500) * mPS[,i:(i+n.bins-1)]}

# Motifs found within +/-1kb
mo            <- read.table( "./xenTro71/supplementary_files/dhs_motif_matrix.txt.gz", header=TRUE, row.names="Gene" )

# Distance from TSS
tss           <- read.table( "./xenTro71/supplementary_files/dhs_closest_ActiveGenes.txt.gz", header=FALSE )
colnames(tss) <- c("scaffold","start","end","dhs","fdr","strand","scaffold.tss","start.tss","end.tss","gene","empty","strand.tss","dist")


# Define breaks and colour for heatmap
bre1 <- c(seq(0,50,length=100),seq(50.5,200,length=60),seq(100000,max(mPS.scale),length=40))
col1 <- colorRampPalette(c("#EBEBEB","black","#FF0000","#FFFFCC"))(199)
bre2 <- c(1,max(mo))
col2 <- "black"
bre3 <- c(min(tss$dist),-10000,-501,500,10000,max(tss$dist))
col3 <- c("lightblue","#4575B4","black","yellow","lightyellow")
bre4 <- c(seq(0,8,length=90),seq(9,max(pol2.pval),length=10))
col4 <- colorRampPalette(c("#FFFFFF","#FF0000"))(99)


# Colour bars

# scaled DNA occupancy
z <- matrix(1:200,nrow=1)
pdf("ColourBar_mPS_LOF_DNAocc_scaled.pdf", width=2, height=4)
  par(mar=c(0,0,0,0) + 4)
  image(1,bre1,z,col=col1,axes=F, ylim=c(0,200), xlab="", ylab="DNA occupancy (scaled)")
  axis(2, at=c(0,50,100,150,200), labels=c("0","0.5","1","1.5","2"))
dev.off()

# TSS proximity
z <- matrix(1:60,nrow=1)
pdf("ColourBar_mPS_LOF_TSSproximity.pdf", width=2, height=4)
  par(mar=c(0,0,0,0) + 4)
  image(1,seq(-14750,14750,length=60),z,col=c(rep("lightblue",10),rep("#4575B4",19),rep("black",2),rep("yellow",19),rep("lightyellow",10)),axes=F,ylim=c(-15000,15000),xlab="",ylab="TSS proximity (kb)")
  axis(2, at=c(-15000,-10000,-5000,0,5000,10000,15000), labels=c("-15","-10","-5","TSS","5","10","15"))
dev.off()

# p-values
z <- matrix(1:100,nrow=1)
pdf("ColourBar_mPS_LOF_diffPval.pdf", width=2, height=4)
  par(mar=c(0,0,0,0) + 4)
  image(1,bre4,z,col=col4,axes=F, ylim=c(0,10), xlab="", ylab="-log(p)")
  axis(2, at=c(0,5,10), labels=c("0","5","10"))
dev.off()


png(file="dhs_mPS_LOF_heatmap.png", width=n.bins*28, height=nrow(mPS))

par(mar=c(0,0,0,0) + 0.1)

layout(matrix(c(rep(1,1),rep(2,2),rep(3,1),rep(4,2),rep(5,1),rep(6,2),rep(7,1),rep(8,2),rep(9,1),rep(10,2),rep(11,1),rep(12,2),rep(13,8)), 1, 26, byrow = TRUE))

image( t(mPS.scale[,1:n.bins]), breaks=bre1, col=col1, axes=FALSE ) #1: DNA occupancy of Sox3

image( t(mPS.scale[,(n.bins+1):(3*n.bins)]), breaks=bre1, col=col1, axes=FALSE ) #2: Chromatin accessibility
image( t(dnase.pval), breaks=bre4, col=col4, axes=FALSE ) #3: Significance of differential accessibility (p-value)

image( t(mPS.scale[,(3*n.bins+1):(5*n.bins)]), breaks=bre1, col=col1, axes=FALSE ) #4: H3K4me1 deposition
image( t(h3k4me1.pval), breaks=bre4, col=col4, axes=FALSE ) #5: Significance of differential H3K4me1 (p-value)

image( t(mPS.scale[,(7*n.bins+1):(9*n.bins)]), breaks=bre1, col=col1, axes=FALSE ) #6: DNA occupancy of Smad2
image( t(smad2.pval), breaks=bre4, col=col4, axes=FALSE ) #7: Significance of differential Smad2 (p-value)

image( t(mPS.scale[,(9*n.bins+1):(11*n.bins)]), breaks=bre1, col=col1, axes=FALSE ) #8: DNA occupancy of b-catenin
image( t(bcat.pval), breaks=bre4, col=col4, axes=FALSE ) #9: Significance of differential b-catenin (p-value)

image( t(mPS.scale[,(5*n.bins+1):(7*n.bins)]), breaks=bre1, col=col1, axes=FALSE ) #10: DNA occupancy of RNAPII
image( t(pol2.pval), breaks=bre4, col=col4, axes=FALSE ) #11: Significance of differential RNAPII (p-value)

image( t(mPS.scale[,(11*n.bins+1):(13*n.bins)]), breaks=bre1, col=col1, axes=FALSE ) #12: (non-coding) RNA levels

image( t(mo), breaks=bre2, col=col2, axes=FALSE )  #13: DNA motifs in the following order: POU.pou5f3, SOX.sox3, POU-SOX.pou5f3-sox3, FOXH.foxH1, T.box, ZF.3xC2H2.Sp, bZIP.max, CCAAT.NFY
dev.off()

# import into Photoshop and reduce height to 4% (bicubic Smoother, ideal for enlargments)


# Plot read densities (mean + standard deviation, logarithmic scale) across affected / non-affected accessible chromatin regions
# ---------------------------------------------------------------------------------------------------------------------------------------

thres   <- length(which(-log10(d.counts$FDR)>1))
nr      <- nrow(mPS.scale)

pdf("metaSox3_ctrl_logstd.pdf", width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[1:thres,3:(n.bins-2)],"#FF0000", 0.1, 75)
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[(thres+1):nr,3:(n.bins-2)],"#4C62A4", 0.1, 75)
dev.off()

pdf("metaDNase_ctrl_logstd.pdf", width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[1:thres,(n.bins+3):(2*n.bins-2)],"#FF0000", 0.1, 75)
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[(thres+1):nr,(n.bins+3):(2*n.bins-2)],"#4C62A4", 0.1, 75)
dev.off()

pdf("metaDNase_mPSkd_logstd.pdf", width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[1:thres,(2*n.bins+3):(3*n.bins-2)],"#FF0000", 0.1, 75)
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[(thres+1):nr,(2*n.bins+3):(3*n.bins-2)],"#4C62A4", 0.1, 75)
dev.off()

pdf("metaH3K4me1_ctrl_logstd.pdf", width=2, height=1.66666)
  par(mfg=c(1,1),mar=c(2,2,2,2))
  plotMetaLogNorm(mPS.scale[1:thres,(3*n.bins+3):(4*n.bins-2)],"#FF0000", 0.1, 75)
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[(thres+1):nr,(3*n.bins+3):(4*n.bins-2)],"#4C62A4", 0.1, 75)
dev.off()

pdf("metaH3K4me1_mPSkd_logstd.pdf", width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[1:thres,(4*n.bins+3):(5*n.bins-2)],"#FF0000", 0.1, 75)
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[(thres+1):nr,(4*n.bins+3):(5*n.bins-2)],"#4C62A4", 0.1, 75)
dev.off()

pdf("metaPol2_ctrl_logstd.pdf", width=2, height=1.66666)
  par(mfg=c(1,1),mar=c(2,2,2,2))
  plotMetaLogNorm(mPS.scale[1:thres,(5*n.bins+3):(6*n.bins-2)],"#FF0000", 0.1, 75)
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[(thres+1):nrow(mPS.scale),(5*n.bins+3):(6*n.bins-2)],"#4C62A4", 0.1, 75)
dev.off()

pdf("metaPol2_mPSkd_logstd.pdf", width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[1:thres,(6*n.bins+3):(7*n.bins-2)],"#FF0000", 0.1, 75)
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[(thres+1):nrow(mPS.scale),(6*n.bins+3):(7*n.bins-2)],"#4C62A4", 0.1, 75)
dev.off()

pdf("metaSmad2_ctrl_logstd.pdf", width=2, height=1.66666)
  par(mfg=c(1,1),mar=c(2,2,2,2))
  plotMetaLogNorm(mPS.scale[1:thres,(7*n.bins+3):(8*n.bins-2)],"#FF0000", 0.1, 75)
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[(thres+1):nrow(mPS.scale),(7*n.bins+3):(8*n.bins-2)],"#4C62A4", 0.1, 75)
dev.off()

pdf("metaSmad2_mPSkd_logstd.pdf", width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[1:thres,(8*n.bins+3):(9*n.bins-2)],"#FF0000", 0.1, 75)
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[(thres+1):nrow(mPS.scale),(8*n.bins+3):(9*n.bins-2)],"#4C62A4", 0.1, 75)
dev.off()

pdf("metaBcat_ctrl_logstd.pdf", width=2, height=1.66666)
  par(mfg=c(1,1),mar=c(2,2,2,2))
  plotMetaLogNorm(mPS.scale[1:thres,(9*n.bins+3):(10*n.bins-2)],"#FF0000", 0.1, 75)
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[(thres+1):nrow(mPS.scale),(9*n.bins+3):(10*n.bins-2)],"#4C62A4", 0.1, 75)
dev.off()

pdf("metaBcat_mPSkd_logstd.pdf", width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[1:thres,(10*n.bins+3):(11*n.bins-2)],"#FF0000", 0.1, 75)
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[(thres+1):nrow(mPS.scale),(10*n.bins+3):(11*n.bins-2)],"#4C62A4", 0.1, 75)
dev.off()

pdf("metaRNA_ctrl_logstd.pdf", width=2, height=1.66666)
  par(mfg=c(1,1),mar=c(2,2,2,2))
  plotMetaLogNorm(mPS.scale[1:thres,(11*n.bins+3):(12*n.bins-2)],"#FF0000", 0.01, 2000)
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[(thres+1):nrow(mPS.scale),(11*n.bins+3):(12*n.bins-2)],"#4C62A4", 0.01, 2000)
dev.off()

pdf("metaRNA_mPSkd_logstd.pdf", width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[1:thres,(12*n.bins+3):(13*n.bins-2)],"#FF0000", 0.01, 2000)
  par(mfg=c(1,1))
  plotMetaLogNorm(mPS.scale[(thres+1):nrow(mPS.scale),(12*n.bins+3):(13*n.bins-2)],"#4C62A4", 0.01, 2000)
dev.off()



# Plot DNA motif densities (mean + standard deviation) across affected / non-affected accessible chromatin regions
# ------------------------------------------------------------------------------------------------------------------------------

pdf("metaPOU.pdf", width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaMotif(mo[1:thres,3:(n.bins-2)],"#FF0000", 0.2)
  par(mfg=c(1,1))
  plotMetaMotif(mo[(thres+1):nr,3:(n.bins-2)],"#4C62A4", 0.2)
dev.off()

pdf("metaSOX.pdf", width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaMotif(mo[1:thres,(n.bins+3):(2*n.bins-2)],"#FF0000", 0.2)
  par(mfg=c(1,1))
  plotMetaMotif(mo[(thres+1):nr,(n.bins+3):(2*n.bins-2)],"#4C62A4", 0.2)
dev.off()

pdf("metaPOU-SOX.pdf", width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaMotif(mo[1:thres,(2*n.bins+3):(3*n.bins-2)],"#FF0000", 0.2)
  par(mfg=c(1,1))
  plotMetaMotif(mo[(thres+1):nr,(2*n.bins+3):(3*n.bins-2)],"#4C62A4", 0.2)
dev.off()

pdf("metaFOX.pdf", width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaMotif(mo[1:thres,(3*n.bins+3):(4*n.bins-2)],"#FF0000", 0.2)
  par(mfg=c(1,1))
  plotMetaMotif(mo[(thres+1):nr,(3*n.bins+3):(4*n.bins-2)],"#4C62A4", 0.2)
dev.off()

pdf("metaT.pdf", width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaMotif(mo[1:thres,(4*n.bins+3):(5*n.bins-2)],"#FF0000", 0.2)
  par(mfg=c(1,1))
  plotMetaMotif(mo[(thres+1):nr,(4*n.bins+3):(5*n.bins-2)],"#4C62A4", 0.2)
dev.off()

pdf("metaZnF.pdf", width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaMotif(mo[1:thres,(5*n.bins+3):(6*n.bins-2)],"#FF0000", 0.2)
  par(mfg=c(1,1))
  plotMetaMotif(mo[(thres+1):nr,(5*n.bins+3):(6*n.bins-2)],"#4C62A4", 0.2)
dev.off()

pdf("metaBZIP.pdf", width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaMotif(mo[1:thres,(6*n.bins+3):(7*n.bins-2)],"#FF0000", 0.2)
  par(mfg=c(1,1))
  plotMetaMotif(mo[(thres+1):nr,(6*n.bins+3):(7*n.bins-2)],"#4C62A4", 0.2)
dev.off()

pdf("metaNFY.pdf", width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaMotif(mo[1:thres,(7*n.bins+3):(8*n.bins-2)],"#FF0000", 0.2)
  par(mfg=c(1,1))
  plotMetaMotif(mo[(thres+1):nr,(7*n.bins+3):(8*n.bins-2)],"#4C62A4", 0.2)
dev.off()


pie.col   <- c( "lightblue", "#4575B4", "black", "yellow", "lightyellow" )
brk.dist  <- c(min(tss$dist),-10000,-501,500,10000,max(tss$dist))
af        <- hist( tss$dist[1:thres],breaks=brk.dist,xlim=c(-12000,12000))
uf        <- hist( htss$dist[(nr-thres+1):nr],breaks=brk.dist,xlim=c(-12000,12000))

pdf("mPS_LOF_tss_dist_pie_chart.pdf", height=1, width=7)
  par(mfrow=c(1,2), mar=c(0,0,1,0) + 0.1)
  pie(af$counts, clockwise=T, labels = NA, col=pie.col, main="affected accessibility")
  pie(uf$counts, clockwise=T, labels = NA, col=pie.col, main="unaffected accessibility")
dev.off()



# Effect of chromatin pioneering (loss of chromatin accessibility) on chromatin conformation (NG capture-C)
# -----------------------------------------------------------------------------------------------------------

d.capC        <- dba( sample="./xenTro71/supplementary_files/dnase_mPS_LOF_affected_capC.csv" )
d.capC.ctrl   <- dba( sample="./xenTro71/supplementary_files/dnase_mPS_LOF_NOTaffected_capC.csv" )

d.capC        <- dba.count( d.capC, minOverlap=0, score=DBA_SCORE_READS, bRemoveDuplicates=TRUE, mapQCth=10, fragmentSize=0, summits=600, filter=5, bParallel=TRUE )
d.capC.ctrl   <- dba.count( d.capC.ctrl, minOverlap=0, score=DBA_SCORE_READS, bRemoveDuplicates=TRUE, mapQCth=10, fragmentSize=0, summits=600, filter=5, bParallel=TRUE )

# set contrast (none versus mPouV/Sox3 LOF; DBA_TREATMENT)
d.capC        <- dba.contrast( d.capC, categories=DBA_TREATMENT, minMembers=3 )
d.capC.ctrl   <- dba.contrast( d.capC.ctrl, categories=DBA_TREATMENT, minMembers=3 )

d.capC        <- dba.analyze( d.capC, method=DBA_DESEQ2, bFullLibrarySize=TRUE, bTagwise=FALSE, bSubControl=FALSE )
d.capC.ctrl   <- dba.analyze( d.capC.ctrl, method=DBA_DESEQ2, bFullLibrarySize=TRUE, bTagwise=FALSE, bSubControl=FALSE )

# report saved as csv file (save in Excel as Windows Formatted Text File)
d.capC.counts       <- dba.report( d.capC, th=1, bCounts=TRUE, file="dhs_mPS_LOF_capC_affected" )
d.capC.ctrl.counts  <- dba.report( d.capC.ctrl, th=1, bCounts=TRUE, file="dhs_mPS_LOF_capC_NOTaffected" )

dhs.capC            <- read.csv( "./xenTro71/supplementary_files/DBA_dhs_mPS_LOF_capC_affected.csv", header=TRUE )
dhs.capC.ctrl       <- read.csv( "./xenTro71/supplementary_files/DBA_dhs_mPS_LOF_capC_NOTaffected.csv", header=TRUE )
dhs.capC$capC_mPS_LOF <- log2( ( dhs.capC$capC_mPS_LOF_st8p_1 + dhs.capC$capC_mPS_LOF_st8p_2 + dhs.capC$capC_mPS_LOF_st8p_3 ) / 3 )
dhs.capC$capC_ctrl    <- log2( ( dhs.capC$capC_ctrl_st8p_1 + dhs.capC$capC_ctrl_st8p_2 + dhs.capC$capC_ctrl_st8p_3 ) / 3 )
dhs.capC.ctrl$capC_mPS_LOF   <- log2( ( dhs.capC.ctrl$capC_mPS_LOF_st8p_1 + dhs.capC.ctrl$capC_mPS_LOF_st8p_2 + dhs.capC.ctrl$capC_mPS_LOF_st8p_3 ) / 3 )
dhs.capC.ctrl$capC_ctrl      <- log2( ( dhs.capC.ctrl$capC_ctrl_st8p_1 + dhs.capC.ctrl$capC_ctrl_st8p_2 + dhs.capC.ctrl$capC_ctrl_st8p_3 ) / 3 )

l0.3    <- length(dhs.capC.ctrl$capC_ctrl[ dhs.capC.ctrl$FDR<=0.1 & dhs.capC.ctrl$Fold<0 ])
l0.4    <- length(dhs.capC.ctrl$capC_mPS_LOF[ dhs.capC.ctrl$FDR<=0.1 & dhs.capC.ctrl$Fold<0 ])

v0.1    <- dhs.capC.ctrl$capC_ctrl
v0.2    <- dhs.capC.ctrl$capC_mPS_LOF
v0.3    <- ifelse(l0.3 > 0, dhs.capC.ctrl$capC_ctrl[ dhs.capC.ctrl$FDR<=0.1 & dhs.capC.ctrl$Fold<0 ], 0)
v0.4    <- ifelse(l0.4 > 0, dhs.capC.ctrl$capC_mPS_LOF[ dhs.capC.ctrl$FDR<=0.1 & dhs.capC.ctrl$Fold<0 ], 0)
v1      <- dhs.capC$capC_ctrl
v2      <- dhs.capC$capC_mPS_LOF
v3      <- dhs.capC$capC_ctrl[ dhs.capC$FDR<=0.1 & dhs.capC$Fold<0 ]
v4      <- dhs.capC$capC_mPS_LOF[ dhs.capC$FDR<=0.1 & dhs.capC$Fold<0 ]

pdf("dhs_capC_mPS_LOF_vioplot.pdf", width=2.5, height=3, useDingbats=FALSE, family="Arial")
  vioplot( v0.1, v0.2, v0.3, v0.4, v1, v2, v3, v4, 
    names=c("UnaffCRM.ALLContact.ctrl", "UnaffCRM.ALLContacts.mPSLOF", "UnaffCRM.AffContacts.ctrl", "UnaffCRM.AffContacts.mPSkd",
            "AffCRM.ALLContacts.ctrl", "AffCRM.ALLContacts.mPSLOF", "AffCRM.AffContacts.ctrl", "AffCRM.AffContacts.mPSLOF"),
    col="#C9992B")
  title("chromatin conformation")
dev.off()

wt0     <- wilcox.test( v0.2, v0.1, alternative="two.sided", paired=TRUE )
r0      <- abs( qnorm( wt0$p.value / 2 ) / sqrt( length( c( v0.2, v0.1 ) ) ) )
wt1     <- wilcox.test( v2, v1, alternative="two.sided", paired=TRUE )
r1      <- abs( qnorm( wt1$p.value / 2 ) ) / sqrt( length( c( v1, v2 ) ) )
wt2     <- wilcox.test( v4, v3, alternative="two.sided", paired=TRUE )
r2      <- abs( qnorm( wt2$p.value / 2 ) ) / sqrt( length( c( v4, v3 ) ) )

# generate file with statitics
sink(file="DHS_capC_statistics.txt")

cat("# EFFECT OF CHROMATIN ACCESSIBILITY ON CHROMATIN CONFORMATION\n\n\n")

cat("# How significant is the effect of unaffected chromatin accessibility on all detected chromatin contacts?\n")
cat("Two-sided Wilcoxon Test (p-value)","Effect size estimate (r)\n", sep="\t")
cat(wt0$p.value, r0, "\n\n", sep="\t")

cat("# How significant is the effect of reduced chromatin accessibility (caused by mPouV/Sox3 LOF) on all detected chromatin contacts?\n")
cat("Two-sided Wilcoxon Test (p-value)","Effect size estimate (r)\n", sep="\t")
cat(wt1$p.value, r1, "\n\n", sep="\t")

cat("# How significant is the effect of reduced chromatin accessibility (caused by mPouV/Sox3 LOF) on affected chromatin contacts?\n")
cat("Two-sided Wilcoxon Test (p-value)","Effect size estimate (r)\n", sep="\t")
cat(wt2$p.value, r2, sep="\t")

sink()

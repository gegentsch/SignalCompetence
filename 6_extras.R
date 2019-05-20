# Title: The role of maternal pioneer factors in predefining first zygotic responses to inductive signals
# Author: G. Gentsch

# Extras generated for revised manuscript

# Installing (if required) and uploading several R packages
required.pkg <- c( "seriation", "extrafont", "ggplot2", "limma","doBy","lattice","latticeExtra","RColorBrewer","plyr","scales")
for ( pkg in required.pkg ) {
    if ( pkg %in% rownames(installed.packages()) == FALSE) {
        install.packages( pkg ) }
    if ( pkg %in% rownames(.packages()) == FALSE ) {
        library(pkg, character.only = TRUE) }
}


# ---- Supplementary Figure 1a: RNAPII+ pCRMs ----

# Temporal dynamics of RNAPII-engaged (RNAPII+) pCRMs (≥1 ChIP tag/million) from the 32-cell to the late gastrula stage

# Normalised RNAPII tag counts
pol2 <- read.table("./xenTro71/pol2_tagCount_ismara.txt.gz")

# Threshold ≥10 tags (per 10 million) per pCRM
pol2.total <- function(x) {
  length(which(x>=10))
  }

# Total RNAPII+ pCRMs per developmental stage
pol2.stage <- apply(pol2,2,pol2.total)

# Dynamics of RNAPII+ pCRMs 
pol2.bin <- pol2
pol2.bin[pol2.bin < 10] <- 0
pol2.bin[pol2.bin >= 10] <- 1

# Generate file with RNAPII statistics
sink( file = "RNAPII_pCRMs_dynamics.txt" )
cat( "## STATISTICS: RNAPII+ pCRMs (Template for Supplementary Fig. 1a) ##\n\n", sep="" )
cat( "## Counts of RNAPII+ pCRMs per developmental stage #\n" )
print( pol2.stage )
cat( "\n" )
cat( "## Counts of RNAPII+ pCRMs with specific engagement pattern over time (0, not engaged; 1, engaged) #\n" )
print( vennCounts( pol2.bin ))
sink()





# ---- Supplementary Figure 1d: RNAPII/H3K4me1 levels at accessible pCRMs ----

tags <- read.csv("./xenTro71/supplementary_files/dhs_pol2_h3k4me1_st8p_input_tags_1Kbp.csv",header=T)

col.h3k4me1 <- colorRampPalette(c("#EBEBEB","#CFA6A9","#AD6367"))(179)
brk.h3k4me1 <- seq(min(log10(tags$H3K4me1)),max(log10(tags$H3K4me1)),length=180)
cl.h3k4me1 <- ifelse( tags$RNAPII-tags$Input>0, col.h3k4me1 [ as.numeric( cut( log10(tags$H3K4me1), breaks=brk.h3k4me1 )) ], "black" )

col.DNase <- colorRampPalette(c("#EBEBEB","#5B7EC6","#233964"))(179)
brk.DNase <- seq(min(log10(tags$DNase)),max(log10(tags$DNase)),length=180)
cl.DNase <- ifelse( ( (tags$RNAPII-tags$Input)>0 & (tags$H3K4me1-tags$Input)>0 ), col.DNase [ as.numeric( cut( log10(tags$DNase), breaks=brk.DNase )) ], "black" )

# Biplot shows the DNA occupancy levels of RNAPII and H3K4me1 at accessible pCRMs. pCRMs (dots) are colored according to their chromatin accessibility
# except for the pCRMs that showed RNAPII and/or H3K4me1 signals piled up over a distance of 1 kb below background.
pdf( "dhs_pol2_h3k4me1_st8p.pdf", width=4.5, height=4.5, useDingbats=FALSE )
  plot( log10(tags$RNAPII), log10(tags$H3K4me1), ann=T, axes=T, cex=1.4, pch=20, lwd=.25, col = alpha(cl.DNase, .7) )
dev.off()

# Color scale for chromatin accessibility (log scale)
pdf("chromatin_access_colorbar.pdf", width=2, height=5)
  z <- matrix(1:180,nrow=1)
  x <- 1
  image(x,brk.DNase,z,col=col.DNase,axes=FALSE,xlab="",ylab="")
  axis(2, at=c( brk.DNase[1], brk.DNase[90], brk.DNase[180]), labels=c( brk.DNase[1], brk.DNase[90], brk.DNase[180]))
dev.off()





# ---- Figure 1f: Chromatin accessibility versus H3K4me1 and RNAPII ----

# FUNCTION: plot mean +/- standard deviation of read density (logarithmic scale) across a 2-kb genomic region
plotMetaLogNorm <- function (df, colour, min, max) {
  
  lstd 	<- sapply(log10(df+0.1),sd)
  dist 	<- seq(-950,950,by=25)
  lmean 	<- colMeans(log10(df+0.1))
  plot(dist, lmean, xlim=c(-1000,1000),ylim=c(log10(min),log10(max)),type="n", axes=F, ann=F)
  polygon(c(dist,rev(dist)),c(lmean+lstd,rev(lmean-lstd)),col=alpha( colour, 0.5 ), border = F)
  lines(dist, lmean, col=colour)
  axis( 1, at=c(-1000,-500,0,500,1000), labels=c("-1","","0","","+1") )
  axis( 2, at=c(log10(0.1),log10(1),log10(10)), labels=c("0.1","1","10") )
  
}

occ.matrix   <- read.table( "./xenTro71/supplementary_files/dhs_st8p_tagDensity1K25bp_DNase_pol2_h3k4me1.txt.gz", header=T, row.names="Gene" )
n.profiles   <- as.numeric( unlist( strsplit( colnames( occ.matrix )[ ncol(occ.matrix) ], ".", fixed=TRUE ))[2] ) + 1
n.bins       <- ncol(occ.matrix)/n.profiles
nr           <- nrow(occ.matrix)

pdf("metaDHSst8ptop_DNase_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#4c77bb", 0.1, 20)
dev.off()

pdf("metaDHSst8ptop_RNAPII_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(n.bins+3):(2*n.bins-2)],"#468453", 0.1, 20)
dev.off()

pdf("metaDHSst8ptop_H3K4me1_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(2*n.bins+3):(3*n.bins-2)],"#b8696e", 0.1, 20)
dev.off()





# ---- Figure 3c: Maternal TFs/signal mediators versus RNAPII ----

# Meta-plots summarize the level of RNAPII engagement across 2,000 pCRMs most frequently occupied by the indicated TFs or signal mediators
# at the 1,024-cell and early gastrula stage, respectively.

directory    <- "./xenTro71/supplementary_files/_pol2"
tagDensity   <- grep("top2K_tagDensity2K25bp",list.files(directory),value=T)
colorTable   <- data.frame(tf=c("sox3","foxh1","vegT","smad1","smad2","bcat"), tf_color=c("#ffd600","#70cde6","#346633","#ec008b","#1b75bb","#00a79c"))
splitter.fun <- function(x, symbol = ".", field = 1) {
  y <- strsplit(as.character(x), symbol, fixed=T)
  z <- y[[1]][field]
  return(z)
}

for (i in tagDensity) { 
  fileDir      <- paste0("./xenTro71/supplementary_files/_pol2/",i)
  fileLabel    <- splitter.fun(i,".",1)
  tf.name      <- splitter.fun(i,"_",1)
  tf.color     <- colorTable[which(colorTable$tf==tf.name),"tf_color"]
  occ.matrix   <- read.table( fileDir, header=T, row.names="Gene" )
  n.profiles   <- as.numeric( unlist( strsplit( colnames( occ.matrix )[ ncol(occ.matrix) ], ".", fixed=TRUE ))[2] ) + 1
  n.bins       <- ncol(occ.matrix)/n.profiles
  nr           <- nrow(occ.matrix)
  pdfname      <- paste0("meta",fileLabel,"_logstd.pdf")
  pdf(pdfname, width=2, height=1.66666)
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  par(mfg=c(1,1))
  plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],as.character(tf.color), 0.1, 10)
  par(mfg=c(1,1))
  plotMetaLogNorm(occ.matrix[,(n.bins+3):(2*n.bins-2)],"#468453", 0.1, 10)
  dev.off()
}





# ---- Figure 3f: TFs versus signal mediators (at TF+ pCRMs) ----

# Meta-plots summarizes the level of signal mediator engagement across 2,000 pCRMs most frequently occupied by the indicated TF 
# at the 1,024-cell and early gastrula stage, respectively.

occ.matrix   <- read.table( "./xenTro71/supplementary_files/_signal_mediators/sox3_st8_top2K_tagDensity2K25bp_sox3_bcat_smad1_smad2.txt.gz", header=T, row.names="Gene" )
n.profiles   <- as.numeric( unlist( strsplit( colnames( occ.matrix )[ ncol(occ.matrix) ], ".", fixed=TRUE ))[2] ) + 1
n.bins       <- ncol(occ.matrix)/n.profiles
nr           <- nrow(occ.matrix)

pdf("metaSox3st8top2K_sox3_bcat_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#ffd600", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(n.bins+3):(2*n.bins-2)],"#00a79c", 0.1, 10)
dev.off()

pdf("metaSox3st8top2K_sox3_smad1_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#ffd600", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(2*n.bins+3):(3*n.bins-2)],"#ec008b", 0.1, 10)
dev.off()

pdf("metaSox3st8top2K_sox3_smad2_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#ffd600", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(3*n.bins+3):(4*n.bins-2)],"#1b75bb", 0.1, 10)
dev.off()

occ.matrix   <- read.table( "./xenTro71/supplementary_files/_signal_mediators/foxh1_st8_top2K_tagDensity2K25bp_foxh1_bcat_smad1_smad2.txt.gz", header=T, row.names="Gene" )
n.profiles   <- as.numeric( unlist( strsplit( colnames( occ.matrix )[ ncol(occ.matrix) ], ".", fixed=TRUE ))[2] ) + 1
n.bins       <- ncol(occ.matrix)/n.profiles
nr           <- nrow(occ.matrix)

pdf("metaFoxH1st8top2K_foxH1_bcat_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#70cde6", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(n.bins+3):(2*n.bins-2)],"#00a79c", 0.1, 10)
dev.off()

pdf("metaFoxH1st8top2K_foxH1_smad1_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#70cde6", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(2*n.bins+3):(3*n.bins-2)],"#ec008b", 0.1, 10)
dev.off()

pdf("metaFoxH1st8top2K_foxH1_smad2_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#70cde6", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(3*n.bins+3):(4*n.bins-2)],"#1b75bb", 0.1, 10)
dev.off()

occ.matrix   <- read.table( "./xenTro71/supplementary_files/_signal_mediators/vegT_st8_top2K_tagDensity2K25bp_vegT_bcat_smad1_smad2.txt.gz", header=T, row.names="Gene" )
n.profiles   <- as.numeric( unlist( strsplit( colnames( occ.matrix )[ ncol(occ.matrix) ], ".", fixed=TRUE ))[2] ) + 1
n.bins       <- ncol(occ.matrix)/n.profiles
nr           <- nrow(occ.matrix)

pdf("metaVegTst8top2K_vegT_bcat_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#346633", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(n.bins+3):(2*n.bins-2)],"#00a79c", 0.1, 10)
dev.off()

pdf("metaVegTst8top2K_vegT_smad1_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#346633", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(2*n.bins+3):(3*n.bins-2)],"#ec008b", 0.1, 10)
dev.off()

pdf("metaVegTst8top2K_vegT_smad2_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#346633", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(3*n.bins+3):(4*n.bins-2)],"#1b75bb", 0.1, 10)
dev.off()

occ.matrix   <- read.table( "./xenTro71/supplementary_files/_signal_mediators/sox3_st10_top2K_tagDensity2K25bp_sox3_bcat_smad1_smad2.txt.gz", header=T, row.names="Gene" )
n.profiles   <- as.numeric( unlist( strsplit( colnames( occ.matrix )[ ncol(occ.matrix) ], ".", fixed=TRUE ))[2] ) + 1
n.bins       <- ncol(occ.matrix)/n.profiles
nr           <- nrow(occ.matrix)

pdf("metaSox3st10top2K_sox3_bcat_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#ffd600", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(n.bins+3):(2*n.bins-2)],"#00a79c", 0.1, 10)
dev.off()

pdf("metaSox3st10top2K_sox3_smad1_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#ffd600", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(2*n.bins+3):(3*n.bins-2)],"#ec008b", 0.1, 10)
dev.off()

pdf("metaSox3st10top2K_sox3_smad2_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#ffd600", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(3*n.bins+3):(4*n.bins-2)],"#1b75bb", 0.1, 10)
dev.off()

occ.matrix   <- read.table( "./xenTro71/supplementary_files/_signal_mediators/foxh1_st10_top2K_tagDensity2K25bp_foxh1_bcat_smad1_smad2.txt.gz", header=T, row.names="Gene" )
n.profiles   <- as.numeric( unlist( strsplit( colnames( occ.matrix )[ ncol(occ.matrix) ], ".", fixed=TRUE ))[2] ) + 1
n.bins       <- ncol(occ.matrix)/n.profiles
nr           <- nrow(occ.matrix)

pdf("metaFoxH1st10top2K_foxH1_bcat_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#70cde6", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(n.bins+3):(2*n.bins-2)],"#00a79c", 0.1, 10)
dev.off()

pdf("metaFoxH1st10top2K_foxH1_smad1_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#70cde6", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(2*n.bins+3):(3*n.bins-2)],"#ec008b", 0.1, 10)
dev.off()

pdf("metaFoxH1st10top2K_foxH1_smad2_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#70cde6", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(3*n.bins+3):(4*n.bins-2)],"#1b75bb", 0.1, 10)
dev.off()

occ.matrix   <- read.table( "./xenTro71/supplementary_files/_signal_mediators/vegT_st10_top2K_tagDensity2K25bp_vegT_bcat_smad1_smad2.txt.gz", header=T, row.names="Gene" )
n.profiles   <- as.numeric( unlist( strsplit( colnames( occ.matrix )[ ncol(occ.matrix) ], ".", fixed=TRUE ))[2] ) + 1
n.bins       <- ncol(occ.matrix)/n.profiles
nr           <- nrow(occ.matrix)

pdf("metaVegTst10top2K_vegT_bcat_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#346633", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(n.bins+3):(2*n.bins-2)],"#00a79c", 0.1, 10)
dev.off()

pdf("metaVegTst10top2K_vegT_smad1_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#346633", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(2*n.bins+3):(3*n.bins-2)],"#ec008b", 0.1, 10)
dev.off()

pdf("metaVegTst10top2K_vegT_smad2_logstd.pdf", width=2, height=1.66666)
par(mfrow=c(1,1),mar=c(2,2,2,2))
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,3:(n.bins-2)],"#346633", 0.1, 10)
par(mfg=c(1,1))
plotMetaLogNorm(occ.matrix[,(3*n.bins+3):(4*n.bins-2)],"#1b75bb", 0.1, 10)
dev.off()




# ---- Figure 8g and Supplementary Figures 15-17: Changes to chromatin accessibility ----

dhs <- read.table("./xenTro71/supplementary_files/dhs_fdr_genes_50kb_distanceTSS.bed.gz", header=F)
colnames(dhs) <- c("scaffold","dhs.centre","dhs","dhs.fdr","dhs.level","tss","gene","strand","distance.tss")

zga <- read.csv("./xenTro71/supplementary_files/zga_lof.csv")
colnames(zga)[1] <- "gene"

act.genes <- read.table("./xenTro71/supplementary_files/pol2_st8to8p_rzRNA_av3tpm.bed", header=F)
colnames(act.genes) <- c("scaffold","start","end","gene","dot","strand")

mo <- read.table( "./xenTro71/supplementary_files/dhs_motif_matrix.txt.gz", header=TRUE, row.names="Gene" )
mo_ <- mo[,c(39:43,120:124,201:205)]
mo_$mPS.motif.sum <- rowSums(mo_)
mo_$dhs <- row.names(mo_)

dhs.zga        <- plyr::join( dhs, zga, by="gene", type="inner" )
dhs.zga.act8p  <- plyr::join( dhs.zga, act.genes, by="gene", type="inner" )
dhs.zga.act8p  <- dhs.zga.act8p[,1:29]
dhs.zga.act8p.mo <- plyr::join( dhs.zga.act8p, mo_, by="dhs", type="inner" )

# Settings of Trellis graphics
d1 <- trellis.par.get("dot.line") 
d1$lwd <- 0  ## set line width to 0 
trellis.par.set("dot.line",d1) 
trellis.par.set( 
  axis.line=list(col="black"), 
  axis.text=list(col="black", cex=0.6),
  panel.background=list(col="transparent"), 
  par.xlab.text= list(col="black")	
) 

brk.rna         <- c( seq(-.001,40,length=30), seq(40.001,66.666,length=30), seq(66.667,100,length=30), seq(100.001,150,length=30), seq(150.001,250,length=30), seq(250.001,2500,length=30) )
col.rna         <- colorRampPalette(c("#3953A4", "#AFBADA", "#EBEBEB", "#EBEBEB", "#FDAE61", "#FAA41A"))(179)




# ---- Figure 8g: Effect of mPouV/Sox3 LOF on chromatin accessibility and RNAPII-mediated gene expression ----

# Localization of accessible pCRMs (affected, dot colored in orange to red with FDR decreasing from 10%; and unaffected, grey dot) 
# relative to the zygotic TSSs that are active by the MBT and produce enough RNA transcripts to show significant two-fold reductions
# upon alfa-amanitin injection). Gene loci are sorted by mPouV/Sox3 LOF-induced transcript fold changes as shown in the log-scaled bar
# graph.

subset <- dhs.zga.act8p.mo[which(abs(dhs.zga.act8p.mo$distance.tss) <= 20000),]
subset <- subset[order(-subset$mPSLOF.st8to10.uni.mean),]

# Localization of accessible pCRMs replative to the zygotic TSSs
pdf( "AllZygoticGenes_DHS_mPSLOF.pdf", width=5, height=6, useDingbats=FALSE )
with(subset,
     dotplot(reorder(gene,-mPSLOF.st8to10.uni.mean) ~ distance.tss, 
             groups=gene, 
             scales=list(x=list(rot=0, cex=.5), y=list(rot=0, cex=.5)),
             col.var1=c("#969696","#FEB24C","#FC4E2A","#B10026")[findInterval(-log10(dhs.fdr),c(0,1,5,10,20))],
             cex.var=c(.01,.03,.05,.07,.09)[findInterval(dhs.level,c(3,4,5,6,7,8))],
             panel=function(x,y,col.var1,col.var2,cex.var){
               panel.abline(h=y, col="grey90")
               panel.dotplot(x,y,col=col.var1, cex=cex.var, pch=19)
             }
     ))
dev.off()

mPSLOF.df   <- unique(subset[,c("gene","mPSLOF.st8to10.uni.mean")])
col.bars    <- c("#3953A4", "#AFBADA", "#C4C4C4","#FDAE61", "#FAA41A")[findInterval(mPSLOF.df$mPSLOF.st8to10.uni.mean,c(0,25,66.666,150,400,1000))]

# Log-scaled bar graph
pdf( "AllZygoticGenes_mPSLOF_log2.pdf", width=2, height=6, useDingbats=F )
  barplot( log2(mPSLOF.df$mPSLOF.st8to10.uni.mean/100), col=col.bars, border=NA, horiz=T )
dev.off()

pdf("FDR_colorbar.pdf", width=2, height=5)
z <- matrix(1:5,nrow=1)
x <- 1
image(x,c(0,1,5,10,20),z,col=c("#969696","#FEB24C","#FC4E2A","#B10026"),axes=FALSE,xlab="",ylab="")
axis(2, at=c(0,1,5,10,20), labels=c("0","1","5","10","20"))
dev.off()





# ---- Supplementary Figure 15: mPouV/Sox3-induced chromatin accessibility is required for the expression of Nodal-responsive genes ----

# All genes listed here are active by the MBT and their transcript levels are significantly reduced (≥two-fold; FDR ≤10%) upon alfa-amanitin 
# injection. The plot aligned shows the localization and DNase sensitivity (bubble size) of accessible pCRMs (affected, dot colored in orange 
# to red with FDR decreasing from 10%; and unaffected, grey dot) relative to zygotic TSSs. Gene loci are sorted by mPouV/Sox3 LOF-induced transcript 
# fold changes as shown in the heat map.

subset <- dhs.zga.act8p.mo[which(dhs.zga.act8p.mo$NodalLOF.st8to10.uni.mean <= 66.666666 & abs(dhs.zga.act8p.mo$distance.tss) <= 20000),]
subset <- subset[order(-subset$mPSLOF.st8to10.uni.mean),]

pdf( "NodalLOF_down_DHS_mPSLOF.pdf", width=4, height=20, useDingbats=FALSE )
with(subset,
     dotplot(reorder(gene,-mPSLOF.st8to10.uni.mean) ~ distance.tss, 
             groups=gene, 
             scales=list(x=list(rot=0, cex=.5), y=list(rot=0, cex=.5)),
             col.var1=c("#969696","#FEB24C","#FC4E2A","#B10026")[findInterval(-log10(dhs.fdr),c(0,1,5,10,20))],
             cex.var=c(.1,.3,.5,.7,.9)[findInterval(dhs.level,c(3,4,5,6,7,8))],
             panel=function(x,y,col.var1,col.var2,cex.var){
               panel.abline(h=y, col="grey90")
               panel.abline(v=c(-20000,-10000,0,10000,20000), col="grey90")
               panel.dotplot(x,y,col=col.var1, cex=cex.var, pch=19)
             }
     ))
dev.off()

# Heat map: Transcript levels
mPSLOF.df <- unique(subset[,c("gene","mPSLOF.st8to10.uni.mean","WntLOF.st8to10.uni.mean","NodalLOF.st8to10.uni.mean","BmpLOF.st8to10.uni.mean")])
png( filename="NodalLOF_down_DHS_mPSLOF_RNA.png", width=4, height=nrow(mPSLOF.df) )
  par( mar=c(0,0,0,0) )
  image( t( mPSLOF.df[,2:5] ), breaks = brk.rna, col = col.rna, axes=FALSE )
dev.off()





# ---- Supplementary Figure 16a: mPouV/Sox3-induced chromatin accessibility is required for the expression of Wnt-responsive genes ----

# All genes listed here are active by the MBT and their transcript levels are significantly reduced (≥two-fold; FDR ≤10%) upon alfa-amanitin 
# injection. The plot aligned shows the localization and DNase sensitivity (bubble size) of accessible pCRMs (affected, dot colored in orange 
# to red with FDR decreasing from 10%; and unaffected, grey dot) relative to zygotic TSSs. Gene loci are sorted by mPouV/Sox3 LOF-induced transcript 
# fold changes as shown in the heat map.

subset <- dhs.zga.act8p.mo[which(dhs.zga.act8p.mo$WntLOF.st8to10.uni.mean <= 66.666666 & abs(dhs.zga.act8p.mo$distance.tss) <= 20000),]
subset <- subset[order(-subset$mPSLOF.st8to10.uni.mean),]

pdf( "WntLOF_down_DHS_mPSLOF.pdf", width=4, height=10, useDingbats=FALSE )
with(subset,
     dotplot(reorder(gene,-mPSLOF.st8to10.uni.mean) ~ distance.tss, 
             groups=gene, 
             scales=list(x=list(rot=0, cex=.5), y=list(rot=0, cex=.5)),
             col.var1=c("#969696","#FEB24C","#FC4E2A","#B10026")[findInterval(-log10(dhs.fdr),c(0,1,5,10,20))],
             cex.var=c(.1,.3,.5,.7,.9)[findInterval(dhs.level,c(3,4,5,6,7,8))],
             panel=function(x,y,col.var1,col.var2,cex.var){
               panel.abline(h=y, col="grey90")
               panel.abline(v=c(-20000,-10000,0,10000,20000), col="grey90")
               panel.dotplot(x,y,col=col.var1, cex=cex.var, pch=19)
             }
     ))
dev.off()

mPSLOF.df <- unique(subset[,c("gene","mPSLOF.st8to10.uni.mean","WntLOF.st8to10.uni.mean","NodalLOF.st8to10.uni.mean","BmpLOF.st8to10.uni.mean")])
png( filename="WntLOF_down_DHS_mPSLOF_RNA.png", width=4, height=nrow(mPSLOF.df) )
par( mar=c(0,0,0,0) )
image( t( mPSLOF.df[,2:5] ), breaks = brk.rna, col = col.rna, axes=FALSE )
dev.off()





# ---- Supplementary Figure 16b: mPouV/Sox3-induced chromatin accessibility is required for the expression of Bmp-responsive genes ----

# All genes listed here are active by the MBT and their transcript levels are significantly reduced (≥two-fold; FDR ≤10%) upon alfa-amanitin 
# injection. The plot aligned shows the localization and DNase sensitivity (bubble size) of accessible pCRMs (affected, dot colored in orange 
# to red with FDR decreasing from 10%; and unaffected, grey dot) relative to zygotic TSSs. Gene loci are sorted by mPouV/Sox3 LOF-induced transcript 
# fold changes as shown in the heat map.

subset <- dhs.zga.act8p.mo[which(dhs.zga.act8p.mo$BmpLOF.st8to10.uni.mean <= 66.666666 & abs(dhs.zga.act8p.mo$distance.tss) <= 20000),]
subset <- subset[order(-subset$mPSLOF.st8to10.uni.mean),]

pdf( "BmpLOF_down_DHS_mPSLOF.pdf", width=4, height=3, useDingbats=FALSE )
with(subset,
     dotplot(reorder(gene,-mPSLOF.st8to10.uni.mean) ~ distance.tss, 
             groups=gene, 
             scales=list(x=list(rot=0, cex=.5), y=list(rot=0, cex=.5)),
             col.var1=c("#969696","#FEB24C","#FC4E2A","#B10026")[findInterval(-log10(dhs.fdr),c(0,1,5,10,20))],
             cex.var=c(.1,.3,.5,.7,.9)[findInterval(dhs.level,c(3,4,5,6,7,8))],
             panel=function(x,y,col.var1,col.var2,cex.var){
               panel.abline(h=y, col="grey90")
               panel.abline(v=c(-20000,-10000,0,10000,20000), col="grey90")
               panel.dotplot(x,y,col=col.var1, cex=cex.var, pch=19)
             }
     ))
dev.off()

mPSLOF.df <- unique(subset[,c("gene","mPSLOF.st8to10.uni.mean","WntLOF.st8to10.uni.mean","NodalLOF.st8to10.uni.mean","BmpLOF.st8to10.uni.mean")])
png( filename="BmpLOF_down_DHS_mPSLOF_RNA.png", width=4, height=nrow(mPSLOF.df) )
par( mar=c(0,0,0,0) )
image( t( mPSLOF.df[,2:5] ), breaks = brk.rna, col = col.rna, axes=FALSE )
dev.off()






# ---- Supplementary Figure 17: mPouV/Sox3-induced chromatin accessibility is required for the expression of signal non-responsive genes ----

# All genes listed here are active by the MBT and their transcript levels are significantly reduced (≥two-fold; FDR ≤10%) upon alfa-amanitin 
# injection. The plot aligned shows the localization of accessible pCRMs (affected, dot colored in orange to red with FDR decreasing from 10%; 
# and unaffected, grey dot) relative to zygotic TSSs. Gene loci are sorted by mPouV/Sox3 LOF-induced transcript fold changes as shown in the heat map.

subset <- dhs.zga.act8p.mo[which(dhs.zga.act8p.mo$WntLOF.st8to10.uni.mean > 66.666666 & 
                                   dhs.zga.act8p.mo$NodalLOF.st8to10.uni.mean > 66.666666 & 
                                   dhs.zga.act8p.mo$BmpLOF.st8to10.uni.mean > 66.666666 & 
                                   abs(dhs.zga.act8p.mo$distance.tss) <= 20000),]
subset <- subset[order(-subset$mPSLOF.st8to10.uni.mean),]

pdf( "SignalLOF_NOTdown_DHS_mPSLOF.pdf", width=4, height=10, useDingbats=FALSE )
with(subset,
     dotplot(reorder(gene,-mPSLOF.st8to10.uni.mean) ~ distance.tss, 
             groups=gene, 
             scales=list(x=list(rot=0, cex=.5), y=list(rot=0, cex=.5)),
             col.var1=c("#969696","#FEB24C","#FC4E2A","#B10026")[findInterval(-log10(dhs.fdr),c(0,1,5,10,20))],
             cex.var=c(.15,.2,.25,.3,.35)[findInterval(dhs.level,c(3,4,5,6,7,8))],
             panel=function(x,y,col.var1,col.var2,cex.var){
               panel.abline(h=y, col="grey90")
               panel.abline(v=c(-20000,-10000,0,10000,20000), col="grey90")
               panel.dotplot(x,y,col=col.var1, cex=cex.var, pch=19)
             }
     ))
dev.off()

mPSLOF.df <- unique(subset[,c("gene","mPSLOF.st8to10.uni.mean","WntLOF.st8to10.uni.mean","NodalLOF.st8to10.uni.mean","BmpLOF.st8to10.uni.mean")])
col.bars <- c("#3953A4", "#AFBADA", "#C4C4C4","#FDAE61", "#FAA41A")[findInterval(mPSLOF.df$mPSLOF.st8to10.uni.mean,c(0,25,66.666,150,400,1000))]

mPSLOF.df <- unique(subset[,c("gene","mPSLOF.st8to10.uni.mean","WntLOF.st8to10.uni.mean","NodalLOF.st8to10.uni.mean","BmpLOF.st8to10.uni.mean")])
png( filename="SignalLOF_NOTdown_mPSLOF.png", width=4, height=nrow(mPSLOF.df) )
par( mar=c(0,0,0,0) )
image( t( mPSLOF.df[,2:5] ), breaks = brk.rna, col = col.rna, axes=FALSE )
dev.off()

pdf( "SignalLOF_NOTdown_DHS_mPSLOF_highRes.pdf", width=4, height=120, useDingbats=FALSE )
with(subset,
     dotplot(reorder(gene,-mPSLOF.st8to10.uni.mean) ~ distance.tss, 
             groups=gene, 
             scales=list(x=list(rot=0, cex=.5), y=list(rot=0, cex=.5)),
             col.var1=c("#969696","#FEB24C","#FC4E2A","#B10026")[findInterval(-log10(dhs.fdr),c(0,1,5,10,20))],
             cex.var=c(.1,.3,.5,.7,.9)[findInterval(dhs.level,c(3,4,5,6,7,8))],
             panel=function(x,y,col.var1,col.var2,cex.var){
               panel.abline(h=y, col="grey90")
               panel.abline(v=c(-20000,-10000,0,10000,20000), col="grey90")
               panel.dotplot(x,y,col=col.var1, cex=cex.var, pch=19)
             }
     ))
dev.off()
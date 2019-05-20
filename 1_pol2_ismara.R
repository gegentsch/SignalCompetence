# DETECTION OF EARLY RNAPII-ENGAGED CIS-REGULATORY MODULES DURING ZYGOTIC GENOME ACTIVATION
# Author: G.E. Gentsch
#
#
# RNAPII ChIP-Seq reads were aligned and peaks detected as described in the Methods section of the paper.
# ISMARA was run in expert mode
# https://ismara.unibas.ch/fcgi/mara
# HOMER motif matrices (to find enriched motifs in the genome; germ_layer_motifs.txt) and ISMARA input files are provided under ./xenTro71/supplementary_files:
# :: pol2_motifCount_ismara.txt.gz (motif site count table), pol2_tagCount_ismara.txt.gz (expression profile)
#
# BAM files and input csv file are required for running the second part of the script as described in the vignette of R package DiffBind.
# Uncomment once you have the files ready. Sample csv file (pol2.csv) and BED file containing RNAPII CRM binding coordinates are provided under ./xenTro71/supplementary_files. The BAM files need to be generated de novo from the raw FASTQ files deposited to GEO archive: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113186
#
# Installing (if required) and uploading several R packages
required.pkg <- c( "seriation", "extrafont", "ggplot2", "DiffBind" )
#
for ( pkg in required.pkg ) {
    if ( pkg %in% rownames(installed.packages()) == FALSE) {
        install.packages( pkg ) }
    if ( pkg %in% rownames(.packages()) == FALSE ) {
        library(pkg, character.only = TRUE) }
}

#
# Auxiliary functions
source( "./R_utils/colour_schemes.R" )
#
# load tables from ISMARA output
mota                <- t( read.table( "./xenTro71/pol2_ismara_activity.txt", header=TRUE ) )
mota.z              <- read.table( "./xenTro71/pol2_ismara_matrices.txt", header=FALSE )
colnames(mota.z)    <- c( "motif", "zval" )

# load tables from HOMER output
p6                  <- read.table( "./xenTro71/pol2_homer_st6to8_motifs.txt", header=TRUE )
p9                  <- read.table( "./xenTro71/pol2_homer_st9to12_motifs.txt", header=TRUE )

cln                 <- c("consensus","pval","log.pval","fdr","occ.motif.top2K","frac.pos.peak","occ.motif.bkg","frac.pos.bkg")
colnames(p6)        <- cln
colnames(p9)        <- cln


convert <- function(df,prefix) {
    df <- within(df, {
        frac.pos.peak   <- as.character( frac.pos.peak )
        frac.pos.bkg    <- as.character( frac.pos.bkg )
        rad.peak        <- sqrt( as.numeric( substr( frac.pos.peak, 1, nchar(frac.pos.peak ) - 1 ) ) / pi )
        rad.bkg         <- sqrt( as.numeric( substr( frac.pos.bkg, 1, nchar(frac.pos.bkg ) - 1 ) ) / pi )
        motif.enriched  <- as.numeric( substr( frac.pos.peak, 1, nchar( frac.pos.peak ) - 1 ) ) / as.numeric( substr( frac.pos.bkg, 1, nchar( frac.pos.bkg ) - 1 ) ) }
    )
    colnames(df)[ 2:11 ] <- paste( prefix, colnames(df)[2:11], sep="." )
    return(df)
}

p6              <- convert(p6,"p6to8")
p9              <- convert(p9,"p9to12")

m               <- merge( mota, mota.z, by.x="row.names", by.y="motif", all=TRUE )
p               <- merge( p6, p9, by="row.names", all=TRUE )
mp              <- merge( m, p, by="Row.names", all=TRUE )
row.names(mp)   <- mp$Row.names
mp              <- mp[-1]
row.names(p)    <- p$Row.names
p               <- p[-1]

mp.fdr          <- mp[ which( mp$p6to8.fdr <= 0.01 | mp$p9to12.fdr <= 0.01 | row.names(mp)=="CCAAT.NFY" | row.names(mp)=="TATA" ), ]
mp.fdr          <- mp.fdr[-which(row.names(mp.fdr)=="T.T"),]

distfunc        <- function(x) dist(x,method="euclidean")

dm              <- distfunc(mp.fdr[,1:6])

s               <- seriate(dm, method="ARSA")
s.order         <- get_order(s)

pdf( "pol2_ismara_homer_motif_ARSA.pdf",  width=5, height=5 )
    par(mfrow=c(1,2))
    layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
    par(mfg=c(1,1),mar=c(0,0,0,15),las=1)
    image(data.matrix(t(mp.fdr[s.order,1:6])), breaks = brk.z, col = col.z, axes=F)
    axis(2, at=seq(0,1,by=1/(nrow(mp.fdr)-1)), labels=row.names(mp.fdr[s.order,]),side=4)
    par(mfg=c(1,2),mar=c(0,0,0,0))
    image(data.matrix(t(log10(-1*mp.fdr[s.order,c(10,21)]))), breaks = brk.p.motif, col = col.p.motif, axes=F)
dev.off()

select.p <- p[row.names(mp.fdr[s.order,]),]


# bubble plot (peak vs background)
rad.bkg         <- select.p[, grep("rad.bkg", colnames(select.p)) ]
rad.peak        <- select.p[, grep("rad.peak", colnames(select.p)) ]

n.cells         <- ( ncol( select.p ) - 2 ) / 10 * nrow( select.p )
tf.stage        <- substr( colnames(rad.peak), 1, nchar( colnames(rad.peak)) - 9 )

grid            <- data.frame( matrix(0, ncol=4, nrow=(2*n.cells)) )
grid$X1         <- rep( tf.stage, each=nrow(select.p) )
grid$X2         <- row.names(rad.peak)
grid$X3         <- unlist( list( rad.peak, rad.bkg ) )
grid$X4         <- rep( c("peak","background"), each=n.cells )

# factorise X1 (TFs) and X2 (motifs) to ensure order does not change when calling ggplot()
grid$X5         <- factor( grid$X1, tf.stage )
grid$X6         <- factor( grid$X2, rev(row.names(select.p)) )

pdf( file="pol2_homer_peakVsBgd.pdf", width=9, height=9, useDingbats=FALSE, family="Arial" )
    ggplot( grid, aes(X5,X6,color=X4 ) ) +
    geom_point( aes( size = X3^2*pi ) ) +
    scale_size_continuous( range = range( c( 0.001, 15 ) ) ) +
    scale_color_manual( values = c("black","skyblue1") ) +
    theme_bw() +
    theme( panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black") )
dev.off()

pdf("Zscore_colorbar.pdf", width=2, height=5)
z <- matrix(1:100,nrow=1)
x <- 1
image(x,brk.z,z,col=col.z,axes=FALSE,xlab="",ylab="")
axis(2, at=c(-8,-6,-4,-2,0,2,4,6,8), labels=c("-8","","-4","","0","","+4","","+8"))
dev.off()

pdf("pval_colorbar.pdf", width=2, height=5)
z <- matrix(1:length(brk.p.motif),nrow=1)
x <- 1
image(x,brk.p.motif,z,col=col.p.motif,axes=FALSE,xlab="",ylab="")
axis(2, at=c(log10(5),log10(50),log10(500),log10(3000)), labels=c("5","50","500","3000"))
dev.off()


#
# Dynamics of RNAPII recruitment to CRMs
# --------------------------------------------------------
#
# dba       <- dba( sampleSheet="./xenTro71/supplementary_files/pol2.csv" )
# d.bam     <- dba.count( dba, minOverlap=0, score=DBA_SCORE_READS, bRemoveDuplicates=TRUE, mapQCth=10, fragmentSize=175, summits=250, bParallel=TRUE )
# 
# d.spearman <- dba.plotHeatmap( d.bam, distMethod="spearman" )
#
# pdf("pol2_CRM_spearman.pdf", width=6, height=6)
#   heatmap.2(
#    t(d.spearman), Rowv=F,Colv=T, dendrogram=c("column"),trace="none", symm=T,scale="none", breaks=br.corr, col=col.corr, key=TRUE, symkey=FALSE, density.info="none", cellnote=round(d.spearman,1), notecol="black", notecex=1.2, cexRow=1.2, cexCol=1.2, margin=c(13,13) )
# dev.off()


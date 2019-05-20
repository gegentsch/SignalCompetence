# ENDOGENOUS BINDING DYNAMICS OF TRANSCRIPTION FACTORS AND SIGNALLING MEDIATORS
# Author: G.E. Gentsch
#
#
#
# ChIP-Seq reads were aligned and peaks detected as described in the Methods section of the paper.
#
# The input templates (.csv) for the DiffBind R package and the genome coordinates in BED format are provided under directory ./xenTro71/supplementary_files. The BAM files need to be generated de novo from raw FASTQ files deposited to GEO archive:
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113186
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48560 (VegT and Eomes, stage 12+; and Xbra, stage 12+ and 20)
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85273 (FoxH1, stage 8)
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53654 (FoxH1 and biological replicate for Smad2, stage 10+)
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72657 (biological replicate for b-catenin, stage 10+)
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30146 (biological replicate for Smad2, stage 10+)

# Installing (if required) and uploading several R packages
required.pkg <- c( "seriation", "DESeq2", "extrafont", "gplots", "DiffBind" )

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




# Endogenous dynamics of TF binding
# -----------------------------------------------------------------------------------------

# The genomic coordinates were selected based on:
# (1)   Top 2,000 peaks for each TF or signalling mediator (FDR <= 0.1%)
# (2)   >= 20 normalised tags per ChIP-Seq peak
# (3)   NOT USED: VegT at stage 10+ (detection of both maternal and zygotic isoform in endoderm und nascent mesoderm, respectively)
# (4)   NOT USED: Spatial binding profiles of Sox3 (separate experiment)
# (5)   NOT USED: b-catenin at stage 6 and 8 because of low DNA occupancy levels, which more frequently caused incorrect peak calling!
# This BED file is provided under ./xenTro71/supplementary_files: tf_top2K_peakScoreMin20.bed
# As a control, pooled input (stage 6 to 20) was used for all ChIP-Seq samples.


# all (maternal and zygotic) TFs and signal mediators
tf.dba            <- dba( sampleSheet = "./xenTro71/supplementary_files/tf.csv" )

# RNAPII, Sox3, Smad1 and Smad2
pol2.dba          <- dba( sampleSheet = "./xenTro71/supplementary_files/tf_pol2.csv")

tf.tme            <- dba.count( tf.dba, minOverlap=0, score=DBA_SCORE_TMM_MINUS_EFFECTIVE, bRemoveDuplicates=TRUE, fragmentSize=175 )
pol2.tme          <- dba.count( pol2.dba, minOverlap=0, score=DBA_SCORE_TMM_MINUS_EFFECTIVE, bRemoveDuplicates=TRUE, fragmentSize=175 )


# Correlation matrix
tf.tme.spearman   <- dba.plotHeatmap( tf.tme, distMethod="spearman", margin=10 )

write.csv( tf.tme.spearman, "tf_spearman.csv" )
# swap branches so that temporal order is upheld. Modified file provided under ./xenTro71/supplementary_files.


ms                <- as.matrix( read.csv( "./xenTro71/supplementary_files/tf_spearman.csv", header=TRUE, row.names="X" ) )

pdf( "tf_spearman.pdf", width=40, height=40, useDingbats=FALSE, family="Arial" )
heatmap.2( ms, Rowv=F, Colv=F, dendrogram=c( "none" ), trace="none", symm=TRUE, scale="none", breaks=brk.corr, col=col.corr, key=TRUE, symkey=TRUE, density.info="none", cellnote=round( ms, 1 ), notecol="black", notecex=5, cexRow=5, cexCol=5, margin=c(40,40), lwid=c(1,10), lhei=c(1,10) )
dev.off()

pdf( "tf_spearman_noCellNote.pdf", width=40, height=40, useDingbats=FALSE, family="Arial" )
heatmap.2( ms, Rowv=F, Colv=F, dendrogram=c( "none" ), trace="none", symm=TRUE, scale="none", breaks=brk.corr, col=col.corr, key=TRUE, symkey=TRUE, density.info="none", cexRow=5, cexCol=5, margin=c(40,40), lwid=c(1,10), lhei=c(1,10) )
dev.off()



# PCA plots
pdf( file="tf_pca.pdf", width=10, height=10, useDingbats=FALSE, family="Arial" )
dba.plotPCA( tf.tme, DBA_TISSUE, label=DBA_ID )
dev.off()

pdf( file="pol2_tf_pca.pdf", width=10, height=10, useDingbats=FALSE, family="Arial" )
dba.plotPCA( pol2.tme, DBA_TISSUE, label=DBA_ID )
dev.off()



# Colour bar for DNA motif p-values
z           <- matrix(1:100,nrow=1)
pdf("ColourBar_motif_pval.pdf", width=2, height=4)
par(mar=c(0,0,0,0) + 4)
    image(1,brk.p.motif2, z, col=col.p.motif, axes=FALSE, ylim=log( c( 1, 10000 ) ), xlab="", ylab="-log(p)" )
    axis(2, at=log( c( 1, 10, 100, 1000, 10000 ) ), labels=c( 1, 10, 100, 1000, 10000 ) )
dev.off()


# Read DNA motif occurence / enrichment
v.8         <- read.motif( "./xenTro71/supplementary_files/homer_motifs/vegT_st8_motifs.txt", "v.8" )
f.8         <- read.motif( "./xenTro71/supplementary_files/homer_motifs/foxH1_st8_motifs.txt", "f.8" )
s2.8        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/smad2_st8_motifs.txt", "s2.8" )
s.8         <- read.motif( "./xenTro71/supplementary_files/homer_motifs/sox3_st8_motifs.txt", "s.8" )
s1.8        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/smad1_st8p_motifs.txt", "s1.8" )
f.10        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/foxH1_st10p_motifs.txt", "f.10" )
s2.10       <- read.motif( "./xenTro71/supplementary_files/homer_motifs/smad2_st10p_motifs.txt", "s2.10" )
s1.10       <- read.motif( "./xenTro71/supplementary_files/homer_motifs/smad1_st10p_motifs.txt", "s1.10" )
b.10        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/bcat_st10p_motifs.txt", "b.10" )
s.10        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/sox3_st10p_motifs.txt", "s.10" )
s.12        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/sox3_st12p_motifs.txt", "s.12" )
e.12        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/eomes_st12p_motifs.txt", "e.12" )
v.12        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/vegT_st12p_motifs.txt", "v.12" )
x.12        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/xbra_st12p_motifs.txt", "x.12" )
x.20        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/xbra_st20_motifs.txt", "x.20" )
t.20        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/tbx6_st20_motifs.txt", "t.20" )

list.tf     <- list(v.8, f.8, s2.8, s.8, s1.8, f.10, s2.10, s1.10, b.10, s.10, s.12, e.12, v.12, x.12, x.20, t.20)
all.tf      <- plyr::join_all( list.tf )
row.names(all.tf) <- row.names( v.8 )

select.tf <- all.tf[ c( "SOX.sox3", "POU.SOX.pou5f3.sox3", "POU.pou5f3", "FOXH.foxH1", "pairedHD.otx", "T.box", "T.T", "E.box.grA", "SMAD.smad12", "bHSH.tfap2"), ]


grid <- plot.motif( dataframe=select.tf, exp.name="tf" )


pdf( file="tf_motifs_peakVsBgd.pdf", width=9, height=9, useDingbats=FALSE, family="Arial" )
    ggplot( grid, aes(X5,X6,color=X4)) +
    geom_point(aes(size=X3^2*pi)) +
    scale_size_continuous(range = range( c( 0.001, 15 ))) +
    scale_color_manual(values = c("black","skyblue1")) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()




# The effectx of ectopic expression of MyoD-HA on chromatin recruitment of Sox3 and RNAPII
# -------------------------------------------------------------------------------------------------------------------

# Ectopic overexpression of MyoD-HA in animal hemisphere
# Sox3/RNAPII from uninjected embryos at stage 8, 10+, 12+  vs. Sox3/RNAPII/MyoD from MyoD-HA-injected embryos at stage 10.
# Top 10,000 peaks of every ChIP-Seq sample merged. This BED file is provided under ./xenTro71/supplementary_files: myod_peaksScoreMin20.bed


m.dba           <- dba( sampleSheet="./xenTro71/supplementary_files/myod.csv" )
m.dba           <- dba.count( m.dba, minOverlap=0, score=DBA_SCORE_TMM_MINUS_EFFECTIVE, bRemoveDuplicates=TRUE, mapQCth=10, fragmentSize=175, summits=250, bParallel=TRUE )


# Spearman correlation matrix
pdf( file="myoD_spearman.pdf", width=10, height=10, useDingbats=FALSE, family="Arial" )
m.tme.spearman  <- dba.plotHeatmap( m.dba, distMethod="spearman", margin=10 )
dev.off()


# PCA plot
pdf( file="myoD_pca.pdf", width=6, height=6, useDingbats=FALSE, family="Arial" )
dba.plotPCA( m.dba.merge, attributes=c( DBA_TISSUE, DBA_TREATMENT ), label=DBA_ID )
dev.off()


# reordering without breaking linkage
reorder         <- c( 6, 5, 9, 4, 8, 3, 1, 7, 2)
m.dba.alt       <- m.tme.spearman[ reorder, reorder ]


# plot re-ordered correlation matrix
print.cor.pdf( filename="myod_spearman.pdf", data.matrix=m.dba.alt )
print.cor.pdf( filename="myod_spearman_noCellNote.pdf", data.matrix=m.dba.alt, display.cor=FALSE )


# DNA binding and motif matrices are provided under ./xenTro71/supplementary_files
win.sox3        <- read.table( "./xenTro71/supplementary_files/sox3_myodOE_matrix.txt.gz", header=TRUE )
win.pol2        <- read.table( "./xenTro71/supplementary_files/pol2_myodOE_matrix.txt.gz", header=TRUE )
sox3.pval       <- pValueMatrix( win.sox3 )
pol2.pval       <- pValueMatrix( win.pol2 )


# DNA motif occurence +/-1 kb from peak center (25 bp bins)
mo.vis          <- read.table("./xenTro71/supplementary_files/myod_motifs_matrix.txt.gz", header=TRUE, row.names="Gene")


# DNA occupancy matrix +/-1 kb from peak center (25 bp bins)
tf              <- read.table("./xenTro71/supplementary_files/sox3_pol2_myod_matrix.txt.gz", header=TRUE, row.names="Gene")


# Calculate average of top 100 peaks (choose 25bp-bin with highest read count) to normalise profiles
n.profiles      <- as.numeric(unlist(strsplit(colnames(tf)[ncol(tf)],".",fixed=T))[2])+1
n.bins          <- ncol(tf)/n.profiles
n.mid           <- ceiling(n.bins/2)
tf.scale        <- tf


# Scaling DNA occupancy
for (i in seq(1,ncol(tf),n.bins)) {
    
    tf.scale[ , i:( i + n.bins - 1 ) ] <- 100 / mean( sort( apply( tf[ , i:( i + n.bins - 1 ) ], 1, max ), decreasing=T )[ 1:100 ] ) * tf[ ,i:( i + n.bins - 1 ) ]

}


# DNA occupancy
brk.uni         <- c( seq( 0, 50, length=100 ), seq( 50.5, max(tf.scale), length=100 ) )
col.uni         <- colorRampPalette( c( "#EBEBEB", "black", "#FF0000", "#FFFFCC" ) )(199)
# DNA motif occurence
brk.mot         <- c( 1, max(mo.vis) )
col.mot         <- "black"
# Differential binding (p-value)
brk.p           <- c( seq(0,8,length=90),seq(9,max(sox3.pval),length=10) )
col.p           <- colorRampPalette(c("#FFFFFF","#FF0000"))(99)


png(file="myoD_heatmap.png", width=n.bins*11, height=nrow(tf))
    par(mar=c(0,0,0,0))
    layout( matrix(c(rep(1,1),rep(2,2),rep(3,1),rep(4,2),rep(5,1),rep(6,4)), 1, 11, byrow = TRUE ) )
    image( t( tf.scale[ , 1:n.bins ] ), breaks=brk.uni, col=col.uni, axes=FALSE ) #1
    image( t( tf.scale[,(n.bins+1):(3*n.bins)] ), breaks=brk.uni, col=col.uni, axes=FALSE ) #2
    image( t( sox3.pval ), breaks=brk.p, col=col.p, axes=FALSE ) #3
    image( t( tf.scale[,(3*n.bins+1):(5*n.bins)] ), breaks=brk.uni, col=col.uni, axes=FALSE ) #4
    image( t( pol2.pval ), breaks=brk.p, col=col.p, axes=FALSE ) #5
    image( t( mo.vis ), breaks=brk.mot, col=col.mot, axes=FALSE )  #6
dev.off()


# Colour bars

# scaled DNA occupancy
z           <- matrix( 1:200, nrow=1 )
pdf("ColourBar_DNAocc_scaled.pdf", width=2, height=4)
par(mar=c(0,0,0,0) + 4)
image( 1, brk.uni, z, col=col.uni, axes=FALSE, ylim=c(0, max(tf.scale) ), xlab="", ylab="DNA occupancy (scaled)")
axis(2, at=c(0,50,100,150,200), labels=c( "0", "0.5", "1", "1.5", "2" ) )
dev.off()

# p-values
z           <- matrix( 1:100, nrow=1 )
pdf( "ColourBar_diffPval.pdf", width=2, height=4)
par(mar=c(0,0,0,0) + 4)
image( 1, brk.p, z, col=col.p, axes=FALSE, ylim=c(0,10), xlab="", ylab="-log(p)")
axis( 2, at=c( 0, 5, 10 ), labels=c( "0", "5", "10" ) )
dev.off()



# Reading DNA motif occurence / enrichment
s.8         <- read.motif( "./xenTro71/supplementary_files/homer_motifs/sox3_st8_top10K_motifs.txt", "s.8" )
s.10        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/sox3_st10p_top10K_motifs.txt", "s.10" )
s.12        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/sox3_st12p_top10K_motifs.txt", "s.12" )
p.8         <- read.motif( "./xenTro71/supplementary_files/homer_motifs/pol2_st8_top10K_motifs.txt", "p.8" )
p.11        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/pol2_st11_top10K_motifs.txt", "p.11" )
p.12        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/pol2_st12_top10K_motifs.txt", "p.12" )
s.10.m      <- read.motif( "./xenTro71/supplementary_files/homer_motifs/sox3_st10p_myodOE_top10K_motifs.txt", "s.10.m" )
m.10.m      <- read.motif( "./xenTro71/supplementary_files/homer_motifs/myod_st10p_myodOE_top10K_motifs.txt", "m.10.m" )
p.10.m      <- read.motif( "./xenTro71/supplementary_files/homer_motifs/pol2_st10p_myodOE_top10K_motifs.txt", "p.10.m" )

list.tf     <- list(s.8, p.8, s.10, p.11, s.12, p.12, m.10.m, s.10.m, p.10.m )
all.tf      <- plyr::join_all(list.tf)

row.names(all.tf) <- row.names(s.8)

select.tf   <- all.tf[c("SOX.sox3","FOXH.foxH1","E.box.grA"),]

grid        <- plot.motif( dataframe=select.tf, exp.name="myod" )

pdf( file="myod_motifs_peakVsBgd.pdf", width=9, height=9, useDingbats=FALSE, family="Arial" )
ggplot( grid, aes(X5,X6,color=X4)) +
geom_point(aes(size=X3^2*pi)) +
scale_size_continuous(range = range( c( 0.001, 15 ))) +
scale_color_manual(values = c("black","skyblue1")) +
theme_bw() +
theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()



# Sox3 chromatin engagement along the anterior-posterior axis: head, trunk and tailbud
# -------------------------------------------------------------------------------------------------------------------

sox3.ap             <- dba( sampleSheet="./xenTro71/supplementary_files/sox3_ap.csv" )
sox3.ap             <- dba.count( sox3.ap, minOverlap=1, score=DBA_SCORE_TMM_MINUS_EFFECTIVE, bRemoveDuplicates=TRUE, mapQCth=10, fragmentSize=175, summits=250, bParallel=TRUE )


# Spearman correlation matrix
sox3.ap.spearman <- dba.plotHeatmap( sox3.ap, distMethod="spearman", margin=10 )


# PCA plots
pdf( file="sox3_ap_pca.pdf", width=6, height=6, useDingbats=FALSE, family="Arial" )
    dba.plotPCA( sox3.ap, DBA_TISSUE, label=DBA_ID )
dev.off()


# Reordering
reorder              <- c(5,1,4:2,8:6)
sox3.ap.spearman.alt <- sox3.ap.spearman[ reorder, reorder ]


# Plot re-ordered correlation matrix
print.cor.pdf( filename="sox3_ap_spearman.pdf", data.matrix=sox3.ap.spearman.alt )
print.cor.pdf( filename="sox3_ap_spearman_noCellNote.pdf", data.matrix=sox3.ap.spearman.alt, display.cor=FALSE )


# DNA motif occurence +/-1 kb from peak center (25 bp bins)
mo.vis               <- read.table( "./xenTro71/supplementary_files/sox3_ap_motifs_matrix.txt.gz", header=TRUE, row.names="Gene")


# DNA occupancy matrix +/-1 kb from peak center (25 bp bins)
tf                   <- read.table( "./xenTro71/supplementary_files/sox3_ap_matrix.txt.gz", header=TRUE, row.names="Gene")


# Calculate average of top 100 peaks (choose 25bp-bin with highest read count) to normalise profiles
n.profiles           <- as.numeric( unlist( strsplit( colnames(tf)[ ncol( tf ) ],".",fixed=TRUE ) )[2]) + 1
n.bins               <- ncol( tf ) / n.profiles
n.mid                <- ceiling( n.bins / 2 )
tf.scale             <- tf


# Scaling DNA occupancy
for (i in seq(1,ncol(tf),n.bins)) {
    
    tf.scale[ , i:( i + n.bins - 1 ) ] <- 100 / mean( sort( apply( tf[ , i:( i + n.bins - 1 ) ], 1, max ), decreasing=T )[ 1:100 ] ) * tf[ ,i:( i + n.bins - 1 ) ]
    
}


# DNA occupancy
brk.uni         <- c( seq( 0, 50, length=100 ), seq( 50.5, max(tf.scale), length=100 ) )
col.uni         <- colorRampPalette( c( "#EBEBEB", "black", "#FF0000", "#FFFFCC" ) )(199)
# DNA motif occurence
brk.mot         <- c( 1, max(mo.vis) )
col.mot         <- "black"
# Differential binding (p-value)
brk.p           <- c( seq(0,8,length=90),seq(9,max(sox3.pval),length=10) )
col.p           <- colorRampPalette(c("#FFFFFF","#FF0000"))(99)



# Use average +/-50 bp from peak center for clustering
# Empty matrix to populate with DNA occupancy levels scaled to 100 (see below)
tf.matr             <- matrix( , nrow=nrow(tf), ncol=n.profiles )
colnames(tf.matr)   <- c( "sox3_st12p", "eomes_st12p", "xbra_st12p","xbra_st20","tbx6_st20","sox3_st20_head","sox3_st20_trunk","sox3_st20_tailbud")
row.names(tf.matr)  <- row.names(tf)



# calculate mean DNA occupancy level (+/- 50 bp from peak center), normalised to mean of top 1000 peaks
j               <- 1
for (i in seq(1,ncol(tf),n.bins)) {
    tf.matr[,j]     <- apply( tf.scale[, (i+n.mid-3):(i+n.mid+1)], 1, mean )
    tf.matr[,j]     <- 100 / mean( sort( tf.matr[ ,j ], decreasing=TRUE )[1:100]) * tf.matr[,j]
    j <- j+1
}

d                   <- dist(tf.matr, method="euclidean")
s.d                 <- seriate(d, method="HC_complete")
s.d.order           <- get_order(s.d)

png(file="sox3_ap_heatmap.png", width=n.bins*15, height=nrow(tf))
par( mar=c(0,0,0,0) )
layout( matrix( c(rep(1,1),rep(2,1),rep(3,2),rep(4,1),rep(5,3),rep(6,7)), 1, 15, byrow=TRUE ))
image( t( tf.scale[ s.d.order, 1:n.bins]), breaks=brk.uni, col=col.uni, axes=FALSE )
image( t( tf.scale[ s.d.order, (n.bins+1):(2*n.bins)]), breaks=brk.uni, col=col.uni, axes=FALSE )
image( t( tf.scale[ s.d.order, (2*n.bins+1):(4*n.bins)]), breaks=brk.uni, col=col.uni, axes=FALSE )
image( t( tf.scale[ s.d.order, (4*n.bins+1):(5*n.bins)]), breaks=brk.uni, col=col.uni, axes=FALSE )
image( t( tf.scale[ s.d.order, c((7*n.bins+1):(8*n.bins),(6*n.bins+1):(7*n.bins),(5*n.bins+1):(6*n.bins)) ]), breaks=brk.uni, col=col.uni, axes=FALSE )
image( t( mo.vis[s.d.order,]), breaks=brk.mot, col=col.mot, axes=FALSE )
dev.off()

# Read DNA motif occurences / enrichment
s.12        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/sox3_st12p_top10K_motifs.txt", "s.12" )
e.12        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/eomes_st12p_top10K_motifs.txt","e.12")
x.12        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/xbra_st12p_top10K_motifs.txt", "x.12")
x.20        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/xbra_st20_top10K_motifs.txt", "x.20")
t.20        <- read.motif( "./xenTro71/supplementary_files/homer_motifs/tbx6_st20_top10K_motifs.txt", "t.20")
s.20.bud    <- read.motif( "./xenTro71/supplementary_files/homer_motifs/sox3_st20_tailbud_top10K_motifs.txt","s.20.bud")
s.20.trunk  <- read.motif( "./xenTro71/supplementary_files/homer_motifs/sox3_st20_trunk_top10K_motifs.txt","s.20.trunk")
s.20.head   <- read.motif( "./xenTro71/supplementary_files/homer_motifs/sox3_st20_head_top10K_motifs.txt","s.20.head")

list.tf     <- list( s.12, e.12, x.12, x.20, t.20, s.20.bud, s.20.trunk, s.20.head )
all.tf      <- plyr::join_all( list.tf )

row.names(all.tf) <- row.names(s.12)

select.tf   <- all.tf[c("SOX.sox3","POU.SOX.pou5f3.sox3","POU.pou5f3","FOXH.foxH1","pairedHD.otx","T.box","T.T","E.box.grA","antpHD.cdx"),]

grid        <- plot.motif( dataframe=select.tf, exp.name="sox3_ap" )

pdf( file="sox3_ap_motifs_peakVsBgd.pdf", width=9, height=9, useDingbats=FALSE, family="Arial" )
    ggplot( grid, aes(X5,X6,color=X4)) +
    geom_point( aes(size=X3^2*pi) ) +
    scale_size_continuous( range = range( c( 0.001, 15 ))) +
    scale_color_manual( values = c("black","skyblue1")) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()


# ----------------------------------------------------------------




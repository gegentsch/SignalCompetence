# ACCESSIBLE CIS-REGULATORY MODULES (CRM) AT MBT
# Author: G.E. Gentsch
#
#
#
# DNase (accessibility), H3K4me1 and RNAPII reads were aligned as described in the Methods section of the paper.
#
# Matrix files (+/-1kb from peak center, 25 bp bins) were generated using HOMER software. The BED (coordinates of accessible conserved and non-conserved CRMs) used to generate these files are provided under ./xenTro91/supplementary_files/
# :: MBT_accessible_cons.bed, MBT_accessible_NOTconserved.bed
#
# Installing (if required) and uploading several R packages
required.pkg <- c( "extrafont", "ggplot2" )
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
# matrix +/-1kb from peak center (25 bp bins): phastCons & phyloP
p1              <- read.table( "./xenTro91/MBT_cons_phastCons_phyloP_matrix.txt.gz", header=TRUE )
p2              <- read.table( "./xenTro91/MBT_NOTcons_phastCons_phyloP_matrix.txt.gz", header=TRUE )
tf.p            <- rbind( p1, p2 )

# matrix +/-1kb from peak center (25 bp bins): DNase, H3K4me1, RNAPII and DNase input
p3              <- read.table("./xenTro91/MBT_cons_dnase_h3k4me1_pol2_input_matrix.txt.gz", header=TRUE )
p4              <- read.table("./xenTro91/MBT_NOTcons_dnase_h3k4me1_pol2_input_matrix.txt.gz", header=TRUE )
tf.ch           <- rbind( p3, p4 )

# 'right' join: keeps order of first table
tf              <- plyr::join( tf.ch, tf.p, by="Gene", type="right" )
row.names(tf)   <- tf$Gene
tf              <- tf[-1]
tf[is.na(tf)]   <- 0

n.profiles      <- 6
n.bins          <- ncol(tf)/n.profiles

# calculate average of top 1000 peaks (choose 25bp-bin with highest read count) to normalise profiles
tf.scale        <- tf

# differential chromatin
for (i in seq(1,ncol(tf),n.bins)) {
    
    tf.scale[ , i:( i + n.bins - 1 ) ] <- 100 / mean( sort( apply( tf[ , i:( i + n.bins - 1 ) ], 1, max ), decreasing=TRUE )[1:1000] ) * tf[ , i : ( i + n.bins - 1 ) ]

}

# Change settings for RNAPII (mean of top 1000 divided by 5)
for (i in seq((2*n.bins+1),(2*n.bins+2*n.bins),n.bins)) {
  
    tf.scale[ , i:( i + n.bins - 1 ) ] <- 100 / mean( sort( apply( tf[ , i:( i + n.bins - 1 )], 1, max ), decreasing=TRUE)[1:1000]/5 ) * tf[ , i:( i + n.bins - 1 ) ]

}

# Breaks and colour scheme
# PhyloP conservation
brk.pp          <- c( seq( min( tf.p[,2:ncol(tf.p)]),-1, length=5 ), seq(-0.999,1,length=40), seq(1.001,max(tf.p[,2:ncol(tf.p)]), length=5 ) )
col.pp          <- colorRampPalette(c("lightblue","lightblue","black","yellow","yellow"))(49)

# Accessibility (DNase hypersensitivity)
brk.dhs         <- c( seq( 0, 50, length=100 ), seq( 50.001, 200, length=60 ), seq( 200.001, max( tf.scale ), length=40 ) )
col.dhs         <- colorRampPalette(c("#EBEBEB","#4C77BB","#2D4770","#2D4770"))(199)

# H3K4me1 occupancy
brk.h3k         <- c( seq( 0, 50, length=100 ), seq( 50.001, 200, length=60 ), seq( 200.001, max( tf.scale ), length=40 ) )
col.h3k         <- colorRampPalette( c( "#EBEBEB", "#F5979A", "#7A4B4D", "#7A4B4D" ) )(199)

# RNAPII occupancy
brk.pol         <- c( seq( 0, 50, length=100 ), seq( 50.001, 200, length=60 ), seq( 200.001, max( tf.scale ), length=40 ) )
col.pol         <- colorRampPalette(c("#EBEBEB","#468453","#34633e","#34633e"))(199)

# Legends
z <- matrix(1:200,nrow=1)

pdf("legend_dhs_scaled.pdf", width=2, height=4)
par(mar=c(0,0,0,0) + 4)
image(1,brk.dhs,z,col=col.dhs,axes=F, ylim=c(0,200), xlab="", ylab="DNA occupancy (scaled)")
axis(2, at=c(0,50,100,150,200), labels=c("0","0.5","1","1.5","2"))
dev.off()

pdf("legend_h3k_scaled.pdf", width=2, height=4)
par(mar=c(0,0,0,0) + 4)
image(1,brk.h3k,z,col=col.h3k,axes=F, ylim=c(0,200), xlab="", ylab="DNA occupancy (scaled)")
axis(2, at=c(0,50,100,150,200), labels=c("0","0.5","1","1.5","2"))
dev.off()

pdf("legend_pol2_scaled.pdf", width=2, height=4)
par(mar=c(0,0,0,0) + 4)
image(1,brk.pol,z,col=col.pol,axes=F, ylim=c(0,200), xlab="", ylab="DNA occupancy (scaled)")
axis(2, at=c(0,50,100,150,200), labels=c("0","0.5","1","1.5","2"))
dev.off()


#
# Heatmap of chromatin accessibility, H3K4me1 and RNAPII occupancy and PhastCons DNA sequence conservation
# ----------------------------------------------------------------------------------------------------------------------------------------------

png(file="dhs_h3k4me1_pol2_conservation.png", width = ncol( tf ) - n.bins, height = nrow( tf ) )

    par( mar=c( 0, 0, 0, 0 ) )
    layout( matrix( 1:4, 1, 4, byrow = TRUE ) )
    image( t( tf.scale[ nrow( tf ):1, 1:( n.bins ) ] ), breaks = brk.dhs, col = col.dhs, axes = FALSE )
    image( t( tf.scale[ nrow( tf ):1, ( n.bins + 1 ):(n.bins*2)]), breaks = brk.h3k, col = col.h3k, axes = FALSE )
    image( t( tf.scale[ nrow( tf ):1, ( n.bins * 2 + 1 ):(n.bins*3)]), breaks = brk.pol, col=col.pol, axes = FALSE )
    image( t( tf[nrow(tf):1,(n.bins*4+1):( n.bins *5 ) ]), breaks=brk.pc, col = col.pc, axes = FALSE )

dev.off()


# Enrichment of DNA motifs in conserved and non-conserved CRMs
cons                    <- read.table("./xenTro91/MBT_cons_homer_motifs.txt", header=T)
ncon                    <- read.table("./xenTro91/MBT_NOTcons_homer_motifs.txt", header=T)

coln                    <- c("consensus","pval","log.pval","fdr","occ.motif.top10K","frac.pos.peak","occ.motif.bkg","frac.pos.bkg")

colnames(cons)          <- coln
colnames(ncon)          <- coln

cons <- within(cons, {
  frac.pos.peak         <- as.character(frac.pos.peak)
  frac.pos.bkg          <- as.character(frac.pos.bkg)
  rad.peak              <- sqrt(as.numeric(substr(frac.pos.peak, 1, nchar(frac.pos.peak)-1))/pi)
  rad.bkg               <- sqrt(as.numeric(substr(frac.pos.bkg, 1, nchar(frac.pos.bkg)-1))/pi)
  motif.enriched        <- as.numeric(substr(frac.pos.peak, 1, nchar(frac.pos.peak)-1)) / as.numeric(substr(frac.pos.bkg, 1, nchar(frac.pos.bkg)-1))}
)
colnames(cons)[2:11]    <- paste("cons",colnames(cons)[2:11],sep=".")

ncon <- within(ncon, {
  frac.pos.peak         <- as.character(frac.pos.peak)
  frac.pos.bkg          <- as.character(frac.pos.bkg)
  rad.peak              <- sqrt(as.numeric(substr(frac.pos.peak, 1, nchar(frac.pos.peak)-1))/pi)
  rad.bkg               <- sqrt(as.numeric(substr(frac.pos.bkg, 1, nchar(frac.pos.bkg)-1))/pi)
  motif.enriched        <- as.numeric(substr(frac.pos.peak, 1, nchar(frac.pos.peak)-1)) / as.numeric(substr(frac.pos.bkg, 1, nchar(frac.pos.bkg)-1))}
)
colnames(ncon)[2:11] <- paste("ncon",colnames(ncon)[2:11],sep=".")

list.dnase              <- list( cons, ncon )
all.dnase               <- plyr::join_all( list.dnase )
row.names(all.dnase)    <- row.names( cons )
nr                      <- nrow( all.dnase )

# bubble plot (occurence peak/background)
rad.bkg                 <- all.dnase[ , c(10,20) ]
rad.peak                <- all.dnase[ , c(11,21) ]

grid                    <- data.frame( matrix( 0, ncol = 4, nrow = 4 * nr ) )
grid$X1                 <- rep( substr( colnames(rad.peak), 1, nchar( colnames(rad.peak))-9), each=nr )
grid$X2                 <- row.names( rad.peak )
grid$X3                 <- unlist( list( rad.peak, rad.bkg ))
grid$X4                 <- rep( c("peak","background"), each = 2 * nr )

# factorise X1 (TFs) and X2 (motifs) to ensure order does not change when calling ggplot()
grid$X5                 <- factor( grid$X1, c("cons", "ncon") )
grid$X6                 <- factor( grid$X2, rev(row.names(all.dnase)) )


#
# Comparison of DNA motif enrichments between conserved and non-conserved accessible CRMs
# ------------------------------------------------------------------------------------------------------
pdf("MBT_consVsNonCons_accessible_crm_motifs.pdf", width=5, height=8, useDingbats=FALSE, family="Arial" )
ggplot(grid, aes(X5,X6,color=X4)) +
  geom_point(aes(size=X3^2*pi)) +
  scale_color_manual(values = c("black","#4C77BB")) +
  scale_size_continuous(range = range(c(0.1,10))) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

# -------------------------------------------------------------------------------------------------------



# BREAKS AND COLOR SCHEMES
# ----------------------------------------------------

required.pkg <- c( "RColorBrewer", "gplots" )

for ( pkg in required.pkg ) {
    if ( pkg %in% rownames(installed.packages()) == FALSE) {
        install.packages( pkg ) }
    if ( pkg %in% rownames(.packages()) == FALSE ) {
        library(pkg, character.only = TRUE) }
}

# Differential expression
brk.rna         <- c( seq(0,40,length=30), seq(40.001,66.666,length=30), seq(66.667,100,length=30), seq(100.001,150,length=30), seq(150.001,250,length=30), seq(250.001,2500,length=30) )
col.rna         <- colorRampPalette(c("#3953A4", "#AFBADA", "#EBEBEB", "#EBEBEB", "#FDAE61", "#FAA41A"))(179)

# Spatial expression
brk.spatio      <- c( seq(-10, -2.0001, length=20), seq( -2, 2, length=60 ), seq( 2.0001, 10, length=20 ) )
col.spatio.1    <- colorRampPalette( c( "yellow","yellow3","black","deepskyblue3","lightblue" ))(99)
col.spatio.2    <- colorRampPalette( c( "#739CD2","#648AB2","black","#C68990","#F6A3A9" ))(99)

# Maternal contribution
brk.mat         <- c( -.0001, .0999, .9999, 9.9999, 99.9999, 999.9999, 9999.9999 )
col.mat         <- brewer.pal( 7, "YlGnBu" )[2:7]

# Fold change ratios
brk.fc          <- seq( log(.249), log(4.001), length=50 )
col.fc          <- colorRampPalette( brewer.pal( 11, "Spectral" ))(49)[49:1]

# Z-score
brk.z           <- seq( -9, 9, length=100)
col.z           <- colorRampPalette( c(rep("#A2BAD9",1),rep("#4575B4",2),"black",rep("#FFFF00",2),rep("#FFFFB2",1)))(99)

# DNA motif P-value
brk.p.motif     <- seq( log10(5), log10(3000), length=100 )
brk.p.motif2    <- c( 0, seq( log(10), log(10000), length=99 ) )
col.p.motif     <- colorRampPalette( c( "white", "black" ) )(99)

# DNA motif enrichment
brk.motif.e      <- seq( log2(1/25), log2(25), length=100 )
col.motif.e      <- colorRampPalette(c("lightsalmon","white","royalblue3"))(99)

# Correlation factor
brk.corr         <- c( seq( -1, 1, length.out=20 ) )
col.corr         <- colorRampPalette( c( "red", "grey98", "royalblue4" ) )(19)
col.corr2        <- colorpanel( 79, "grey98", "royalblue4" )

# PhastCons conservation
brk.pc          <- seq( -0.001, 1.001, length=50 )
col.pc          <- colorRampPalette( c( "black", "yellow" ) )(49)

# -----------------------------------------------------

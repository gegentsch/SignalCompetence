# Read and plot motif DNA motif occurence and enrichment in peak and background regions

source("./R_utils/colour_schemes.R")

# FUNCTION: Read HOMER motif table
read.motif <- function( path, label ) {
    
    motif.tab           <- read.table( path, header=TRUE )
    
    colnames(motif.tab) <- c( "consensus", "pval", "log.pval", "fdr", "occ.motif.top2K", "frac.pos.peak", "occ.motif.bkg", "frac.pos.bkg" )
    
    motif.tab <- within(motif.tab, {
        frac.pos.peak   <- as.character( frac.pos.peak )
        frac.pos.bkg    <- as.character( frac.pos.bkg )
        rad.peak        <- sqrt( as.numeric( substr( frac.pos.peak, 1, nchar( frac.pos.peak ) -1 ) ) / pi )
        rad.bkg         <- sqrt( as.numeric( substr( frac.pos.bkg, 1, nchar( frac.pos.bkg ) -1 ) ) / pi )
        motif.enriched  <- as.numeric( substr( frac.pos.peak, 1, nchar( frac.pos.peak ) - 1 ) ) / as.numeric( substr( frac.pos.bkg, 1, nchar( frac.pos.bkg ) - 1 ) ) }
    )
    colnames( motif.tab )[ 2:11 ] <- paste( label, colnames( motif.tab )[ 2:11 ], sep="." )
    return( motif.tab )
}


# FUNCTION: Generate plots for motif occurence/enrichment
plot.motif <- function( dataframe, exp.name ) {
    
    # p-values (heatmap)
    m.pval              <- data.matrix( -dataframe[, grep( "log.pval", colnames(dataframe))] )
    m.pval[m.pval < 1]  <- 1
    
    pdf( file=paste0( exp.name, "_motifs_pval.pdf" ), width=10, height=10, useDingbats=FALSE, family="Arial" )
        par( mar=c(0,0,0,0) + 10 )
        image( t( log( m.pval[ nrow(m.pval):1, ]) ), breaks=brk.p.motif2, col=col.p.motif, axes=FALSE )
        axis( 2, seq( 0, 1, length=nrow( m.pval ) ), rev( row.names( dataframe ) ), las=1 )
        axis( 3, seq( 0, 1, length=ncol( m.pval ) ), colnames( m.pval ), las=2 )
    dev.off()
    
    # motif enrichment (heatmap)
    m.enrich            <- data.matrix( dataframe[, grep( "motif.enriched", colnames(dataframe))] )
    m.enrich[m.enrich == 0] <- NA
    
    pdf( file=paste0( exp.name, "_motifs_enrichment.pdf" ), width=10, height=8, useDingbats=FALSE, family="Arial" )
        heatmap.2( log2( m.enrich ), Rowv=FALSE, Colv=FALSE, dendrogram="none", col=col.motif.e, breaks=brk.motif.e, trace="none", cellnote=format( m.enrich, digits=1 ), key=TRUE,
        keysize=1.5, density.info="none", symkey=FALSE, notecol="black", margin=c(15,15) )
    dev.off()
    
    # bubble plot (peak vs background)
    rad.bkg             <- dataframe[, grep("rad.bkg", colnames(dataframe)) ]
    rad.peak            <- dataframe[, grep("rad.peak", colnames(dataframe)) ]
    
    n.cells             <- ( ncol( dataframe ) - 1 ) / 10 * nrow( dataframe )
    tf.stage            <- substr( colnames( rad.peak ), 1, nchar( colnames( rad.peak ) ) - 9 )
    
    grid                <- data.frame( matrix(0, ncol=4, nrow=(2*n.cells)) )
    grid$X1             <- rep( tf.stage, each=nrow( dataframe ) )
    grid$X2             <- row.names( rad.peak )
    grid$X3             <- unlist( list( rad.peak, rad.bkg ) )
    grid$X4             <- rep( c( "peak", "background" ), each=n.cells )
    
    # factorise X1 (TFs) and X2 (motifs) to ensure order does not change when calling ggplot()
    grid$X5             <- factor( grid$X1, tf.stage )
    grid$X6             <- factor( grid$X2, rev( row.names( dataframe ) ) )
    
    return(grid)
}

# ---------------------------------------------------------------------------------------------
# FUNCTION: calculate p-value of differential binding (in form of DNA binding matrix) using DESeq2
pValueMatrix <- function( win ) {
    
    row.names(win)      <- win$Gene
    win                 <- win[-1]
    
    n.profiles          <- as.numeric( unlist( strsplit( colnames( win )[ ncol(win) ], ".", fixed=TRUE ) )[2] ) + 1
    n.bins              <- ncol( win ) / n.profiles
    dm                  <- matrix( rep( 0, n.bins * nrow( win ) * n.profiles ), nrow=n.bins * nrow(win), ncol=n.profiles )
    row.names(dm)       <- paste0( rep( row.names( win ), each=n.bins ), ".", 1:n.bins )
    colnames(dm)        <- c( "c1", "c2", "oe1", "oe2" )
    
    for( i in 1:nrow( win ) ) {
        for( j in 1:n.bins ) {
            dm[ ( (i-1) * n.bins + j), ] <- unlist( win[ i, seq(j, ncol(win), n.bins) ] )
        }
    }
    
    col <- read.csv("./xenTro71/supplementary_files/colData.csv", header=T)
    
    dds <- DESeqDataSetFromMatrix(countData = round( dm ), colData = col, design = ~ condition)
    
    colData(dds)$condition <- factor(colData(dds)$condition, levels=c( "ctrl", "oe" ) )
    colData(dds)$replicate <- factor(colData(dds)$replicate, levels=c( "A", "B" ) )
    
    
    # tests for significance of coefficients in a Negative Binomial GLM, using previously calculated sizeFactors and dispersion estimates.
    # switch off log2 fold change shrinkage (betaPrior = FALSE) for this experiment
    # This function performs a default analysis through the steps:
    #   1. estimation of size factors: estimateSizeFactors
    #   2. estimation of dispersion: estimateDispersions
    #   3. Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
    
    dds                 <- DESeq( dds, betaPrior=FALSE, fitType="local" )
    res.dds             <- results( dds, contrast=c( "condition", "oe", "ctrl" ), cooksCutoff=FALSE, independentFiltering=FALSE )
    
    dmp                 <- matrix( rep( 1, nrow(win) * n.bins ), nrow=nrow( win ), ncol=n.bins )
    dmp                 <- matrix( res.dds$pvalue, nrow=nrow(win), ncol=n.bins, byrow=TRUE )
    dmp[ is.na( dmp ) ] <- 1
    dmp.log10           <- -log10( dmp )
    
    return( dmp.log10 )
}
# Installing (if required) and uploading several R packages
required.pkg <- c( "extrafont", "gplots" )

for ( pkg in required.pkg ) {
    if ( pkg %in% rownames(installed.packages()) == FALSE) {
        install.packages( pkg ) }
    if ( pkg %in% rownames(.packages()) == FALSE ) {
        library(pkg, character.only = TRUE) }
}

source( "./R_utils/colour_schemes.R" )


print.cor.pdf <- function( filename, data.matrix, display.cor=TRUE, display.size=10 ) {
    
    if(display.cor) {
        pdf( file=filename, width=40, height=40, useDingbats=FALSE, family="Arial" )
        heatmap.2( data.matrix, Rowv=FALSE, Colv=FALSE, dendrogram="none", trace="none", symm=TRUE, scale="none", breaks=brk.corr, col=col.corr, key=TRUE, symkey=TRUE, density.info="none", cellnote=round( data.matrix, 1 ), notecol="black", notecex=display.size, cexRow=display.size/2, cexCol=display.size/2, margin=c(40,40), lwid=c(1,10), lhei=c(1,10) )
        dev.off()
    } else {
        pdf( file=filename, width=40, height=40, useDingbats=FALSE, family="Arial" )
        heatmap.2( data.matrix, Rowv=FALSE, Colv=FALSE, dendrogram="none", trace="none", symm=TRUE, scale="none", breaks=brk.corr, col=col.corr, key=TRUE, symkey=TRUE, density.info="none", cexRow=display.size/2, cexCol=display.size/2, margin=c(40,40), lwid=c(1,10), lhei=c(1,10) )
        dev.off()
    }
    
}
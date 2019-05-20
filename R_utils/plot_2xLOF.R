# FUNCTION: Plot two LOFs and color dots according to another feature
# --------------------------------------------------------------------

plotDifSp <- function( fc.x, fc.y, select, cl, x.lim.min=0, x.lim.max=1000, y.lim.min=0, y.lim.max=1000, pdf.name, log=FALSE, label=FALSE, height=5.5 ) {
    
    fc.x_ <- if(log) log( fc.x / 100 +.01 ) else fc.x
    fc.y_ <- if(log) log( fc.y / 100 +.01 ) else fc.y
    
    s <- ( fc.x <= x.lim.max & fc.y <= y.lim.max )
    
    x.lim.min_ <- ifelse( log, log(x.lim.min/100 +0.01), x.lim.min )
    x.lim.max_ <- ifelse( log, log(x.lim.max/100), x.lim.max )
    y.lim.min_ <- ifelse( log, log(y.lim.min/100 +0.01), y.lim.min )
    y.lim.max_ <- ifelse( log, log(y.lim.max/100), y.lim.max )
    
    thsld <- if(log) c(log(2/3),log(1),log(1.5)) else c(t66,100,150)
    
    labels <- if(label) paste0( 1:length(select),":",sapply(rn[select],f4) ) else 1:length(select)
    
    at.1 <- if(log) c( x.lim.min_,log(.1),log(2/3),log(1),log(1.5),x.lim.max_ ) else c( x.lim.min_,t66,100,150,x.lim.max_)
    at.2 <- if(log) c( y.lim.min_,log(.1),log(2/3),log(1),log(1.5),y.lim.max_ ) else c( y.lim.min_,t66,100,150,y.lim.max_)
    
    labels.1 <- if(log) c( exp(1)^x.lim.min_,0.1,"",1,"",exp(1)^x.lim.max_ ) else c( x.lim.min_,"",100,"",x.lim.max_)
    labels.2 <- if(log) c( exp(1)^y.lim.min_,0.1,"",1,"",exp(1)^y.lim.max_ ) else c( y.lim.min_,"",100,"",y.lim.max_)
    
    pdf ( pdf.name, width=5, height=height, useDingbats=F, family="Arial" )
    plot( fc.x_[s], fc.y_[s], xlim=c( x.lim.min_, x.lim.max_ ), ylim=c( y.lim.min_, y.lim.max_ ), ann=F, axes=F, cex=1.4, pch=20, lwd=.25, col = alpha(cl[s], .7) )
    abline( h=thsld[1],v=thsld[1],lty="22",col="grey50",lwd=.5 ); abline(h=thsld[2],v=thsld[2],col="grey50",lwd=.5); abline(h=thsld[3],v=thsld[3],lty="22",col="grey50",lwd=.5)
    points( fc.x_[select], fc.y_[select], cex=1.6, pch=21, lwd=.25, col="black", bg=alpha(cl[select], 1 ) )
    text( fc.x_[select], fc.y_[select], labels=labels, pos=4, offset=.3, cex=.9 )
    title( main=pdf.name, xlab=strsplit(pdf.name, "_", fixed=T)[[1]][1], ylab=strsplit(pdf.name, "_", fixed=T)[[1]][3], cex.main=.8 )
    axis( 1, at=at.1, labels=labels.1, lwd=.5 )
    axis( 2, at=at.2, labels=labels.2, lwd=.5 )
    dev.off()
    
}
# FUNCTION: plot mean +/- standard deviation of read density (logarithmic scale) across a 2-kb genomic region

plotMetaLogNorm <- function (df, colour, min, max) {
    
    lstd 	<- sapply(log10(df+0.1),sd)
    dist 	<- seq(-950,950,by=25)
    lmean 	<- colMeans(log10(df+0.1))
    plot(dist, lmean, xlim=c(-1000,1000),ylim=c(log10(min),log10(max)),type="n", axes=F, ann=F)
    polygon(c(dist,rev(dist)),c(lmean+lstd,rev(lmean-lstd)),col=alpha( colour, 0.5 ), border = F)
    lines(dist, lmean, col=colour)
    axis( 1, at=c(-1000,-500,0,500,1000), labels=c("-1","","0","","+1") )
    axis( 2, at=c(log10(0.1),log10(1),log10(10),log10(100),log10(1000)), labels=c("0.1","1","10","100","1000") )

}




# FUNCTION: plot mean +/- standard error of DNA motif density across a 2-kb genomic region

plotMetaMotif <- function (df, colour, max) {

    sem 	<- sapply(df,function(x) sd(x)/sqrt(length(x)))
    dist 	<- seq(-950,950,by=25)
    mean	<- colMeans(df)
    plot(dist, mean, xlim=c(-1000,1000),ylim=c(0,max),type="n", axes=F, ann=F)
    polygon(c(dist,rev(dist)),c(mean+sem,rev(mean-sem)),col=alpha( colour, 0.5 ), border = F)
    lines(dist, mean, col=colour)
    axis( 1, at=c(-1000,-500,0,500,1000), labels=c("-1","","0","","+1") )
    axis( 2, at=c(0,max), labels=c("0",paste(max)) )

}
rm( list=ls() )
load( 'Rdata/compiled_results.Rdata' )
load( 'Rdata/setup.Rdata' )
apply( !is.na(allpvals[1,,,,]), 1:3, sum )

pdf( 'fig.pdf', width=29/2, height=12/2 )
layout( cbind( matrix(1:8,2,4), 9 ), widths=c( 6, 6, 6, 6, 5 ) )

par(mar=c(5,5,4,4)) 
for( zzz in 1:4 )
	for( sig2.i in 1:2 )
{
	if( zzz == 2 ){
		ys	<- allpvals['Hom'													,,'1e4_K1',,]
	} else {
		ys	<- allpvals[c('Null','','Hom','Het')[zzz]	,,'N1e4'  ,,]
	}

	lims  <- c( 0, 6 )
	plot( lims, lims, type='n', xlab='Expected -log10(p)', ylab='Observed -log10(p)', main=c( 'Null Effects', 'Hom Effects (K=1)', 'Hom Effects', 'Het Effects' )[zzz], cex.lab=1.5 )
	abline( a=0, b=1 ) 

	addline <- function(y,col){
		y <- sort( y[!is.na(y)] )
		n <- length(y)
		x <- sort( -log10( 1:n/(n+1) ) )
		lines(	x, y, col=col, lty=1  )
		points(	x, y, col=col, pch=16 )
	}

	for( meth in rev(methods) )
		addline(ys[meth,sig2.i,],cols[meth])
}

par( mar=c( 2, 0, 2, 0 ) )
plot.new()
legend( 'left', fill=cols, leg=nicemethods, cex=1.8, border=F, bty='n' )

dev.off() 

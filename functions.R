ldprune	<- function( y1, maf=.05, r2=.1, bfile='parse_gala_data/data/GALAII_all_freeze_041112' ){
	if( length(y1) == 0 ) return(y1)
	write.table( names(y1), file='Rout/tmp.extract', quote=F, row.names=F, col.names=F )
	system( paste0( '~/misc/plink --threads 1 --bfile ', bfile, ' --extract Rout/tmp.extract --indep-pairwise 1000 1 ', r2, ' --out Rout/tmp --maf ', maf, ' --chr 1-22 --geno .1' ) )
	rsids	<- read.table( 'Rout/tmp.prune.in' )[,1]
	y1[rsids]
}

add_y_qq	<- function(y,col,add_labs=TRUE){
	y	<- sort( y[!is.na(y)] )
	n	<- length( y )
	x	<- sort( -log10( 1:n/(n+1) ) )
	points( x, y, pch=16, col=col )

	fdrs	<- p.adjust( 10^-y, 'BH' )

	if( length(y) > 1e4 ){
		cexs	<- 2.5
	} else {
		cexs	<- 2.5 + ( fdrs < .2 ) * 4
		if( add_labs ){
			labs	<- names(y)
			labs[ fdrs >= .2 ]	<- ''
			delta	<- rep(0,length(x))
			if( sum( fdrs < .2 ) == 3 ) delta[1:3 + length(x) - 3] <- c( -.21, -.02, .04 )
			text( x+delta, y+1.2, labs, pch=16, col=col, cex=3.8, srt=90 )
		}
	}
	points( x, y, pch=16, col=col, cex=cexs )
	return(names(y)[( fdrs < .2 )])
}

my_qq	<- function(y,xlab='Expected -log10(p)', ylab='Observed -log10(p)',maftol=.05,misstol=.10,cand_rs,legloc='bottomright',miny=8,add_labs=TRUE, prune_lists=FALSE, r2=.1,...){

	y[mafs < maftol	]	<- NA
	y[miss > misstol]	<- NA

	y	<- sort( y[!is.na(y)] )

	lam_gc<- lam_gc_fxn( y )

	lims	<- range( c( y , miny ) )
	plot( lims, lims, type='n', xlab='', ylab='', main='', axes=F, ... )
	box()
	axis( 1, cex.axis=3.8, padj=.9 )
	axis( 2, cex.axis=3.8, padj=-.5 )
	mtext( side=1, line=8, xlab, cex=3.0 )
	mtext( side=2, line=8, ylab, cex=3.0 )
	abline( a=0, b=1 )
	abline( h=-log10( 5e-8 ), col=1 )

	topsnps <- list( add_y_qq( y, col=1, add_labs=add_labs ) ) 
	base	<- NULL
	if( !missing( cand_rs ) )
		for( i in 1:length(cand_rs) )
	{
		y1	<- y[cand_rs[[i]]]
		y1	<- y1[!is.na(y1)]
		y1	<- ldprune( y1, r2=r2 )
		if( prune_lists & length(y1) > 0 & length(base) > 0 ){
			manual_prune	<- function( y, base ){
				write.table( c( y, base ), file='Rout/tmp.r2list.extract', quote=F, row.names=F, col.names=F )
				system( paste0( '~/misc/plink --threads 1 --bfile ', bfile, ' --extract Rout/tmp.r2list.extract --r2 square --out Rout/tmp.r2s' ) )
				x	<- read.table( file='Rout/tmp.r2s.ld', head=F )[1,]#, quote=F, row.names=F, col.names=F )
				( max(x[-1]) < .2 )
			}
			snpsub	<- which(sapply( names(y1), function(y) manual_prune( y, base=base ) ))
			y1	<- y1[snpsub]
		}
		base		<- c( base, cand_rs[[i]] )
		topsnps <- c( topsnps, list( add_y_qq( y1, col=i+1, add_labs=add_labs ) ) )
	}

	legend( legloc, bty='n', leg=paste0( 'Lambda_GC=', round( lam_gc, 3 ) ), cex=3.8 )
	list( topsnps=topsnps, lam_gc=lam_gc )
}

ssapply	<- function(x,f)
	sapply( x, function(x.i)
{
	out	<- NA
	try( out	<- f(x.i), silent=T )
	out
})

lam_gc_fxn	<- function( log10p ){
	chi2s		<- qchisq(1-10^-log10p,df=1)
	median(chi2s,na.rm=T)/qchisq(0.5,1)
}

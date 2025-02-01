rm( list=ls() )
load( 'Rdata/compiled_results.Rdata' )
load( 'Rdata/setup.Rdata' )
Ns <- unique(Ns)
x			<- log10( Ns )

pvals	<- apply( allpvals, 1:4, function(x) mean( x > -log10(.05), na.rm=T )	)

all_pvals	<- array( NA,
	dim=c(					4													,Me			, length(x)	, 2 ),
	dimnames=list( c('null','K=1','hom','het'),methods, x					, sig2hets )
)

all_pvals[c('null','hom','het'),,,]	<- pvals[			,,c( 'N1e3'  ,'N3e3'  ,'N1e4'  ),]
#all_pvals[c('null','hom','het'),,,]	<- pvals[			,,c( '1e3_noasc'  ,'3e3_noasc'  ,'1e4_noasc'  ),]
all_pvals['K=1',,1:3,]							<- pvals['Hom',,c( '1e3_K1','3e3_K1','1e4_K1'),]
#types		<- c( '1e3_noasc'	,'3e3_noasc','1e4_noasc')
if( 'oracle' %in% methods )
all_pvals['K=1','oracle',,]	<- NA


pdf( 'figure.pdf', width=23, height=5.9 )
layout( cbind( rbind(
	c( 1, 3, 4, 5, 2, 6 ),
	c( 7, 14+c(1:3,5,4) ),
	8+c( 1, 3, 4, 5, 2, 6 )
	),8), widths=c( 1.7, rep( 6, 3 ), 2.4, 6, 4.8 ), heights=c( 1.0, 5.3, 1.8 ) ) 

##### top
par( mar=rep(0,4) )
plot.new()
plot.new()
for( zzz in 1:4 ){
plot.new()
text( .5, .45, cex=4, lab=c( 'Null Effects', 'Hom Effects (K=1)', 'Hom Effects', 'Het Effects' )[zzz] )
}
plot.new() ##### left

##### right
par( mar=c( 10, 0, 10, 0 ) )
plot.new()
legend( 'center', fill=cols, leg=nicemethods, cex=3, border=F, bty='n' )#, horiz=T )

##### bottom
par( mar=rep(0,4) )
plot.new()
plot.new()
par( mar=c( 1.0, 2, 1.0, 2.0 )/2 )
for( zzz in 1:4 ){
plot(range(x)+c(-.1,.1),0:1,type='n',axes=F,ylab='',xlab='')
text( mean(range(x)), .25, cex=3.3, lab='Sample Size (N)' )
axis( 1, cex.axis=2.4,  at=x, line=-9.9, lab=F )
axis( 1, cex.axis=2.4,  at=x, line=- 9.7, lab=c( '1,000', '3,000', '10,000' ), padj=.9, tick=F )
#axis( 1, cex.axis=2.4,  at=x, line=- 9.7, lab=c( '1,000', '3,000', '10,000', '30,000' ), padj=.9, tick=F )
}
##### meat
par( mar=c( 1.0, 1, 1.0, 1.0 )/2 )
for( zzz in 1:4 )
{
	allys	<- all_pvals[zzz,,,]
	if( zzz != 4 ){
	ylim	<- c( -2.2, 0 )
	ats	<- c(.01	,.05	,.2		,1)
	ats	<- log10(ats)
	labs<- c('.01','.05','.2'	,'1' )
	allys	<- log10(allys)
	} else {
	ylim	<- c( 0, 1 )
	ats		<- c(0,.25,.5	,.75,1)
	labs	<- ats
	}

	plot( range(x)+c(-.1,.1), ylim, type='n', axes=F, ylab='', main='', xlab='' )
	box()
	legend( 'topleft', bty='n', cex=4, leg=letters[zzz], adj=c(2.5,-.3) )

	#abline( a=log10(.05)	, b=0, col='lightgrey', lty=3, lwd=8 )
	abline( a=log10(.01)	, b=0, col='lightgrey', lty=3, lwd=8 )
	if( zzz %in% c(1,4) ){
		axis( 2, cex.axis=2.7, padj=-.3	, at=ats			, lab=labs )
		mtext( side=2, line=6.2, cex=2.3 , ifelse( zzz == 1, 'False Positive Rate', 'True Positive Rate' ) )
	}

	for( meth in rev(methods) )
		for( sig2.i in 1:2 )
	{
		lines(	x, allys[meth,,sig2.i], col=cols[meth], lty=c(1,2)[sig2.i], lwd=ifelse( meth %in% c( 'oracle', 'mfmr3' ), 5, 3 ) )
		points(	x, allys[meth,,sig2.i], col=cols[meth], pch=18-sig2.i, cex=4 )
	}
}
plot.new()

dev.off()

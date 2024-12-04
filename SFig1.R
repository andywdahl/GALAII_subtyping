rm( list=ls() )
load('Rdata/setup.Rdata')

Kmax		<- 6
n.folds	<- 100
type		<- 'imp'

ps	<- array( NA, dim=c( Kmax, n.folds ), dimnames=list( 1:Kmax, 1:n.folds ) )
for( K in 1:Kmax )
{
	savefile	<-	paste0( 'Rdata/findK_ll_'	, K, '_', type, '.Rdata' )
	load( savefile )
	ps	[K,]	<- out[,1]
	rm( out )
}
print( t( apply( is.na(ps), 1, which ) ) ) ### 2/100 folds fail, but this is consistent across all choices for K and driven by rare binary phenotypes and covariates giving singularities in test set--not actually mfmr failing

for( k in Kmax:1 )
	ps[k,]	<- ps[k,]-ps[1,]

pdf( 'figs/SFig1.pdf', width=6.8, height=6.5 )
par( mar=c(5.2,5.2,1,1) )

plot( c(1,Kmax), range(ps,na.rm=T)*.95, type='n', xlab='Number of Latent Subtypes (K)', ylab='Cross-Validated Log Likelihood Ratio', main='', axes=F, cex.lab=1.7 )
axis(2,cex=1.1)
axis(1,cex=1.1,at=1:Kmax)

for( j in 1:n.folds )
	lines(	1:K, ps[,j], col=1, lwd=.7, lty=3 )

ysd	<- apply( ps, 1, sd		, na.rm=T ) / sqrt(n.folds)
ps	<- apply( ps, 1, mean	, na.rm=T )
arrows(2:K, (ps-ysd)[-1], 2:K, (ps+ysd)[-1], length=0.05, angle=90, code=3, col=2, lwd=2.4)

points(	1:K, ps, col=1, pch=16, cex=2 )
lines(	1:K, ps, col=1, lwd=7, lty=1 )
k0	<- which.max(ps)
points(	k0, ps[k0], col=3, pch=16, cex=2.5 )

legend( 'bottomleft', bty='n', lty=c(1,1,3,NA), pch=c(16,NA,NA,16), col=c(1,2,1,3), leg=c( 'Mean over folds', 'Std. Err. over folds', '100 individual folds', 'Optimal K' ), lwd=c(4,2,1.5,NA), pt.cex=c(1.5,NA,NA,1.8), cex=1.3 )

dev.off()

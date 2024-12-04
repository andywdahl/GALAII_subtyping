rm( list=ls() )
library(GxEMM)
load( 'Rdata/setup.Rdata' )
sink( 'Rout/Fig4D.Rout' )

phens	<- c( "asthma_onset_frombirth", "Pre.FEV1.FVC.perc.pred", "deltafev", "tIGE" ) 
ylabs	<- c( 'Age-of-Onset Heritability', 'FEV1/FVC Heritability', '', '' ) 
type	<- 'ctr'
qn		<- FALSE ### was already quantnormed, so TRUE automatically

h2s		<- array( NA, dim=c( length(phens), 6				), dimnames=list( phens, c( 'hom', 'hom1', 'iid'	, paste0( 'h2', 1:3 ) ) ) )
h2se	<- array( NA, dim=c( length(phens), 6				), dimnames=list( phens, c( 'hom', 'hom1', 'iid'	, paste0( 'h2', 1:3 ) ) ) )
ps		<- array( NA, dim=c( length(phens), 5				), dimnames=list( phens, c( 'hom', 'iid', 'iid_lr', 'free', 'free_lr' ) ) )
ps_e	<- array( NA, dim=c( length(phens), 5				), dimnames=list( phens, c( 'hom', 'iid', 'iid_lr', 'free', 'free_lr' ) ) )
ps0		<- array( NA, dim=c( length(phens), 5, 1001	), dimnames=list( phens, c( 'hom', 'iid', 'iid_lr', 'free', 'free_lr' ), 1:1001 ) )
h2Covmats	<- lapply( phens, function(x) NA )
names(h2Covmats)	<- phens

for( phen in phens ){
	load( paste0( 'Rdata/gxemm/', phen, '_', type, '_qn=', qn, '.Rdata' ) ) 
	
	ps		[phen,]	<- out$pvals[c( 'hom', 'iid', 'iid_lr', 'free', 'free_lr' )]
	h2s		[phen,]	<- c( out$fits$hom$h2[1], out$fits$iid$h2[1:2], out$fits$free$h2 )
	h2se	[phen,]	<- sqrt(c( out$fits$hom$h2Covmat[1], diag(out$fits$iid$h2Covmat), diag(out$fits$free$h2Covmat) ))
	h2Covmats[[phen]] <- out$fits$free$h2Covmat

	for( k in 1:1000 ){
		load( paste0( 'Rdata/gxemm/', phen, '_rand', k, '_', type, '_qn=', qn, '.Rdata' ) )
		ps0	[phen,,k]	<- out$pvals[c( 'hom', 'iid', 'iid_lr', 'free', 'free_lr' )]
		rm(out)
	}
	print(apply( ps0, 1:2, function(x) sum(is.na(x)) ) )
	ps0[phen,,1001]	<- ps[phen,] 
	ps0[is.na(ps0)] <- 0 ### conservatively, when GxEMM crashes on permuted data, assume worst possible result (p=0)
	ps_e[phen,]	<- sapply( colnames(ps), function(x) mean( ps[phen,x] >= ps0[phen,x,], na.rm=T ) )
} 
print(round(ps	[,c('iid','iid_lr','free_lr')],12))
print(round(ps_e[,c('iid','iid_lr','free_lr')],3))

h2s_all	<- h2s
h2se_all<- h2se
rm( h2s, h2se ) 

load( 'Rdata/mfmrx_3_imp.Rdata' )
perm<- 3:1
for( phen in c( "Pre.FEV1.FVC.perc.pred", "asthma_onset_frombirth" ) ){
	print( phen )
	cat( 'H2:   ', round( colMeans( out$pmat ) %*% t(h2s_all[,paste0('h2',1:3)]), 3 ), '\n' )
	cat( 'H2se: ', round( sqrt(colMeans( out$pmat ) %*% (h2Covmats[[which(phens==phen)]] %*% colMeans( out$pmat ))), 3 ), '\n' )
	cat( 'h2ses:', round( diag(h2Covmats[[which(phens==phen)]]), 3 ), '\n' )

	h2s	<- h2s_all	[phen,c('hom',paste0( 'h2',perm ))]
	h2se<- h2se_all	[phen,c('hom',paste0( 'h2',perm ))]

	pdf( paste0( 'figs/Fig4D_', phen, '.pdf' ), width=3.2, height=5.5 )
	par( mar=c(.5,4.5,.5,1) )

	xs		<- c( -.27, -.08+.29*seq(0,1,length=3) )
	ylim	<- c( -.2, 1.22 )
	xlim	<- range(xs)			+ c(-1,1)*.05 

	plot( range(xs), ylim, axes=F, type='n', xlab='', ylab='', main='', xlim=xlim ) #, 
	axis(2, cex.axis=1.2, at=0:5/5, lab=paste0( 0:5/5*100, '%') )
	mtext(2,text=ylabs[which(phens==phen)], cex=1.5, line=3.1 )

	yy	<- 1.2
	lines( xs[c(2,4)]+c(-1,1)*.05, c(yy,yy), lwd=2 )
	text( xs[3], yy+.04, paste0( 'p=', format( ps[phen,'free_lr'], dig=2 ) ), cex=1.3 ) 

	abline( h=0, col=2, lty=3, lwd=2 )
	points(xs, h2s, pch=16, cex=4																, col=rep(cols,length(phens)) )
	arrows(xs, h2s-h2se, xs, h2s+h2se, length=0.09, angle=90, code=3, col=rep(cols,length(phens)), lwd=3.5 )

	dev.off() 
}
sink()

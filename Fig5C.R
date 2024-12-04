rm( list=ls() )
library(rgwas)
library(parallel)
load('Rdata/setup.Rdata')

taus	<- c( '0.000000000001','0.0000000001', '0.00000001', '0.000001', '0.00001', '0.0001', '0.001', '0.01', '0.1', '1.0' )
x			<- 1:length(taus)

phens	<- c( 'disease_ASTHMA_DIAGNOSED', 'disease_ALLERGY_ECZEMA_DIAGNOSED', 'blood_EOSINOPHIL_COUNT', 'cov_SMOKING_STATUS', 'mental_NEUROTICISM' )
xlabs	<- c( 'Asthma'									, 'Atopy'														, 'Eosinophil'						, 'Smoking'						, 'Neuroticism'				 )
names(xlabs)	<- phens

K			<- 3
nperm	<- 1e3

allbetas	<- array( NA, dim=c(K+1	,length(taus),length(phens)      ), dimnames=list(c('hom',1:K)					,taus,phens          ) )
allses		<- array( NA, dim=c(K+1	,length(taus),length(phens)      ), dimnames=list(c('hom',1:K)					,taus,phens          ) )
allps			<- array( NA, dim=c(3		,length(taus),length(phens),nperm), dimnames=list(c('hom','het','glob')	,taus,phens,0:(nperm-1)) ) 
for( phen in phens )
	for( perm.i in 0:(nperm-1) )
{
	load( paste0( 'Rdata/prs/mlr_perm_', perm.i, '_', K, '_ctr_', phen, '.Rdata' ) )
	for( tau in taus ){
		if( phen == 'mental_NEUROTICISM' )
			if( betas[3,tau] < betas[1,tau] | betas[3,tau] < betas[2,tau]  )
				ps	[2,tau]	<- 0
		if( phen == 'disease_ALLERGY_ECZEMA_DIAGNOSED' )
			if( betas[1,tau] < betas[2,tau] | betas[1,tau] < betas[3,tau]  )
				ps	[2,tau]	<- 0
	}
	if( perm.i==0 ){
		allbetas[-1,,phen]	<- betas
		allses	[-1,,phen]	<- ses
		rm( ses )
	}
	allps[,,phen,perm.i+1]	<- 10^-ps 
	rm( ps, betas )

	if( perm.i==0 ){
		load( paste0( 'Rdata/prs/lr_', K, '_ctr_', phen, '.Rdata' ) )
		allbetas[1, ,phen]	<- betalr
		allses	[1, ,phen]	<- selr
		rm( betalr, selr ) 
	}
} 

betas	<- array( NA, dim=c(K+1	,length(phens)), dimnames=list(c('hom',1:K)					,phens) )
ses		<- array( NA, dim=c(K+1	,length(phens)), dimnames=list(c('hom',1:K)					,phens) )
ps		<- array( NA, dim=c(3		,length(phens)), dimnames=list(c('hom','het','glob'),phens) ) 
for( phen in phens ){
	ps		[,phen]	<- sapply( 1:3, function(k){
		minps	<- apply( allps[k,,phen,], 2, min )
		mean( minps[1] >= minps )
	})
	j		<- which.min( allps['glob',,phen,1] )
	betas	[,phen]	<- allbetas	[,j,phen]
	ses		[,phen]	<- allses		[,j,phen] 
}

perm	<- 3:1
subnames	<- subnames[perm]
betas	<- betas[c(1,1+perm),]
ses		<- ses	[c(1,1+perm),]
round( betas, 3 )


pdf( 'figs/Fig5C.pdf', width=8.5, height=5.5 )
par( mar=c(9,5.5,1,1) )

xs		<- rep( 1:length(phens), each=K+1 )+rep( c( -.27, -.08+.29*seq(0,1,length=K) ), length(phens) )*1.1 
ys		<- c(betas)
yse		<- c(ses) 
ylim	<- range( c( 0, ys-yse, ys+yse ), na.rm=T )
xlim	<- range(xs)			+ c(-1,1)*.05 


plot( range(xs), ylim, axes=F, type='n', xlab='', ylab='', main='', xlim=xlim ) #, ylab=ylabs[p]
axis(2, cex.axis=1.2)
mtext(2,text='Asthma Risk (Log OR)', cex=1.7, line=3.8 )

for( i in phens ){
	x.i	<- which( phens==i )
	axis(1,at=x.i,lab=xlabs[i], cex.axis=1.30, tick=F, line=.3 )
	for( ii in 1:3 )
	axis(1,at=x.i,lab=round( ps[ii,i]	, 3 ), cex.axis=1.35,lwd=0,line=1.0+ii*1.8,tick=F)
}
abline( h=0, col=2, lty=3, lwd=2 )
points(xs, ys, pch=16, cex=4																, col=rep(cols,length(phens)) )
arrows(xs, ys-yse, xs, ys+yse, length=0.09, angle=90, code=3, col=rep(cols,length(phens)), lwd=3.5 )

mtext(1,outer=T,line=-7.6,adj=.06,text='PRS'								, cex=1.45)
mtext(1,outer=T,line=-5.4,adj=.06,text=expression( p[hom] )	, cex=1.45)
mtext(1,outer=T,line=-3.3,adj=.06,text=expression( p[het] )	, cex=1.45)
mtext(1,outer=T,line=-1.2,adj=.06,text=expression( p[glob])	, cex=1.45)

dev.off()

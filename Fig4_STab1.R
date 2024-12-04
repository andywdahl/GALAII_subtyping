rm( list=ls() )
load('Rdata/setup.Rdata') 
K				<- 3

### tested covariates and 2ndary phenotypes
locphens<- c( "asthma_onset_frombirth", 'deltafev', "tIGE", '#wheeze', "hosp_er_12mo", "rhinitis", "hayfever", "eczema", "rash", "sinusitis")
loccov	<- c("AFR","NAM",'Male','MX', 'BMI', smokenames) 

##################### get p-values from droptests
load( 'Rdata/droptest_3_imp.Rdata' ) 
drop_ps	<- out$pvals['Het',loccov,locphens]
drop_ps[ is.na( drop_ps ) ]	<- 1
rm(out)

load( 'Rdata/droptest_mlr_3_imp.Rdata' )
drop_ps	<- cbind( sapply( out[loccov], function(x) 10^-x$log10p_het ), drop_ps )
colnames(drop_ps)[1]	<- 'mlr'
rm(out)

fdrmat	<- round(matrix( p.adjust( c(drop_ps), 'BH' ), nrow(drop_ps), ncol(drop_ps) ),2)
dimnames(fdrmat)	<- dimnames(drop_ps)

##################### get effect sizes from post-subtyping fit
load( 'Rdata/mlrx_3_imp.Rdata' )
betas.mlr	<- betas		[,paste0( 'X', loccov )]
ses.mlr		<- betas.se	[,paste0( 'X', loccov )]
rm( betas, betas.se )

load( 'Rdata/mlrx_3_ctr.Rdata' ) ### use ctr-adjusted estimates for adm/ethn effects
betas.mlr	[,'XMX']	<- betas		[,'XMX']
ses.mlr		[,'XMX']	<- betas.se	[,'XMX'] 
betas.mlr	[,'XNAM']	<- betas		[,'XNAM']
ses.mlr		[,'XNAM']	<- betas.se	[,'XNAM']

rm( betas, betas.se )

load( 'Rdata/mfmrxx_3_imp.Rdata' )
out		<- c( out_Yb, out_Yq )[locphens]
betas	<- lapply( out, function(outi){
	mm	<- outi[paste0( 'Z', 1:K, ':G', rep(loccov,each=K) ),1]
	mm	<- matrix( mm, nrow=K, ncol=length(loccov) )
	dimnames(mm)	<- list( 1:K, loccov )
	mm
})

ses	<- lapply( out, function(outi){
	mm	<- outi[paste0( 'Z', 1:K, ':G', rep(loccov,each=K) ),2]
	mm	<- matrix( mm, nrow=K, ncol=length(loccov) )
	dimnames(mm)	<- list( 1:K, loccov )
	mm
})
rm( out )

load( 'Rdata/mfmrxx_3_ctr.Rdata' ) ### use ctr-adjusted estimates for adm/ethn effects
out		<- c( out_Yb, out_Yq )[locphens]
for( p in 1:length(betas) ){
	betas[[p]][,'MX']	<- out[[p]][paste0( 'Z', 1:K, ':G', rep('MX' ,each=K) ),1]
	ses  [[p]][,'MX']	<- out[[p]][paste0( 'Z', 1:K, ':G', rep('MX' ,each=K) ),2]
	betas[[p]][,'NAM']<- out[[p]][paste0( 'Z', 1:K, ':G', rep('NAM',each=K) ),1]
	ses  [[p]][,'NAM']<- out[[p]][paste0( 'Z', 1:K, ':G', rep('NAM',each=K) ),2]
}
rm( out )

### glue together MLR with case-only interaction tests for plotting
betas	<- c( list(betas.mlr)	, betas )
ses		<- c( list(ses.mlr)		, ses )
names(betas)[1]	<- 'mlr'
names(ses)[1]		<- 'mlr'

mnames	<- names(mains)
mains		<- c( 'Asthma', mains )
names(mains)<- c('mlr',mnames)


### misc setup for plots 
mains0	<- mains
mains[ mains=="Exacerbation last year" ]	<- "Exacerbation\nlast year" 
cols	<- cols[-1] 

perm			<- 3:1
betas			<- lapply( betas, function(x) x[perm,] )
ses				<- lapply( ses	, function(x) x[perm,] )
subnames	<- subnames[perm]


pdf( 'figs/Fig4.pdf', width=7.7*1.2, height=8*1.2 )
layout( cbind( 7:8, matrix( 1:6, 2, 3, byrow=T ) ), widths=c(.41,2.2,1,1) ) 
par( mar=c(5,4,8,2) )

tau	<- .05 ### display all nominally-significant effects
locphens2	<- colnames(drop_ps)[ which( apply( drop_ps, 2, function(x) any( x < tau ) ) ) ] 
for( p in locphens2 ){

	is	<- which( drop_ps[,p] < tau )
	ys	<- c( betas[[p]][,is] )
	yse	<- c( ses	[[p]][,is] )
	pvs	<- drop_ps[,p] 
	fdr	<- fdrmat[,p] 

	xs	<- rep( 1:length(is), each=K )+rep( .17*seq(-1,1,length=K), length(is) )

	plot( range(xs), range( c( 0, ys-yse, ys+yse ) ), axes=F, type='n', ylab='', xlab='', main=mains[p], xlim=range(xs)+c(-.2,.2), cex.lab=1.6, cex.main=1.6 )
	axis(2,cex.axis=1.5)
	abline( h=0, col=2, lty=1, lwd=1.5 )

	fonts	<- sapply( fdr[is], function(x) ifelse( x > .2, 1, 2 ) )
	for( i in 1:length(is) ){
		pv	<- format( pvs[is[i]], digits=2, scientific=T )
		axis(1,at=i,lab=mains[loccov[is[i]]], cex.axis=1.4,line=1.1	,tick=F,font=fonts[i])
		axis(1,at=i,lab=pv									, cex.axis=1.4,line=2.8	,tick=F,font=fonts[i])
	}

	points(xs, ys, pch=16, cex=3																, col=rep(cols,length(is)) )
	arrows(xs, ys-yse, xs, ys+yse, length=0.05, angle=90, code=3, col=rep(cols,length(is)), lwd=2.2)

	if( p %in% locphens2[c(1,4)] )
		mtext(2,text='Subtype-Specific\nEffect Sizes', cex=1.4, line=4.1 )
}

mtext(1,outer=T,padj=-33.9,adj=.10,text=expression( p[het] ), cex=1.2)
mtext(1,outer=T,padj=-01.0,adj=.10,text=expression( p[het] ), cex=1.2)

dev.off()

### write output to STab 1 (betas from post-hoc fit, pvalues from droptest)
mytab	<- matrix( NA, nrow=0, ncol=6 )
colnames(mytab)	<- c( 'pheno', 'covariate', 'subtype', 'beta', 'se', 'het-pvalue (from droptest)' )
subnames	<- c( 'T2-', 'T2+_Atopy-', 'T2+_Atopy+' )
mains[mains == "Albuterol\nResponse" ]	<- "Albuterol Response"
mains[mains == "Exacerbation\nLast Year" ]	<- "Exacerbation Last Year"

for( p in c( 'mlr', locphens ) ){
	mytab1	<- cbind( mains0[p], rep( loccov, each=K ), rep( subnames, length(loccov) ), c(betas[[p]]), c(ses[[p]]), rep( drop_ps[,p], each=K ) )
	mytab	<- rbind( mytab, mytab1 )
}
write.table( mytab, quote=F, file='figs/STab1.csv', sep=',', row.names=F )

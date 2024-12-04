rm( list=ls() )
source( 'functions.R' )
load( 'Rdata/setup.Rdata' )
sink( 'Rout/qqplots.Rout' )

mains['deltafev']			<- "Albuterol Response"
mains['hosp_er_12mo']	<- "Exacerbation Last Year"
mains		<- c( 'Asthma', mains )
names(mains)[1]	<- 'mlr'

locphens	<- c( 
	'deltafev', 'asthma_onset_frombirth', 'tIGE', "Pre.FEV1.FVC.perc.pred", 'Pre.FEV1.perc.pred',
	'#wheeze', "hosp_er_12mo", "symptoms_nocold", "hayfever", "sinusitis", "rhinitis", "symptoms_12mo" 
)
for( K in 2:4 ){

load( 'Rdata/ast_snps.Rdata' )
ast_snps	<- intersect( ast_snps, rsids )

topsnps	<- vector("list",length(locphens)+1)
names(topsnps)	<- c( 'mlr', locphens )

load( paste0( 'Rdata/gwas_ctr.Rdata' ) )
names(out)	<- rsids
phomscc	<- ssapply( out, function(i) -log10(i[4]) )
subcc		<- names(phomscc)[ which( phomscc > -log10( 1e-3 ) ) ]
rm(out) 

tiff( paste0( 'figs/qqplots_', K, '.tiff' ), width=40*100, height=40*100 )
par( mfrow=c(4,4), mar=c(14,14,4,4) )
for( phen in names(topsnps) ) try({

	if( phen == 'mlr' ){
		load( paste0( 'Rdata/gwas_mlr_'	, K, '_ctr.Rdata' ) )
		names(out)	<- rsids
		phoms	<- phomscc 
		phets	<- ssapply( out, function(x) x$log10p_het )
		pglobs<- ssapply( out, function(x) x$log10p_glob )
		cand_rs	<- list( ast_snps, subcc ) 
	} else {
		load( paste0( 'Rdata/gwas_caseonly_'	, phen, '_'	, K, '_ctr.Rdata' ) )
		phoms	<- ssapply( out, function(x) -log10(x$homp0) )
		sub		<- names(phoms)[ which( phoms   > -log10( 1e-3 ) ) ]
		phets	<- ssapply( out, function(x) -log10(x$hetp) )
		pglobs<- ssapply( out, function(x) -log10(x$globp) )
		cand_rs	<- list( ast_snps, subcc, sub )
	}
	rm(out)

	out	<- my_qq( y=phets	, ylab='log10(p) for Het Association', cand_rs=cand_rs )

	cat( K, phen, 'lambda_gc=', round(out$lam_gc,2), '\n' )

	topsnps[[which(names(topsnps)==phen)]]	<- lapply( out$topsnps, function(x){
		if( length(x) == 0 ) return(NA)
		cbind( phomscc[x], phoms[x], phets[x], pglobs[x] )
	})
	rm( out, phoms, phets, pglobs, cand_rs ) 
	legend( 'topleft', bty='n', cex=5.8, leg=mains[phen] )

	if( phen == 'Pre.FEV1.perc.pred' ){ plot.new(); plot.new() }
})

par( mar=rep(0,4) )
plot.new()
legend( 'center', bty='n', pch=16, col=1:4, leg=c( 'Genome-wide SNPs', 'Known Asthma GWAS SNPs', 'Nominal Hom Asthma Effect', 'Nominal Hom Secondary Trait Effect' ), cex=6.6 )
dev.off()

print( 'topsnps' )
print( topsnps[1:6] )

}

sink()

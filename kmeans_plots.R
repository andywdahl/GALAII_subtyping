rm( list=ls() )
source( 'functions.R' )
load( 'Rdata/setup.Rdata' )
sink( 'Rout/qqplots.Rout' )

types			<- c( 'ctr'		, 'kmeans', 'kmeans_noadjust' )
legs			<- c( 'RGWAS'	, 'k-means', 'k-means ignoring pop strat' )
names(legs)	<- types

tiff( 'figs/kmeans_mlr.tiff', width=40*100, height=30*100 )
par( mfcol=c(3,4), mar=c(14,14,4,4) )

for( K in 2:4 )
	for( type in types )
{

	if( K == 4 & type == 'kmeans' ){
		plot.new()
		next
	}

	load( 'Rdata/ast_snps.Rdata' )
	ast_snps	<- intersect( ast_snps, rsids )

	load( paste0( 'Rdata/gwas_ctr.Rdata' ) )
	names(out)	<- rsids
	phomscc	<- ssapply( out, function(i) -log10(i[4]) )
	subcc		<- names(phomscc)[ which( phomscc > -log10( 1e-3 ) ) ]
	rm(out) 

	load( paste0( 'Rdata/gwas_mlr_'	, K, '_', type, '.Rdata' ) )
	names(out)	<- rsids
	phoms	<- phomscc 
	phets	<- ssapply( out, function(x) x$log10p_het )
	pglobs<- ssapply( out, function(x) x$log10p_glob )
	cand_rs	<- list( ast_snps, subcc ) 

	my_qq( y=phets	, ylab='log10(p) for Subtype-specific Effect', cand_rs=cand_rs, adjust_gc=FALSE, legloc='topleft', miny=16, add_labs=FALSE )
	rm( out, phoms, phets, pglobs, cand_rs )

	legend( 'bottomright', bty='n', cex=5.8, leg=paste0( legs[type], ', k=', K ) )
}

par( mar=rep(0,4) )
plot.new()
plot.new()
legend( 'center', bty='n', pch=16, col=1:3, leg=c( 'Genome-wide SNPs', 'Known Asthma GWAS SNPs', 'Nominal Hom Asthma Effect' ), cex=6.6 )
plot.new()

dev.off()

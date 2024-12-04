rm( list=ls() )
library(rgwas)
library(parallel)
load('Rdata/setup.Rdata')

pvs	<- c( '0.000000000001','0.0000000001', '0.00000001', '0.000001', '0.00001', '0.0001', '0.001', '0.01', '0.1', '1.0' )
phens	<- c( 'disease_ASTHMA_DIAGNOSED', 'disease_ALLERGY_ECZEMA_DIAGNOSED', 'blood_EOSINOPHIL_COUNT', 'mental_NEUROTICISM' )

K			<- 3
nperm	<- 1e3
type	<- 'ctr'

for( perm.i in 0:nperm )
	for( phen in phens )
{
	loadfile	<- paste0( 'Rdata/mfmrx_'											, K, '_', 'imp', '.Rdata' )
	savefile	<- paste0( 'Rdata/prs/mlr_perm_', perm.i, '_'	, K, '_', type, '_', phen, '.Rdata' )
	sinkfile 	<- paste0( 'Rout/prs_mlr_perm_'	, perm.i, '_'	, K, '_', type, '_', phen, '.Rout' )
	if( file.exists( savefile ) | file.exists( sinkfile ) | ! file.exists( loadfile ) )	next
	sink(	sinkfile )

	betas	<- array( NA, dim=c(K,length(pvs)), dimnames=list(1:K,pvs) )
	ses		<- array( NA, dim=c(K,length(pvs)), dimnames=list(1:K,pvs) ) ### only done for perm=0
	ps		<- array( NA, dim=c(3,length(pvs)), dimnames=list(1:3,pvs) )
	for( pv in pvs ){

		### load in cc, G, X
		load(paste0( 'parse_gala_data/parsed_data/final_ctr.Rdata' ))
		X[ X[,'centerOK']==1, ]	<- NA
		X	<- X[ ,-which(colnames(X)=='centerOK')]

		G	<- G[,intersect( colnames(G), Gnames_cc )]
		X	<- cbind( 1, apply( cbind(G,X), 2, rgwas:::scale01 ) )

		### load in Z
		load( loadfile )
		Z		<- rep( '0', length(cc) )
		Z[which( cc == 1 )]	<- as.character(apply( out$pmat, 1, function(x) ifelse( max(x) > 0, which.max(x), NA ) ))
		names(Z)	<- rownames(X)
		rm( G, Yb, Yq, cc, out )

		prs		<- read.table( paste0( 'parse_gala_data/parsed_data/', phen, '.', pv, '.profile' ), head=T )
		sub		<- intersect( names(Z), prs[,1] )

		prs		<- prs[ match( sub, prs[,1] ),]
		Z			<- Z[sub]
		X			<- X[sub,]

		prs	<- scale( prs[,'SCORE'] )
		if( perm.i > 0 ){
			set.seed( perm.i )
			prs	<- sample(prs,replace=F)
		}

		out	<- run_mlr( X, prs, Z, scaleXg=FALSE )

		betas	[,pv]	<- out$betas
		ses		[,pv]	<- out$betas.se
		ps		[,pv]	<- c( out$log10p_hom, out$log10p_het, out$log10p_glob )
	}
	save( betas, ses, ps, file=savefile )
	sink()
	rm( out, prs, Z, X )
}

### fit logistic regression for plotting the additive effects
for( phen in phens ){

	loadfile	<- paste0( 'Rdata/mfmrx_'	, K, '_', 'imp', '.Rdata' )
	savefile	<- paste0( 'Rdata/prs/lr_', K, '_', type, '_', phen, '.Rdata' )
	sinkfile 	<- paste0( 'Rout/prs_lr_'	, K, '_', type, '_', phen, '.Rout' )
	if( file.exists( savefile ) | file.exists( sinkfile ) | ! file.exists( loadfile ) )	next
	sink(	sinkfile )

	betalr	<- array( NA, dim=c(length(pvs)), dimnames=list(pvs) )
	selr		<- array( NA, dim=c(length(pvs)), dimnames=list(pvs) )
	for( pv in sample(pvs) ){
		### load in cc, G, X
		load(paste0( 'parse_gala_data/parsed_data/final_ctr.Rdata' )) 
		G	<- G[,intersect( colnames(G), Gnames_cc )]
		X	<- cbind( 1, apply( cbind(G,X), 2, rgwas:::scale01 ) )
		rm( G, Yb, Yq )

		prs		<- read.table( paste0( 'parse_gala_data/parsed_data/', phen, '.', pv, '.profile' ), head=T )
		sub		<- intersect( rownames(X), prs[,1] )

		prs		<- prs[ match( sub, prs[,1]  ),]
		X			<- X[sub,]
		cc		<- cc[sub]

		prs	<- scale( prs[,'SCORE'] )

		out	<- summary( glm( cc ~ X + prs, family=binomial(link='logit') ) )
		betalr[pv]	<- out$coef['prs',1]
		selr	[pv]	<- out$coef['prs',2] 
		rm( out, prs, cc, X )
	}
	save( betalr, selr, file=savefile )
	sink()
}

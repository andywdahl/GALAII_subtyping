rm( list=ls() )
library(rgwas)
library(parallel)
library(BEDMatrix)
load('Rdata/setup.Rdata')

K			<- as.numeric( commandArgs(TRUE)[[1]] )
mc.co	<- as.numeric( commandArgs(TRUE)[[2]] )
locphens	<- c( 'deltafev', 'asthma_onset_frombirth', 'tIGE', "Pre.FEV1.FVC.perc.pred", 'Pre.FEV1.perc.pred', '#wheeze', "hosp_er_12mo", "symptoms_nocold", "hayfever", "sinusitis", "rhinitis", "symptoms_12mo" )
for( phen in locphens ){

	loadfile	<- paste0( 'Rdata/mfmrx_'		, K, '_imp.Rdata' )
	savefile	<- paste0( 'Rdata/gwas_caseonly_'	, phen, '_'	, K, '_ctr.Rdata' )
	sinkfile 	<- paste0( 'Rout/gwas_caseonly_'	, phen, '_'	, K, '_ctr.Rout' )
	if( file.exists( savefile ) | file.exists( sinkfile ) | ! file.exists( loadfile ) )	next
	print( sinkfile )
	sink(	sinkfile ) 

	### load in cc, G, X
	load(paste0( 'parse_gala_data/parsed_data/final_ctr.Rdata' ))
	X[ X[,'centerOK']==1, ]	<- NA
	X	<- X[ ,-which(colnames(X)=='centerOK')]
	G	<- G[,intersect( colnames(G), Gnames_cc )]
	X	<- cbind( 1, apply( cbind(G,X), 2, rgwas:::scale01 ) ) 

	if( phen %in% Ybnames ){
		y		<- Yb[which( cc == 1 ),phen]
	} else {
		y		<- Yq[which( cc == 1 ),phen]
	}
	X			<- X [which( cc == 1 ),,drop=F]
	rm( Yq, Yb, G )

	### load in Z
	load( loadfile )
	Z	<- out$pmat
	rownames(Z)	<- rownames(X)
	rm( out )

	Gsnp		<- BEDMatrix( 'parse_gala_data/data/GALAII_all_freeze_041112' )
	snprows	<- as.character(sapply( rownames(Gsnp), function(y) strsplit( y, '_' )[[1]][1] ))
	sub			<- intersect(snprows,rownames(X))
	rownames( Gsnp )	<- snprows
	y				<- y		[sub]
	X				<- X		[sub,,drop=F]
	Z				<- Z		[sub,,drop=F]

	out	<- mclapply( 1:ncol(Gsnp), mc.cores=mc.co, function(s){
		print(s)
		out.i	<- NA
		try({
			g			<- Gsnp[sub,s]
			subi	<- which(!is.na(g))
			y			<- y [subi]
			X			<- X [subi,]
			Z			<- Z [subi,]
			g			<- g [subi]

			if( phen %in% Ybnames ){
				nonlm0<- glm( y ~ X						, family='binomial' )
				homlm0<- glm( y ~ X + g				, family='binomial' )
				nonlm	<- glm( y ~ X + Z				, family='binomial' )
				homlm	<- glm( y ~ X + Z + g 	, family='binomial' )
				hetlm	<- glm( y ~ X + Z + Z:g , family='binomial' )
			} else if( phen %in% c(Ynames, 'tIGE_raw' ) ){
				nonlm0<- lm( y ~ X     )
				homlm0<- lm( y ~ X + g )
				nonlm	<- lm( y ~ X + Z )
				homlm	<- lm( y ~ X + Z + g )
				hetlm	<- lm( y ~ X + Z  + Z:g )
			} else {
				stop()
			}
			homp0	<- anova( homlm0, nonlm0, test='Chisq')$Pr[2]
			homp	<- anova( homlm	, nonlm , test='Chisq')$Pr[2]
			hetp	<- anova( hetlm	, homlm , test='Chisq')$Pr[2]
			globp	<- anova( hetlm	, nonlm , test='Chisq')$Pr[2]
			out.i	<- list( homp0=homp0, homp=homp, hetp=hetp, globp=globp, coef=summary( hetlm )$coef )
		})
		out.i
	})
	names(out)	<- rsids

	save( out, file=savefile )
	print(warnings())
	print('Done')
	sink()
	rm( out, cc, Z, X )
}

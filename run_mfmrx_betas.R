rm( list=ls() )
library(rgwas)
library(parallel)
load('Rdata/setup.Rdata')

K				<- as.numeric( commandArgs(TRUE)[[1]] )
for( type in c( 'imp', 'ctr' ) ){ ### ctr is important confounder for interpreting the MX ethnicity effect

	loadfile	<- paste0( 'Rdata/mfmrx_'	, K, '_imp.Rdata' )
	savefile	<- paste0( 'Rdata/mfmrxx_', K, '_', type, '.Rdata' )
	if( file.exists( savefile ) )	next

	load(paste0( 'parse_gala_data/parsed_data/final_', type, '.Rdata' ))
	cases	<- which( cc == 1 )
	Yb	<-Yb[cases,,drop=F]
	Yq	<-Yq[cases,,drop=F]
	G		<-G [cases,,drop=F]

	G	<- G[,intersect( colnames(G), Gnames_cc )]
	G	<- scale(G)
	X	<- cbind( 1, apply( cbind(G,X), 2, rgwas:::scale01 ) )

	### load in Z
	load( loadfile )
	Z	<- out$pmat
	rm( out )

	out_Yq	<- lapply( colnames(Yq), function(j){
		yq	<- Yq[,j]
		summary( lm( yq ~ -1 + Z + X:Z ) )$coef
	})
	names( out_Yq )	<- colnames(Yq)

	out_Yb	<- lapply( colnames(Yb), function(j){
		yb	<- Yb[,j]
		summary( glm( yb ~ -1 + Z + X:Z, family='binomial' ) )$coef
	})
	names( out_Yb )	<- colnames(Yb)

	save( out_Yb, out_Yq, file=savefile )
	rm( Yb, Yq, out_Yb, out_Yq, G, X )

}

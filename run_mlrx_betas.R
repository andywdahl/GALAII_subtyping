rm( list=ls() )
library(rgwas)
library(parallel)
load('Rdata/setup.Rdata')

K				<- as.numeric( commandArgs(TRUE)[[1]] ) 
for( type in c( 'imp', 'ctr' ) ){ ### ctr is important confounder for interpreting the MX ethnicity effect

	loadfile	<- paste0( 'Rdata/mfmrx_'	, K, '_imp.Rdata' )
	savefile	<- paste0( 'Rdata/mlrx_'	, K, '_', type, '.Rdata' )
	if( file.exists( savefile ) | ! file.exists( loadfile ) )	next

	### load in cc, G, X
	load(paste0( 'parse_gala_data/parsed_data/final_', type, '.Rdata' ))
	G	<- G[,intersect( colnames(G), Gnames_cc )]
	X		<- cbind( 1, apply( cbind(G,X), 2, rgwas:::scale01 ) )
	rm( G, Yb, Yq )

	### load in Z
	load( loadfile )
	Z		<- rep( '0', length(cc) )
	Z[which( cc == 1 )]	<- as.character(apply( out$pmat, 1, function(x) ifelse( max(x) > 0, which.max(x), NA ) ))
	rm( out )

	suppressMessages( out   <- summary(nnet::multinom( Z ~ -1+X, Hess=T, maxit=1e3, MaxNWtsA=1e4 )) )
	betas     <- out$coefficients
	betas.se  <- out$standard.errors

	save( betas, betas.se, file=savefile )
	print(warnings())
	print('Done')
	rm( out, cc, Z, X )

}

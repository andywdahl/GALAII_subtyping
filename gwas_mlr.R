rm( list=ls() )
library(rgwas)
library(parallel)
library(BEDMatrix)
load('Rdata/setup.Rdata')
types			<- c( 'ctr', 'kmeans', 'kmeans_noadjust' ) 
loadtypes	<- c( 'imp', 'kmeans', 'kmeans_noadjust' ) 
names(loadtypes)	<- types

K			<- as.numeric( commandArgs(TRUE)[[1]] )
mc.co	<- as.numeric( commandArgs(TRUE)[[2]] )
for( type in types ){

	loadfile	<- paste0( 'Rdata/mfmrx_'		, K, '_', loadtypes[type], '.Rdata' )
	savefile	<- paste0( 'Rdata/gwas_mlr_', K, '_', type, '.Rdata' )
	sinkfile 	<- paste0( 'Rout/gwas_mlr_'	, K, '_', type, '.Rout' )
	if( file.exists( savefile ) | file.exists( sinkfile ) | ! file.exists( loadfile ) )	next
	print( sinkfile )
	sink(	sinkfile )

	### load in cc, G, X
	if( type %in% c( 'ctr', 'kmeans' ) ){
		load(paste0( 'parse_gala_data/parsed_data/final_ctr.Rdata' ))
		X[ X[,'centerOK']==1, ]	<- NA
		X	<- X[ ,-which(colnames(X)=='centerOK')]
	}  else {
		load(paste0( 'parse_gala_data/parsed_data/final_imp.Rdata' ))
	}
	G	<- G[,intersect( colnames(G), Gnames_cc )]
	X	<- cbind( 1, apply( cbind(G,X), 2, rgwas:::scale01 ) ) 
	if( type %in% c( 'kmeans_noadjust' ) )
		X	<- X[,-which( colnames(X) %in% c( 'MX', admnames ) )]

	### load in Z
	load( loadfile )
	Z		<- rep( '0', length(cc) )
	Z[which( cc == 1 )]	<- as.character(apply( out$pmat, 1, function(x) ifelse( max(x) > 0, which.max(x), NA ) ))
	names(Z)	<- rownames(X)
	rm( G, Yb, Yq, cc, out )

	Gsnp		<- BEDMatrix( 'parse_gala_data/data/GALAII_all_freeze_041112' )
	snprows	<- as.character(sapply( rownames(Gsnp), function(y) strsplit( y, '_' )[[1]][1] ))
	sub			<- intersect(snprows,rownames(X))
	rownames(Gsnp)	<- snprows
	X				<- X	[sub,]
	Z				<- Z	[sub]
	print( dim( Gsnp ) )

	out	<- mclapply( 1:ncol(Gsnp), mc.cores=mc.co, function(s){
		print(s)
		out.i	<- NA
		tryCatch({ out.i	<- run_mlr( X, Gsnp[sub,s], Z, scaleXg=FALSE ) },error=function(e) print(e))
		out.i
	})
	names(out)	<- colnames(Gsnp)

	save( out, file=savefile )
	print(warnings())
	print('Done')
	sink()
	rm( out, cc, Z, X )
}

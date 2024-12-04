rm( list=ls() )
library(mixtools)
load('Rdata/setup.Rdata')

K				<- as.numeric( commandArgs(TRUE)[[1]] ) 
for( type in c( 'kmeans', 'kmeans_noadjust' ) ){

	savefile	<-	paste0( 'Rdata/mfmrx_', K, '_', type, '.Rdata' )
	sinkfile 	<-	paste0( 'Rout/mfmrx_'	, K, '_', type, '.Rout' )
	if( file.exists( savefile ) | file.exists( sinkfile ) )	stop()
	print( sinkfile )
	sink(	sinkfile )

	load( 'parse_gala_data/parsed_data/final_imp.Rdata' )
	G	<- G[,intersect( colnames(G), Gnames_cc )]

	if( type == 'kmeans_noadjust' )
		G	<- G[,setdiff( colnames(G), c( 'MX', admnames ) )]

	cases	<- which( cc == 1 )
	Yb	<-Yb[cases,,drop=F]
	Yq	<-Yq[cases,,drop=F]
	G		<-G [cases,,drop=F]
	X		<-X [cases,,drop=F]

	Y	<- scale( cbind( Yb, Yq, G, X ) )
	runtime	<- system.time({
		 out  <- kmeans( Y, centers=K )#, iter.max=1e3, nstart=10, trace=T )
		 out$pmat  <- model.matrix(~as.factor(out$cluster)-1)
	})[3]

	Z	<- out$pmat

	out_Yq	<- lapply( colnames(Yq), function(nam){
		yq	<- Yq[,nam]
		summary( lm( yq ~ -1 + Z + G:Z ) )$coef
	})
	names( out_Yq )	<- colnames(Yq)

	out_Yb	<- lapply( colnames(Yb), function(nam){
		yb	<- Yb[,nam]
		summary( glm( yb ~ -1 + Z + G:Z, family='binomial' ) )$coef
	})
	names( out_Yb )	<- colnames(Yb)

	save( out, out_Yb, out_Yq, file=savefile )
		sink()
	rm( Yb, Yq, out, out_Yb, out_Yq, G, X )
}

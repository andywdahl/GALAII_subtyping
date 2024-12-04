rm( list=ls() )
library(rgwas)
library(parallel)
load('Rdata/setup.Rdata')

K				<- as.numeric( commandArgs(TRUE)[[1]] )
mc.cores<- as.numeric( commandArgs(TRUE)[[2]] )

savefile	<-	paste0( 'Rdata/findK_ll_'	, K, '_imp.Rdata' )
sinkfile 	<-	paste0( 'Rout/findK_ll_'	, K, '_imp.Rout' )
if( file.exists( savefile ) | file.exists( sinkfile ) )	next
print( sinkfile )
sink(	sinkfile )

load( 'parse_gala_data/parsed_data/final_imp.Rdata' )
G	<- G[,intersect( colnames(G), Gnames_cc )]

cases	<- which( cc == 1 )
Yb	<- Yb[cases,,drop=F]
Yq	<- Yq[cases,,drop=F]
G		<- G [cases,,drop=F]
X		<- X [cases,,drop=F]

runtime	<- system.time({
	out	<- score_K( Yb=Yb, Yq=Yq, G=cbind( 1, G ), X=X, K, n.folds=100, mc.cores=mc.cores, trace=T )
})[3]

save( out, runtime, file=savefile )
print(warnings())
print('Done')
sink()

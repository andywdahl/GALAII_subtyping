rm( list=ls() )
library(rgwas)
library(parallel)
load('Rdata/setup.Rdata')
type		<- 'imp'

K				<- as.numeric( commandArgs(TRUE)[[1]] )
mc.cores<- as.numeric( commandArgs(TRUE)[[2]] )

savefile	<-	paste0( 'Rdata/mfmrx_', K, '_', type, '.Rdata' )
sinkfile 	<-	paste0( 'Rout/mfmrx_'	, K, '_', type, '.Rout' )
if( file.exists( savefile ) | file.exists( sinkfile ) )	next
print( sinkfile )
sink(	sinkfile )

load(paste0( 'parse_gala_data/parsed_data/final_', type, '.Rdata' ))
G	<- G[,intersect( colnames(G), Gnames_cc )]

cases	<- which( cc == 1 )
Yb	<-Yb[cases,,drop=F]
Yq	<-Yq[cases,,drop=F]
G		<-G [cases,,drop=F]
X		<-X [cases,,drop=F]

runtime	<- system.time({
	out	<- mfmr( Yb=Yb, Yq=Yq, G=cbind(1,G), X=X, K=K, trace=T, nrun=20, mc.cores=mc.cores, init_sd=0.001 )
})[3]

save( out, runtime, mc.cores, file=savefile )
print(warnings())
print('Done')
sink()
rm( Yb, Yq, out, G, X, runtime )


### output subtype responsibilities
if( K == 3 ){
	load(paste0( 'parse_gala_data/parsed_data/final_', type, '.Rdata' ))
	load(paste0( 'Rdata/mfmrx_', K, '_', type, '.Rdata' ) )
	pmat	<- out$pmat
	colnames(pmat)	<- subnames
	rownames(pmat)	<- names(cc)[which( cc == 1 )] 
	write.table( pmat, file='figs/subtypes.csv', row.names=T, col.names=subnames, sep=',' )
}

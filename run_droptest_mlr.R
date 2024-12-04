rm( list=ls() )
library(rgwas)
library(parallel)
load('Rdata/setup.Rdata')
type		<- 'imp'

K				<- as.numeric( commandArgs(TRUE)[[1]] )
mc.cores<- as.numeric( commandArgs(TRUE)[[2]] )

savefile	<- paste0( 'Rdata/droptest_mlr_', K, '_', type, '.Rdata' )
sinkfile 	<- paste0( 'Rout/droptest_mlr_'	, K, '_', type, '.Rout' )

if( file.exists( savefile ) | file.exists( sinkfile ) )	next
print( sinkfile )
sink(	sinkfile )

load(paste0( 'parse_gala_data/parsed_data/final_', type, '.Rdata' ))
G	<- G[,intersect( colnames(G), Gnames_cc )]

runtime	<- system.time({
	out	<- droptest_mlr( cc, Yb, Yq, G=cbind(1,G), X=X, test_inds=1+1:ncol(G), trace=T, K=K, nrun=20, mc.cores=mc.cores, init_sd=0.001 )
})[3]

save( out, runtime, mc.cores, file=savefile )
print(warnings())
print('Done')
sink()
rm( Yb, Yq, out, G, X, runtime )

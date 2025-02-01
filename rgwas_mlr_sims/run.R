rm( list=ls() )
library(rgwas)
source( 'simfxn.R' )
source( 'main.R' )
load( 'Rdata/setup.Rdata' )

it  <- as.numeric( commandArgs(TRUE)[[1]] )
set.seed( it )

for( type in types )
	for( method in sample(methods) )
		for( sig2.i in sample(length(sig2hets)) )
try({
	savefile  <- paste0( 'Rdata/', sig2.i, '_', method, '_', type, '_', it, '.Rdata' )
	sinkfile  <- paste0( 'Rout/' , sig2.i, '_', method, '_', type, '_', it, '.Rout'  )
	if( file.exists(savefile) | file.exists(sinkfile) ) next
	print( sinkfile )
	sink( sinkfile )
	try({

		simdat  <- sim_fxn( Ns[type], K=Ks[type],asc=ascs[type],
		P=P, w=c(.3,.7), sig2E=.1, S=12, S_hom=4, S_het=4, prev=.2,
		sig2homs[sig2.i], sig2hets[sig2.i],  seed=it )

		runtime <- system.time( out <- main( Y=simdat$Y, G=simdat$G, method=method, K=2 ) )[3]

		save( runtime, out, file=savefile )
		rm( simdat, out )
	})
	sink()
})
warnings() 

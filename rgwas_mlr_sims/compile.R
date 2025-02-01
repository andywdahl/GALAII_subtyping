rm( list=ls() )
load( 'Rdata/setup.Rdata' ) 

snptypes	<- c( 'Null', 'Hom', 'Het' )
allpvals	<- array( NA, dim=c(3,Me,length(types),2,4*maxit), dimnames=list(snptypes,methods,types,sig2hets,1:(4*maxit)) )
for( type in types )
	for( method in sample(methods) )
		for( sig2.i in sample(length(sig2hets)) )
			for( it in 1:maxit )
{
	loadfile	<- paste0( 'Rdata/', sig2.i, '_', method, '_', type, '_', it, '.Rdata' )
	if( !file.exists( loadfile ) ) next
	load(	loadfile )
	for( xx in snptypes ){
		if( xx == 'Null' ){
			snps	<- 8+1:4
		} else if( xx == 'Hom' ){
			snps	<- 1:4
		} else {
			snps	<- 4+1:4
		}
		if( method != 'rgwas+' ){
			allpvals[xx,method,type,sig2.i,(it-1)*4+1:4]	<- unlist(out['log10p_het',snps]) #c('log10p_hom','log10p_het','log10p_glob') 
		} else {
			allpvals[xx,method,type,sig2.i,(it-1)*4+1:4]	<- sapply( out[-1], function(x) x$log10p_het )[snps]
		}
	}
	rm(out)
}

save( allpvals, file='Rdata/compiled_results.Rdata' )
source( 'fig.R' )

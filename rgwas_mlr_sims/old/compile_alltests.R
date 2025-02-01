rm( list=ls() )
load( 'Rdata/setup.Rdata' ) 

pvals	<- array( NA, dim=c(S,P,Me, maxit ), dimnames=list(1:S	, 1:P	, methods	, 1:maxit ))
for( sig2.i in 1:2 )
	for( method in methods )
		for( it in 1:maxit )
tryCatch({
#if( mean( is.na( pvals[-1,,method,it] ) ) < .5 ) next 
	loadfile	<- paste0( 'Rdata/', sig2.i, '_', method, '_', type, '_', it, '.Rdata' )
	if( !file.exists( loadfile ) ) next
	load(	loadfile )
	#print( loadfile )
	print( out )
	pvals[,,method,it]	<- out$pvals[1:S + ifelse( method %in% c( 'rgwas', 'rgwas+' ), 1, 0 ),]
	rm(out)
},error=function(e){ print(loadfile); print(e); print( method ) })
save( pvals, file='Rdata/compiled_results_raw.Rdata' )

snptypes	<- c( 'Null', 'Hom', 'Het' )
allp	<- array( NA, dim=c(3,3,Me,length(types),length(sig2hets)), dimnames=list(snptypes,c('cc','quant','bin'),methods,types,sig2hets) )
for( type in sample(types) )
	for( sig2.i in sample(2) )
{ 
	for( xx in snptypes )
		for( pp in c('cc','quant','bin') )
	{
		if( pp == 'cc' ){
			phens	<- 1
		} else if( pp == 'bin' ){
			phens	<- 2:P_bin
		} else {
			phens	<- (P_bin+1):P
		}
		if( xx == 'Null' ){
			snps	<- 8+1:4
		} else if( xx == 'Hom' ){
			snps	<- 1:4
		} else {
			snps	<- 4+1:4
		}
		allp	[dir,xx,pp,,type,sig2.i]	<- apply( pvals[snps,phens,,,drop=F], 3	, function(x) mean( x<.01, na.rm=T )	)
	}
	#nruns	[dir,,type,sig2.i]	<- apply( pvals[1,1,,], 1		, function(x) sum( ! is.na(x) )				)
	rm( pvals )
}
save( pvals, file='Rdata/compiled_results.Rdata' )
source( 'fig1.R' )

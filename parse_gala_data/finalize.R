rm( list=ls() )
library(phenix)
sink( 'parsed_data/make_Y.Rout' )
load( '../Rdata/setup.Rdata' )

### raw contains un-imputed data, used in some plots
### drop is row-wise complete data, not used as the sample size reduces >80%
### ctr additionally keeps assessment center, an additional covariate in genetic tests
for( type in c( 'imp', 'drop', 'raw', 'ctr' ) ){
	print( type )
	load( 'parsed_data/covar.Rdata'	)
	print( colnames(covar) )

	if( type != 'raw' ){
		Yq	<- covar[,setdiff( Ynames, '#wheeze' )]
		Yq	<- quantnorm( Yq )
		colnames(Yq)	<- setdiff( Ynames, '#wheeze' )
		covar[,setdiff( Ynames, '#wheeze' )]	<- Yq
	} else {
		covar[,'tIGE']	<- log10( covar[,'tIGE'] )
	}

	if( type %in% c( 'drop' ) ){
		print( table( cc ) )
		misses	<- colMeans( is.na(covar[which(cc==1),]) )
		covar		<- covar[,misses<.5]
		sub			<- which( rowSums( is.na(covar) ) == 0 )
		sub			<- union( sub, which( cc ==  0 ) )
		covar		<- covar	[sub,]
		cc			<- cc			[sub ]
		print( table( cc ) )
		print( dim( covar ) )
	} else if( type == 'raw' ){
	} else {
		covnames<- colnames(covar)
		covar		<- MVN_impute(covar,trace=T)$Y
		colnames(covar)	<- covnames
	}

	Yq	<- covar[,setdiff( Ynames, '#wheeze' )]
	Yb	<- covar[,Ybnames]

	col0	<- colnames(Yq)
	Yq		<- cbind( Yq, rowSums( covar[,wheezenames] ) )
	colnames(Yq)	<- c( col0, '#wheeze' )
	Yb	<- covar[,setdiff( Ybnames, wheezenames )]

	G		<- covar[,intersect( c( admnames, 'Male', 'age', 'MX', 'smoke', smokenames, 'BMI', drugnames ), colnames( covar ) ) ]
	if( type == 'ctr' ){ 
		X		<- covar[,c( ctrnames ) ]
		print( t(sapply( 1:ncol(G), function(i){
			c( colnames(G)[i], round( summary(lm( G[,i] ~ X  ))$r.squared, 2 ) )
		})) )
	} else {
		X		<- NULL
	}

	Yb	<- apply( Yb, 2, function(y) as.numeric( y >= .5 ) )
	rownames(Yb)<- rownames(Yq)
	names(cc)		<- rownames(Yq) 

	if( type != 'raw' )
		Yq	<- scale( Yq )

	stopifnot( all( apply( Yq, 2, function( x ) length( unique( x ) ) ) > 3 ) )
	stopifnot( all( apply( Yb, 2, function( x ) length( unique( x ) ) ) <= 3) )

	N			<- length(cc)
	cases	<- which(cc==1)

	save( N, cases, cc, Yq, Yb, G, X, file=paste0( 'parsed_data/final_', type, '.Rdata' )  )
	rm(		N, cases, cc, Yq, Yb, G, X, covar )
}
sink()

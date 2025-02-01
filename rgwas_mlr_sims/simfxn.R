sim_fxn <- function(
	N, P, K,
	S, S_hom, S_het,
	sig2E, sig2hom, sig2het,
	prev=.2, asc=TRUE,
	w=rep(1/K,K),
	seed,
	N_univ=N*1e2
){

	set.seed( seed )
	stopifnot( sig2E+sig2hom+sig2het < 1 )

	# Simulate Yl, an N x K matrix of liabilities
	G       <- scale( sapply( 1:S, function(s) rbinom( N_univ, 2, runif(1,.05,.5) ) ) )
	betasl  <- gammafxn( K, S, P=1, S_hom, S_het, sig2hom=sig2hom, sig2het=sig2het, seed ) 
	Yl      <- sapply( 1:K, function(k) G %*% betasl[k,,1] + sqrt(1-sig2hom-sig2het)*rnorm(N_univ)  ) 
	Yl_std  <- sapply( 1:K, function(k) Yl[,k] - quantile(Yl[,k],1-prev*w[k]) )

	# Threshold Yl into K case subtypes
	z <- apply( Yl_std, 1, function(y) ifelse( all(y < 0), 0, which.max(y) ) )

	# ascertainment
	if( asc ){
		sub <- c( sample( which( z == 0 ), N/2 ), # ctrls
		sample( which( z != 0 ), N/2 )) # cases
	} else { 
		sub <- sample( N_univ, N )
	}
	G   <- G[sub,]
	z   <- z[sub]

	# simulate case-only phenotypes Y
	P			<- P-1  # excluding cc 
	gamma <- gammafxn( K, S, P, S_hom, S_het, sig2hom=sig2hom, sig2het=sig2het, seed )
	E     <- sqrt(sig2E) * matrix( rnorm(K*P), K, P ) 
	eps		<- sqrt(1-sig2E-sig2hom-sig2het) * matrix( rnorm( N*P ), N, P )

	Y <- matrix( NA, N, P )
	for( k in 1:K ){ 
		ksub  <- which( z == k )
		Y[ksub,]  <- G[ksub,,drop=F] %*% gamma[k,,] + matrix( E[k,], sum(z==k), P, byrow=T ) + eps[ksub,]
	} 
	Y[,P  ] <- ( Y[,P  ] > quantile(Y[,P  ],.3,na.rm=T) )
	Y[,P-1] <- ( Y[,P-1] > quantile(Y[,P-1],.7,na.rm=T) )
	cc  <- sapply( z, function(z.i) ifelse( z.i == 0, 0, 1 ) )
	list( Y=cbind( cc, Y ), G=G, z=z, E=E ) 
} 

gammafxn  <- function( K, S, P, S_hom, S_het, sig2hom, sig2het, seed ){
	set.seed( seed ) 
	gamma_hom									<- sqrt(sig2hom/S_hom)*rnorm(S_hom*P) 
	gamma <- array( 0, dim=c( K, S, P ) )
	for( k in 1:K )
		gamma[k,1:S_hom,]       <- gamma_hom
	for( k in 1:K )
		gamma[k,S_hom+1:S_het,] <- sqrt(sig2het/S_het)*rnorm(S_het*P)
	gamma
}

caus_set_fxn  <- function( S_null,S_hom,S_het,out_type )
	if( out_type == 1 ){
		if( S_null == 0 )
			return( NULL )
		(1+S_hom+S_het):(S_null+S_hom+S_het)
	} else if( out_type == 2 ){
		if( S_hom == 0 )
			return( NULL )
		1:S_hom
	} else if( out_type == 3 ){
		if( S_het == 0 )
			return( NULL )
		S_hom+1:S_het
	}

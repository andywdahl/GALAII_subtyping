library(phenix)
library(rgwas)

main <- function( Y, G, method, K ){ 
	cases <- which(Y[,1]==1)
	Yq		<- Y[,-c(1,ncol(Y)-c(1,0))]
	Yb		<- Y[,     ncol(Y)-c(1,0)] 

	if( method != 'rgwas+' ){
		Gcase <- scale(G[cases,])
		Ycase <- scale(Y[cases,-1]) 
		if(       method == 'cca-Y' ){
			pvec  <- Ycase %*% cancor( Gcase, Ycase )$ycoef[,1]
		} else if(method == 'cca-G' ){
			pvec  <- Gcase %*% cancor( Gcase, Ycase )$xcoef[,1]
		} else if(method == 'kmeans' ){
			out		<- kmeans( Ycase, centers=K )
			pvec  <- model.matrix(~as.factor(out$cluster)-1)[,-K] 
		} else if(method == 'rgwas' ){
			pvec <- mfmr( Yb=Yb[cases,], Yq=Yq[cases,], G=cbind(1,Gcase), X=matrix(rnorm(nrow(Ycase))), K=K )$pmat[,1]
		} else if(method == 'geno_pc' ){
			pvec  <- svd( Gcase )$u[,1,drop=F]
		} else if(method == 'pheno_pc' ){
			pvec  <- svd( Ycase )$u[,1,drop=F]
		} else {
			stop(method)
		} 

		Z   <- rep( '0', nrow(G) )
		Z[cases]  <- as.factor(sapply( pvec, function(u.i) ifelse( u.i>mean(pvec)-1e-8, 2, 1 ) )) 
		sapply( 1:ncol(G), function(s) run_mlr( matrix(rep(1,nrow(G))), scale(G[,s]), Z, scaleXg=FALSE ))

	} else {
		droptest_mlr( cc=Y[,1], Yb=Yb, Yq=Yq, G=cbind(1,G), X=matrix(rnorm(nrow(Y))), test_inds=1+1:ncol(G), K=K ) # , trace=T, K=K, nrun=20, init_sd=0.001
	} 
}

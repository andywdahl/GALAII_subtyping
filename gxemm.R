rm( list=ls() )
library(phenix)
library(rgwas)
library(GxEMM)
load('Rdata/setup.Rdata')

type	<- 'ctr'
phens	<- c( "asthma_onset_frombirth", "Pre.FEV1.perc.pred", "Pre.FEV1.FVC.perc.pred", "deltafev", "tIGE" ) ### quant traits=traits w reasonable GWAS
K			<- 3
qn		<- FALSE ### was already QN'd in preprocessing

grmfile   <- 'Rdata/grm.grm.raw'
if( ! file.exists(grmfile) ){
	load( 'parse_gala_data/parsed_data/final_imp.Rdata' )
	cases	<- names(cc)[cc==1]
	write.table( cbind( cases, cases ), file='Rout/keeplist.txt', col.names=F, row.names=F, quote=F )
	system( paste0( '~/ldak5.linux --bfile parse_gala_data/data/GALAII_all_freeze_041112 --calc-kins-direct Rdata/grm --keep Rout/keeplist.txt --ignore-weights YES --power -1 --kinship-raw YES' ) )
}
GRM <- as.matrix( read.table( grmfile ) )
sub	<- read.table( 'Rdata/grm.grm.id' )[,1]

for( gxemmtype in c( '', paste0( '_rand', 1:1000 ) ) )
	for( phen in phens )
{

	savefile	<-	paste0( 'Rdata/gxemm/', phen, gxemmtype, '_', type, '_qn=', qn, '.Rdata' )
	sinkfile	<-	paste0( 'Rout/gxemm/'	, phen, gxemmtype, '_', type, '_qn=', qn, '.Rout' )
	if( file.exists( savefile ) | file.exists( sinkfile ) )	next
	sink( sinkfile )
	try({

	load( paste0( 'parse_gala_data/parsed_data/final_ctr.Rdata' ) )

	Yq	<- Yq[sub,,drop=F]
	X0	<- X[sub,,drop=F]
	X0	<- X0[,-which( colSums( X0 ) == 0 ),drop=F]
	rm(X)
	X		<-G [sub,intersect( colnames(G), Gnames_cc ),drop=F]
	rm(G)

	load( 'Rdata/mfmrx_3_imp.Rdata' )
	Z	<- out$pmat
	rownames(Z)	<- names(cc)[ cc == 1 ]
	Z	<- Z[sub,]
	rm( out )

	#### fixed effect of Z + Z*X interactions
	X	<- scale(X)
	X	<- cbind( X, diag(Z[,1]) %*% X, diag(Z[,2]) %*% X )
	X	<- cbind(X,Z[,-ncol(Z)]) 
	X	<- scale(X)

	if( gxemmtype %in% paste0( '_rand', 1:1000 ) )
		Z	<- Z[sample(nrow(Z),replace=F),] 

	X	<- cbind( X, X0 ) 
	X	<- scale(X)
	y	<- scale( as.matrix(Yq[,phen]) )

	out	<- GxEMM:::full_GxEMM( y=y, X=X, K=GRM, Z=Z, binary=FALSE, tmpdir=paste0('Rdata/gxemm/tmp_',phen,gxemmtype,'_',type,'_tmp'), ldak_loc='/wynton/home/ye/andyd/ldak5.linux' )

	save( out, file=savefile )
	rm( out )

	})
	sink()
}

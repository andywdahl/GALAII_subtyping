rm( list=ls() )
load('Rdata/setup.Rdata')
sink( 'Rout/Fig3ab.Rout' )
covar	<- read.csv('parse_gala_data/data/gala2_clean2016_02_21_de_indent.csv', header=T, stringsAsFactors=F, na.strings=c('NA', 'UNK'))
covar	<- covar[-which( !is.na( covar$recontact_status ) ),] #### remove recontact records

myvec	<- c( "WBC" ,  "RBC" ,  "Neutrophils" ,  "Lymphocytes" ,  "Monocytes" ,  "Eosinophils" ,  "Basophils", 'NO.value_ppb' )
mymat	<- sapply( myvec, function(x) as.numeric(covar[,x] ) )
rownames(mymat)	<- covar[,'SubjectID']
rm( covar ) 

load('parse_gala_data/parsed_data/final_imp.Rdata' )
load( 'Rdata/mfmrx_3_imp.Rdata' )
pmat<- out$pmat
rm( out )

perm	<- 3:1
pmat	<- pmat[,perm]
subnames	<- subnames[perm] 

Z		<- rep( '0', length(cc) )
Z[which( cc == 1 )]	<- as.character(apply( pmat, 1, function(x) ifelse( max(x) > 0, which.max(x), NA ) ))
names(Z)	<- rownames(X)
Z		<- as.factor(Z)

pvs	<- rep( NA, length(myvec) )
names(pvs)	<- myvec
for( x in myvec ){ 
	qn	<- function(x) as.numeric(phenix:::quantnorm( as.matrix(x) ))
	y		<- qn(mymat[rownames(Yq),x])
	pvs[x]	<- anova(
		lm( y ~ Z  + G[,intersect( colnames(G), Gnames_cc )], subset=cases ),
		lm( y ~      G[,intersect( colnames(G), Gnames_cc )], subset=cases )
	)$Pr[2]

	pv_noadj	<- anova(
		lm( y ~ Z  , subset=cases ),
		lm( y ~ 1  , subset=cases )
	)$Pr[2]

	y		<- mymat[rownames(Yq),x]
	pv_noqn	<- anova(
		lm( y ~ Z  + G[,intersect( colnames(G), Gnames_cc )], subset=cases ),
		lm( y ~      G[,intersect( colnames(G), Gnames_cc )], subset=cases )
	)$Pr[2]

	cat( round( pvs[x]  , 4 ), x, '\n' )
	cat( round( pv_noadj, 4 ), x, '\n' )
	cat( round( pv_noqn , 4 ), x, '\n' )
	print( table(is.na(y)) ) 
}

pdf( paste0( 'figs/Fig3a.pdf' ), width=5.5, height=4.2 )
par( mar=c(1,4.7,1,1) ) 
eos	<- as.numeric(mymat[ rownames(Yq), 'Eosinophils' ] ) 
boxplot( eos ~ Z, ylab='Eosinophil %', main='', col=cols, cex.main=1.9, axes=F, xlab='', ylim=c(0,34), cex.lab=1.7 ) 
lines( c(2,4), c(32,32), lwd=2 )
text( 3, 33.5, paste0( 'p=', round( pvs['Eosinophils'], 4 ) ), cex=1.3 ) 
axis(2)
dev.off() 

pdf( paste0( 'figs/Fig3b.pdf' ), width=5.5, height=4.2 )
par( mar=c(1,4.7,1,1) )
y	<- as.numeric(mymat[ rownames(Yq), 'NO.value_ppb' ] )
boxplot( y ~ Z, ylab='FeNO', main='', col=cols, cex.main=1.9, axes=F, xlab='', ylim=c(0,172), cex.lab=1.7 )
lines( c(2,4)	, c(164,164), lwd=2 )
text(  3			, 171				, paste0( 'p=', round( pvs['NO.value_ppb'], 4 ) ), cex=1.3 )
axis(2) 
dev.off()
sink()

rm( list=ls() )
library(rgwas)
library(parallel)
library(BEDMatrix)
load('Rdata/setup.Rdata')

type	<- 'ctr'
mc.co	<- as.numeric( commandArgs(TRUE)[[1]] )

savefile	<- paste0( 'Rdata/gwas_', type, '.Rdata' )
sinkfile 	<- paste0( 'Rout/gwas_'	, type, '.Rout' )

if( file.exists( savefile ) | file.exists( sinkfile ) )	next
print( sinkfile )
sink(	sinkfile )

### load in cc, G, X
load(paste0( 'parse_gala_data/parsed_data/final_ctr.Rdata' ))
X[ X[,'centerOK']==1, ]	<- NA
X	<- X[ ,-which(colnames(X)=='centerOK')]
G	<- G[,intersect( colnames(G), Gnames_cc )]
X		<- cbind( 1, apply( cbind(G,X), 2, rgwas:::scale01 ) ) 
rm( G, Yb, Yq )

Gsnp		<- BEDMatrix( 'parse_gala_data/data/GALAII_all_freeze_041112' )
snprows	<- as.character(sapply( rownames(Gsnp), function(y) strsplit( y, '_' )[[1]][1] ))
sub			<- intersect(snprows,rownames(X))
rownames( Gsnp )	<- snprows
X				<- X		[sub,]
cc			<- cc		[sub]

out	<- mclapply( 1:ncol(Gsnp), mc.cores=mc.co, function(s){
	out.i	<- NA
	g			<- Gsnp[sub,s]
	subi	<- which(!is.na(g))
	cc		<- cc[subi]
	X			<- X [subi,]
	g			<- g [subi]
	try( out.i	<- as.numeric( summary( glm( cc ~ X + g, family=binomial(link='logit') ) )$coef['g',] ) )
	out.i
})
names(out)	<- colnames(Gsnp)

save( out, file=savefile )
print(warnings())
print('Done')
sink()
rm( out, cc, X )
}

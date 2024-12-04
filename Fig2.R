rm( list=ls() )
load('Rdata/setup.Rdata')
load( 'parse_gala_data/parsed_data/final_raw.Rdata' ) 

perm	<- 3:1 # subtype labels are arbitrary
load( 'Rdata/mfmrx_3_imp.Rdata' )
subnames<- subnames[perm]
pmat		<- out$pmat[,perm] 
rm( out, perm )

Z					<- rep( '0', length(cc) )
Z[cases]	<- as.character(apply( pmat, 1, which.max ))
names(Z)	<- rownames(X)
Z					<- as.factor(Z) 

EUR <- 1-rowSums(G[,c('AFR','NAM')])
G		<- cbind( G[,1:2], EUR, G[,-(1:2)] )

mains	<- c( 'EUR Admix %', 'Female', mains )
names(mains)[1:2]	<- c( 'EUR' , 'Female' )

mains['MX']				<- 'Mexican Ethnicity'
mains['tIGE']			<- 'IgE (log10)'
mains['deltafev']	<- 'Albuterol Effect\n(FEV Change)'

# recode to center female
G	<- cbind( 1-G[,'Male'], G )
G	<- G[,-which( colnames(G) == 'Male' )] 
colnames(G)[1]	<- 'Female'

bin_vars<- c( 'Female', 'MX', "steroids_12mo", "steroids_more2wks", "tylenol", "OTCallergy", "inhaled.steroid", "short.beta.agonist", "preg_smoke", "adults_smoke_u2yrs", "adults_smoke_3_6yrs", "adults_smoke_7plsyrs", colnames(Yb) ) 

Yq[ which(Z == '0'), c( 'asthma_onset_frombirth', 'deltafev' ) ]	<- 0 # not meaningful/measured in controls

focal_vars	<- c( "symptoms_nocold", "symptoms_12mo", "hosp_er_12mo", "steroids_more2wks", "inhaled.steroid", "asthma_onset_frombirth", "Pre.FEV1.FVC.perc.pred", "deltafev", "OTCallergy", "tIGE", "hayfever", "eczema", "rhinitis" )

cex.axis	<- 1.4
cex.main	<- 1.9
pdf( 'figs/Fig2.pdf', width=12, height=11 )
layout( matrix( c(	1:5,
										9:13,
										6:8, 14,14), 3, 5, byrow=T ) ) 
par( mar=c(3,4,5,1) ) 
for( j in focal_vars ){
	y	<- cbind( G, Yq, Yb )[,j]
	if( j %in% bin_vars ){
		barplot( tapply( y, Z, mean, na.rm=T ), main='', ylab='', col=cols, axes=F, names.arg=rep('',length(levels(Z))), ylim=c(0,1) )
		axis(2,cex.axis=cex.axis,at=0:5/5,lab=paste0(0:5/5*100,'%') )
	} else {
		ylim	<- range(y,na.rm=T)
		vioplot:::vioplot( y ~ Z, main='', col=cols, ylab='', xlab='', names=rep('',4), xaxt="n", yaxt="n", ylim=ylim, areaEqual=F  ) #axes=F, 
		if( j == 'deltafev' ){
			axis(2,cex.axis=cex.axis,at=0:4/4,lab=paste0( 100*0:4/4, '%' ) )
		} else {
			axis(2,cex.axis=cex.axis)
		}
	}
	title( main=mains[j], col=2, cex.main=cex.main, font.main=1 )
} 

par( mar=c(3,0,3,3) ) 
pie( c( sum(cc==0), colSums( pmat ) ), col=cols,lab=c( 'Controls', subnames ), cex=2.2) 
dev.off()

pdf( 'figs/SFig2.pdf', width=19, height=13 )
layout( matrix( c(1:26,26,27), 4, 7, byrow=T ) ) 
par( mar=c(3,4,5,1) ) 
for( j in setdiff( colnames( cbind( G, Yq, Yb )), focal_vars ) )try({
	y	<- cbind( G, Yq, Yb )[,j]
	if( j %in% bin_vars ){
		barplot( tapply( y, Z, mean, na.rm=T ), main='', ylab='', col=cols, axes=F, names.arg=rep('',length(levels(Z))), ylim=c(0,1) )
		axis(2,cex.axis=cex.axis,at=0:5/5,lab=paste0(0:5/5*100,'%') )
	} else {
		ylim	<- range(y,na.rm=T)
		h	<- NULL
		if( j == "current_smokers" ) h	<- 1/50 ## seems to be an edge case for vioplot
		vioplot:::vioplot( y ~ Z, main='', col=cols, ylab='', xlab='', names=rep('',4), xaxt="n", yaxt="n", ylim=ylim, areaEqual=F, h=h  ) #axes=F, 
		axis(2,cex.axis=cex.axis) 
	}
	title( main=mains[j], col=2, cex.main=cex.main, font.main=3 )
})
par( mar=c(3,0,3,3) ) 
pie( c( sum(cc==0), colSums( pmat ) ), col=cols,lab=c( 'Controls', subnames ), cex=2.2) 
dev.off()

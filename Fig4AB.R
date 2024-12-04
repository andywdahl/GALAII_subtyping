rm( list=ls() )
load('Rdata/setup.Rdata')

K				<- 3
locphens<- c( 'mlr', 'tIGE', 'asthma_onset_frombirth', 'Pre.FEV1.FVC.perc.pred', 'deltafev' )
ylabs		<- c( 'Asthma (Log OR)', 'IgE', 'Age of Onset', 'FEV1/FVC', 'Albuterol\nResponse' )
names(ylabs)<- locphens
np	<- length(locphens)

plotdata	<- 'Rdata/Fig4AB.Rdata'
if( file.exists( plotdata ) ){
	load( plotdata )
} else { 
	goodrs		<- c(
		'rs2033784', 'rs992969', 'rs2544523', ## mlr hits
		'rs12568083', 'rs4129267',						## IgE
		'rs13194312',													## AoO
		'rs7593948',													## deltaFEV 
		'rs2396255'														## FEV/FVC 
	)

	coefs	<- 
	ses		<- array(NA,dim=c(np,length(goodrs),K+1),dimnames=list(locphens,goodrs,1:(K+1)))
	pvs		<- array(NA,dim=c(np,length(goodrs),4)	,dimnames=list(locphens,goodrs,c('gwas','hom','het','glob')))

	p	<- 'mlr'
	load( 'Rdata/gwas_ctr.Rdata' )
	names( out )	<- rsids
	out	<- out[goodrs]
	coefs	[p,,1	]	<- sapply( out, function(i) i[1] )
	ses		[p,,1	]	<- sapply( out, function(i) i[2] )
	pvs		[p,,'hom']	<-
	pvs		[p,,'gwas']	<- sapply( out, function(i) i[4] )
	rm(out)

	load( paste0( 'Rdata/gwas_mlr_'	, K, '_ctr.Rdata' ) )
	names( out )	<- rsids
	out	<- out[goodrs]
	coefs	[p,,-1]	<- t(sapply( out, function(x) x$betas ))
	ses		[p,,-1]	<- t(sapply( out, function(x) x$betas.se ))
	pvs		[p,,'het'	]	<- sapply( out, function(x) 10^-x$log10p_het )
	pvs		[p,,'glob']	<- sapply( out, function(x) 10^-x$log10p_glob )
	rm(out)

	for( p in locphens[-1] ){
		load( paste0( 'Rdata/gwas_caseonly_'	, p, '_'	, K, '_ctr.Rdata' ) )
		out	<- out[goodrs]
		coefs	[p,,-1]	<- t(sapply( out, function(x) x$coef[c('Z1:g','Z2:g','Z3:g'),1] ))
		ses		[p,,-1]	<- t(sapply( out, function(x) x$coef[c('Z1:g','Z2:g','Z3:g'),2] ))
		pvs		[p,,'hom']	<- sapply( out, function(x) x$homp0)
		pvs		[p,,'het'	]	<- sapply( out, function(x) x$hetp )
		pvs		[p,,'glob']	<- sapply( out, function(x) x$globp)
		rm(out)
	}
	save( goodrs, pvs, coefs, ses, file=plotdata )
}

perm	<- 3:1
subnames	<- subnames[perm] 
coefs	<- coefs[,,c(1,1+perm)]
ses		<- ses	[,,c(1,1+perm)]

for( i in goodrs ) ### polarize by hom effect on asthma
	coefs[,i,]	<- coefs[,i,] * sign(coefs['mlr',i,1]) 

pdf( 'figs/Fig4AB.pdf', width=9.4, height=5.0 )
layout( matrix( 1:6, 1, 6, byrow=T ), widths=c(1.4,6.5,2.5,2.5,2.5,4.0) )
par( mar=rep(0,4) )
par( mar=c(14.5,0,4,0) )
plot.new()
mtext(2,line=-4.1,text='SNP Effects on Asthma\nand Secondary Traits', cex=1.25)

snptypes	<- c( 'Yes', 'Yes', 'Yes', 'No', 'Yes', 'No', 'No', 'No', 'No', 'No' )
names(snptypes)	<- goodrs
par( mar=c(14.5,2.5,4,1.9) )
for( p in locphens[c(1,3,5,4,2)] ){
	if( p == locphens[1] ){
		is	<- c( 'rs992969', 'rs2033784', 'rs2544523' )
	} else if( p == locphens[2] ){
		is	<- c( 'rs4129267', 'rs12568083' )
	} else if( p == locphens[3] ){
		is	<- c( 'rs13194312' )
	} else if( p == locphens[4] ){
		is	<- c( 'rs2396255' )
	} else if( p == locphens[5] ){
		is	<- c( 'rs7593948' )
	}

	xs	<- rep( 1:length(is), each=K+1 )+rep( c( -.27, -.08+.29*seq(0,1,length=K) ), length(is) )*1.1

	ys	<- c(t(coefs[p,is,]))
	yse	<- c(t(ses	[p,is,]))

	ylim	<- range( c( 0, ys-yse, ys+yse ), na.rm=T )
	if( p != 'mlr' ){
		xlim	<- range(xs[-1])	+ c(-1,1)*.10
	} else {
		xlim	<- range(xs)			+ c(-1,1)*.05 
	}

	plot( range(xs), ylim, axes=F, type='n', xlab='', ylab='', main='', xlim=xlim ) #, ylab=ylabs[p]
	mtext(3,line=1,text=ylabs[p], cex=0.9)
	axis(2)

	for( i in is ){
		x.i	<- which( is==i )
		axis(1,at=x.i,lab=i, cex.axis=1.30, tick=F, line=.3 )
		for( ii in 1:3 )
		axis(1,at=x.i,lab=format( pvs[p,i,1+ii]	, digits=1, scientific=T ), cex.axis=1.35,lwd=0,line=1.0+ii*1.8,tick=F)
		if( p == 'mlr' ){
		axis(1,at=x.i,lab='--'																						, cex.axis=1.35,lwd=0,line=1.0+ 4*1.8,tick=F)
		} else {
		axis(1,at=x.i,lab=format( pvs['mlr',i,'gwas'], digits=1, scientific=T ), cex.axis=1.35,lwd=0,line=1.0+ 4*1.8,tick=F)
		}
		axis(1,at=x.i,lab=snptypes[i], cex.axis=1.35,lwd=0,line=11.1,tick=F)
	}
	abline( h=0, col=2, lty=3, lwd=2 )
	points(xs, ys, pch=16, cex=3																, col=rep(cols,length(is)) )
	arrows(xs, ys-yse, xs, ys+yse, length=0.05, angle=90, code=3, col=rep(cols,length(is)), lwd=2.5 )
	rm(is)
}
mtext(1,outer=T,line=-10.7,adj=.01,text=expression( p[hom] )	, cex=1.05)
mtext(1,outer=T,line=-8.7	,adj=.01,text=expression( p[het] )	, cex=1.05)
mtext(1,outer=T,line=-6.7	,adj=.01,text=expression( p[glob])	, cex=1.05)
mtext(1,outer=T,line=-4.7	,adj=.01,text=expression( p[asthma]), cex=1.05) 
mtext(1,outer=T,line=-1.6	,adj=.003,text='Known GWAS \nHit for Asthma',cex=0.77)

dev.off()

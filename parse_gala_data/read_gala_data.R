rm( list=ls() )
savefile	<- paste0( 'parsed_data/covar.Rdata' )
if( file.exists( savefile ) ) stop('Already built parsed_data/covar.Rdata' )
sink( 'parsed_data/read_gala_data.Rout' )

covar	<- read.csv('data/gala2_clean2016_02_21_de_indent.csv', header=T, stringsAsFactors=F, na.strings=c('NA', 'UNK'))
colnames( covar )
colnames( covar )[ grep( '2wk', colnames( covar ) ) ]
dim( covar )

#### remove recontact records
reconts	<- which( !is.na( covar$recontact_status ) )
covar	<- covar[-reconts,]
dim( covar )

#### remove NA cc
cc_NA	<- which( is.na( covar$asthma ) )
covar	<- covar[-cc_NA,]
dim( covar )

##### remove ever-pregnant
table( covar[,'pregnant'], exclude=NULL )
dim( covar )
covar	<- covar[-which( !is.na( covar[,'pregnant'] ) ),]
dim( covar )

phenonames = c(
	'SubjectID', 'case',
	'hosp_er_12mo', 'asthma_doc',
	'steroids_12mo', 'steroids_more2wks', 'OTCallergy', 'tylenol', 'short.beta.agonist', 'inhaled.steroid',
	'Male', 'height_cm', 'weight_kg', 'bmicatObese', 'bmicatOverweight', 'bmicatUnderweight', 'BMI',
	'tIGE',
	'hayfever', 'rash', 'eczema', 'sinusitis', 'rhinitis',
	'preg_smoke', 'adults_smoke_u2yrs', 'adults_smoke_3_6yrs', 'adults_smoke_7plsyrs', 'current_smokers',
	'wake_asthma_freq_1wk', 'asthma_morning_1wk', 'asthma_lmtd_1wk', 'asthma_shortbreath_1wk', 'wheez_1wk',
	'Pre.FVC.perc.pred', 'Pre.FEV1.perc.pred', 'Pre.FEV1.FVC.perc.pred', "Pre.PEF.perc.pred", 'Pre.FEF.perc.pred', 'maxfvc', 'maxfev1','deltafev',
	'asthma_duration', 'asthma_onset_frombirth',
	'age', 'symptoms_12mo', 'sleep_trouble_2wks', 'symptoms_nocold',
	'wheez_cough_2wks', 'wheez_cough_weather', 'wheez_cough_pollen', 'wheez_cough_cold', 'wheez_cough_activity', 'wheez_cough_housedust', 'wheez_cough_pets', 'wheez_cough_wind', 'wheez_cough_odors',
	'AFR','NAM'
)

### code sex, BMI, ethn, center as factors
center	<- covar[,'Center']
sum(center=='BP')
covar	<- covar[-which(center=='BP'),]

center	<- covar[,'Center']
center	<- model.matrix( ~ -1 + center )[,-1]
apply( center, 2, table )

covar[,'Male']<- as.numeric( covar[,'Male'] == 'Male' )

table( bmicat	<- as.factor( covar[,'bmicat'] ) )
bmifac	<- model.matrix( ~ -1 + bmicat )
nabmi		<- which( colnames(bmifac) %in% c( "bmicatNo BMI Pctile Value", "bmicatNo BMI Value" ) )
bmifac[which( rowMeans( bmifac[,nabmi] ) > 0 ),]	<- NA
bmifac	<- bmifac[,-nabmi]
bmifac	<- bmifac[,-1]
colSums( is.na( bmifac ) )
colSums( bmifac, na.rm=T )

ethn	<- t(sapply( covar[,'child.ethnicity'], function(eth.i)
	if( eth.i == 'Puerto Rican' ){
		c(1,0)
	} else if( eth.i == 'Mexican' ){
		c(0,1)
	} else { ## combines 'Mixed Latino' and 'Other Latino'
		c(0,0)
	}
))
colnames( ethn )	<- c( 'PR', 'MX' )
apply( ethn, 2, table )
table( covar[,'child.ethnicity'] )

### add functionals
maxfev1			<- apply( covar[,c("Pre.FEV1.Meas","Post1.FEV1.Meas","Post2.FEV1.Meas")], 1, function(x) ifelse( all(is.na(x)), NA, max(x,na.rm=T) ))
maxfvc			<- apply( covar[,c("Pre.FVC.Meas","Post1.FVC.Meas","Post2.FVC.Meas")]		, 1, function(x) ifelse( all(is.na(x)), NA, max(x,na.rm=T) ))
maxfev1			<- apply( covar[,c("Pre.FEV1.Meas","Post1.FEV1.Meas","Post2.FEV1.Meas")], 1, function(x) ifelse( all(is.na(x)), NA, max(x,na.rm=T) ))
maxfvc			<- apply( covar[,c("Pre.FVC.Meas","Post1.FVC.Meas","Post2.FVC.Meas")]		, 1, function(x) ifelse( all(is.na(x)), NA, max(x,na.rm=T) ))
deltafev		<- apply( covar[,c("Pre.FEV1.Meas","Post2.FEV1.Meas")], 1, function(x) (x[2] - x[1])/x[1] )

Pre.PEF.perc.pred	<- covar[,'Pre.PEFR.Meas']			/covar[,'PEF.Pred']*100
Pre.FEF.perc.pred	<- covar[,'Pre.FEF.25.75.Meas']	/covar[,'FEF.Pred']*100

asthma_doc	<- apply( covar[,c("asthma_doc_sick","asthma_doc_nosick")]							, 1, function(x) as.numeric( any(x,na.rm=T) ) )
hosp_er_12mo<- apply( covar[,c("asthma_er_12mo","asthma_hospital_12mo")]						, 1, function(x) as.numeric( any(x,na.rm=T) ) )
apply( covar[,c("asthma_er_12mo","asthma_hospital_12mo")]						, 1, function(x) as.numeric( any(x,na.rm=T) ) )
case				<- as.numeric( covar[,'asthma'] )
BMI					<- as.numeric( covar[,'weight_kg'] / ( covar[,'height_cm'] / 100 )^2 )

covar				<- cbind( covar, bmifac, maxfev1, maxfvc, deltafev, asthma_doc, case, hosp_er_12mo, BMI, Pre.PEF.perc.pred, Pre.FEF.perc.pred	)

### subset to phenonames
covar	<- covar[ , phenonames ] 
print( colMeans( is.na( covar ) ))

### add messy factors
cov_rownames	<- covar[,1]
covar					<- covar[,-1]
covar					<- cbind( covar, ethn, center )
covar					<- apply( covar, 2, as.numeric )

### remove measure that are high missing or zero variance
badcols	<- which( apply( covar, 2, function(x) mean(!is.na(x))<.05 ) )
stopifnot( length(badcols) == 0 )

badcols	<- which( apply( covar, 2, function(x) sd(x,na.rm=T)==0 ) )
stopifnot( length(badcols) == 0 )

cc		<- covar[,'case']
cases	<- which( cc == 1 )
round( sort( colMeans( is.na( covar					) ) ), 4 )
round( sort( colMeans( is.na( covar[cases,] ) ) ), 4 )

colnames( covar )
dim( covar )

rownames(covar)	<- cov_rownames

#### subset to 'PR' and 'MX', which loses roughly 1K 'mixed'/'other'
print( dim( covar ) )
cc		<- cc		[ -which( rowSums( covar[,c( 'PR', 'MX' )] ) == 0 ) ]
covar	<- covar[ -which( rowSums( covar[,c( 'PR', 'MX' )] ) == 0 ),]
covar	<- covar[,-which( colnames(covar) == 'PR' 								)	]
print( dim( covar ) )

save( cc, covar, file=savefile	)
warnings()
sink()

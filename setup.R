rm( list=ls() )

drugnames		<- c( "steroids_12mo", "steroids_more2wks", "tylenol", "OTCallergy", "inhaled.steroid", "short.beta.agonist" ) ### all correspond to last 12 months
ctrnames		<- c( "centerCH", "centerCM", 'centerHR', "centerLC", "centerNW", "centerNY", "centerOK", "centerRI", "centerSF", 'centerSJ', "centerTX", "centerVA" ) ### only 3 obs: "centerBP", 
wheezenames	<- c( "wheez_cough_weather", "wheez_cough_pollen", "wheez_cough_cold", "wheez_cough_activity", "wheez_cough_housedust", "wheez_cough_pets", "wheez_cough_wind", "wheez_cough_odors" )
admnames		<- c( "AFR", "NAM" )
smokenames	<- c( "preg_smoke", "adults_smoke_u2yrs", "adults_smoke_3_6yrs", "adults_smoke_7plsyrs", "current_smokers" )
Gnames_cc		<- c( admnames, 'Male', 'age', 'MX', 'smoke', smokenames, 'BMI' )

Ynames	<- c( 
	"asthma_onset_frombirth",
	"wake_asthma_freq_1wk", "asthma_lmtd_1wk", "wheez_1wk", "wheez_cough_2wks", "asthma_morning_1wk", "asthma_shortbreath_1wk",
	"Pre.FEV1.perc.pred", "Pre.FEV1.FVC.perc.pred", "deltafev", 
	"tIGE",
	'#wheeze'
)

Ybnames	<- c(
	"symptoms_12mo",
	"hayfever", "eczema", "rash", "sinusitis", "rhinitis",
	"hosp_er_12mo",
	"symptoms_nocold",
	"sleep_trouble_2wks",
	wheezenames
)

mains	<- c( 'AFR Admix %', 'NAM Admix %', 'Male', 'Age', 'Mexican\nEthnicity', 'Smoke', 'Adult Smoke\nwhile Pregnant', 'Adult Smoke\nwhile <2yo', 'Adult Smoke\nwhile 2-6yo', 'Adult Smoke\nwhile >=7yo', 'Adult Smoke\nCurrent',	'BMI', 
	"Steroids\nLast Year", "Steroids for\n>2 Weeks", "Tylenol", "OTC Allergy", "Inhaled Steroid", "Short Beta Agonist",
	"Age of Onset",
	"Wake from symptoms\n(Freq past week)", "Limited by symptoms\n(Freq past week)", "Wheezing\n(Freq past week)", "Wheeze or cough\n(Freq past 2 weeks)", "Symptoms in morning\n(Freq past week)", "Short of breath\n(Freq past week)",
	'FEV1', "FEV1/FVC", "Albuterol\nResponse", "IgE", '# Wheeze Causes', "Symptoms\nLast Year", "Hay Fever", "Eczema", "Rash", "Sinusitis", "Rhinitis", "Exacerbation\nLast Year",
	"Symptoms w/o\na Cold", "sleep_trouble_2wks", wheezenames
)
names( mains )	<- c( Gnames_cc, drugnames, Ynames, Ybnames ) 

cols	<- c( '#CED7D8', '#FF652F', '#14A76C', '#5680E9' )
subnames	<- c( 'Th2+, Atopy+', 'Th2+, Atopy-', 'Th2-' )

save.image( 'Rdata/setup.Rdata' )


library(BEDMatrix)
bfile		<- 'parse_gala_data/data/GALAII_all_freeze_041112'
Gsnp		<- BEDMatrix( bfile )
longrsids	<- colnames(Gsnp)
rsids	<- sapply( colnames(Gsnp), function(y) strsplit( y, '_' )[[1]][1] )
mafs	<- sapply( 1:ncol(Gsnp), function(i){ print(i); x	<- mean(Gsnp[,i],na.rm=T)/2; ifelse( x > 1/2, 1-x, x ) })
miss	<- sapply( 1:ncol(Gsnp), function(i){ print(i); mean(is.na(Gsnp[,i]) ) })
names( mafs ) <- rsids
names( miss ) <- rsids
rm( Gsnp )

save.image( 'Rdata/setup.Rdata' )

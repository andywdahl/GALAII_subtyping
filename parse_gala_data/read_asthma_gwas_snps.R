x	<- read.csv( 'efotraits_EFO_0000270-associations-2020-03-2.csv' )
pv	<- sapply( as.character(x[,'P.value']), function(ii){
	ii	<- strsplit( ii, ' x ' )[[1]]
	ii1	<- as.numeric( ii[1] )
	ii2	<- as.numeric( strsplit( ii[2], '-' )[[1]][2] )
	ii1 * 10^-ii2
})

snpraw		<- as.character(x[,1])
ast_snps	<- sapply( strsplit( snpraw, '-' ), function(x) x[1] )
ast_snps	<- unique( ast_snps )

length( ast_snps )
save( ast_snps, file='../Rdata/ast_snps_all.Rdata' )

relpops	<- which( x[,'Study.accession'] %in% c( 
'GCST003831',
'GCST005212',
'GCST003832',
'GCST003833',
'GCST007266',
'GCST001182',
'GCST002445',
'GCST009435',
'GCST004107',
'GCST009434',
'GCST002519',
'GCST000548',
'GCST004108',
'GCST002552'
))

head(x[relpops,])
snpraw		<- as.character(x[relpops,1])
ast_snps	<- sapply( strsplit( snpraw, '-' ), function(x) x[1] )
ast_snps_info		<- cbind( ast_snps, x[relpops,-c(1,3,4,5,6,7,8,9,10,12)] )
ast_snps	<- unique( ast_snps )
length(ast_snps)
dim(ast_snps_info)
dim(x[relpops,])

save( ast_snps, relpops, file='../Rdata/ast_snps.Rdata' ) 

load( '../Rdata/setup.Rdata' )
sub	<- which( ast_snps_info[,1] %in% rsids )
dim( ast_snps_info)
sub
ast_snps_info	<- ast_snps_info[sub,]
dim( ast_snps_info)


write.table( ast_snps_info, file='../figs/STabX.txt', quote=F, row.names=F )
unique( ast_snps_info[,3] )
sort(table( ast_snps_info[,3] ))

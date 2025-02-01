rm( list=ls() )
Ns			<- c( 1e3		,3e3		,1e4	,1e3			,3e3			,1e4			) #, 1e3				,3e3				,1e4				)#, 3e4   
types		<- c( 'N1e3','N3e3','N1e4','1e3_K1'	,'3e3_K1'	,'1e4_K1'	) #,'1e3_noasc','3e3_noasc','1e4_noasc')#,'N3e4' 
ascs		<- c( TRUE  ,TRUE  ,TRUE  ,TRUE			,TRUE			,TRUE			) #,FALSE			,FALSE			,FALSE			)#,TRUE		
Ks			<- c( 2			,2     ,2     ,2				,1			  ,1				) #,2					,2					,2					)#,2			
names(Ns)		<- types
names(Ks)		<- types
names(ascs) <- types

sig2homs	<- c( .04, .004 )
sig2hets	<- c( .004, .04 )

S			<- 12
P			<- 30

maxit	<- 1000

methods			<- c( 'rgwas', 'rgwas+', 'kmeans', 'geno_pc', 'pheno_pc', 'cca-Y' )
nicemethods	<- c( 'RGWASX', 'RGWAS', 'k-means', 'Genotype PC1', 'Phenotype PC1', 'Geno-Pheno CC1' )
Me					<- length( methods )

cols	= 1:length(methods)
names( cols )	<- methods

save.image( 'Rdata/setup.Rdata' )

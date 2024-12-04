#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o Rout/out			#-- output directory (fill in)
#$ -e Rout/err                     #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r n                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=6G                  #-- submits on nodes with enough free memory (required)
#$ -l h_rt=4:20:00                #-- runtime limit (see above; this requests 24 hours)
#$ -t 2-4													#-- remove first '#' to specify the number of
#$ -pe smp 12

date
hostname

### find subtypes
Rscript run_mfmr.R $SGE_TASK_ID 10

### estimate covariate effect sizes
Rscript run_mfmrx_betas.R $SGE_TASK_ID 10
Rscript run_mlrx_betas.R  $SGE_TASK_ID 

### test covariate effect sizes
Rscript run_droptest.R $SGE_TASK_ID 12
Rscript run_droptest_mlr.R  $SGE_TASK_ID 12

### fit alternate subtypes with kmeans
Rscript run_kmeans.R $SGE_TASK_ID 10

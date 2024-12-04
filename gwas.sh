#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o Rout/gwas			#-- output directory (fill in)
#$ -e Rout/err                     #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r n                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=4G                  #-- submits on nodes with enough free memory (required)
#$ -l h_rt=46:20:00                #-- runtime limit (see above; this requests 24 hours)
#$ -t 2-4													#-- remove first '#' to specify the number of
#$ -pe smp 30

date
hostname

Rscript gwas.R	30
Rscript gwas_mlr.R			$SGE_TASK_ID 30
Rscript gwas_caseonly.R	$SGE_TASK_ID 30

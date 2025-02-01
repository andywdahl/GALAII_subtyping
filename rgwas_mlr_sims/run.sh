#!/bin/bash
#
#$ -S /bin/bash
#$ -o ./Rout  ## where to put standard output (to screen)
#$ -e ./Rout  ## where to put error messages
#$ -cwd            #-- start in current working directory
#$ -l h_data=1G      #### 5G good for 2K samples,
#$ -l time=18:20:00       #-- runtime limit (this requests 24 hours)
#$ -l highp
#$ -t 1:1000

export OMP_NUM_THREADS=1 

Rscript run.R ${SGE_TASK_ID}

### module load R/3.6.1; qsub -cwd -V run.sh


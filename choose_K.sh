#!/bin/bash
#
#$ -S /bin/bash
#$ -o ./Rout/findK  ## where to put standard output (to screen)
#$ -e ./Rout/findK  ## where to put error messages
#$ -cwd            #-- start in current working directory
#$ -l mem_free=3G                  #-- submits on nodes with enough free memory (required)
#$ -l h_rt=88:20:00                #-- runtime limit (see above; this requests 24 hours)
#$ -t 1-6												#-- remove first '#' to specify the number of
#$ -pe smp 20

Rscript choose_K.R $SGE_TASK_ID 20

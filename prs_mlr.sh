#!/bin/bash
#
#$ -S /bin/bash
#$ -o ./Rout/prs  ## where to put standard output (to screen)
#$ -e ./Rout/prs  ## where to put error messages
#$ -cwd            #-- start in current working directory
#$ -l mem_free=.5G                  #-- submits on nodes with enough free memory (required)
#$ -l h_rt=0:20:00                #-- runtime limit (see above; this requests 24 hours)
#$ -t 1-100												#-- remove first '#' to specify the number of

Rscript prs_mlr.R $SGE_TASK_ID

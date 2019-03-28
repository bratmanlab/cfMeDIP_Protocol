#!/bin/bash
#$ -S /bin/bash
#$ -cwd

module load R

R CMD BATCH ~/qsub/R_files/normal_cfmedip_sample_19_figure.R ~/qsub/R_files/normal_cfmedip_sample_19_figure.Rout

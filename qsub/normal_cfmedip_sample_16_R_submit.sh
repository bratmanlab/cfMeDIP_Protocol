#!/bin/bash
#$ -S /bin/bash
#$ -cwd

module load R

R CMD BATCH ~/qsub/R_files/normal_cfmedip_sample_16_figure.R ~/qsub/R_files/normal_cfmedip_sample_16_figure.Rout

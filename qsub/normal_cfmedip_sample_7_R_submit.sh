#!/bin/bash
#$ -S /bin/bash
#$ -cwd

module load R

R CMD BATCH /cluster/projects/scottgroup/questionable_people/justin/Projects/cfMeDIP_Protocol/qsub/R_files/normal_cfmedip_sample_7_figure.R /cluster/projects/scottgroup/questionable_people/justin/Projects/cfMeDIP_Protocol/qsub/R_files/normal_cfmedip_sample_7_figure.Rout

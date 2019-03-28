#!/bin/bash
#$ -S /bin/bash
#$ -cwd

######################################################
## Bash Script Generator for R Analysis in Figure 4 ##
######################################################

# March 21st 2019
#
# Data used was from 21 representative cfMeDIP-seq libraries generated from healthy donors
#
# BAM files were generated using bwa-mem and duplicates were removed by samtools
#
# This bash script generates .R files for each BAM file in the BAM directory. The .R files are then ran in parallel to generate the RData files used in all described analysis
#
# Run this script as: sh generate_scripts_figure4.sh bam_directory project_directory

## DEFINE DIRECTORIES ## 
########################

bam_dir=$1
project_dir=$2
qsub_dir=$project_dir/qsub
output_dir=$project_dir/output

## CREATE SAMPLE COUNTER ##
###########################

count=1

## CREATE OUTPUT DIRECTORY ##
#############################

mkdir $output_dir
mkdir $qsub_dir

R_dir=$qsub_dir/R_files

mkdir $R_dir

## LIST EACH SAMPLE ##
######################

for i in $bam_dir 
do

identifier=$(echo normal_cfmedip_sample_$count) # i.e. normal_cfmedip_sample_1

(( count++ ))

mkdir $output_dir/$identifier

#####################
## CREATE R SCRIPT ##
#####################

## Load R Packages ##
#####################

echo -e "library(GenomicRanges)" > $R_dir/$identifier"_figure.R"
echo -e "library(GenomicAlignments)" >> $R_dir/$identifier"_figure.R"
echo -e "library(BSgenome.Hsapiens.UCSC.hg19)" >> $R_dir/$identifier"_figure.R"
echo -e "library(DeCarvalho)" >> $R_dir/$identifier"_figure.R"
echo -e "library(dplyr)" >> $R_dir/$identifier"_figure.R" 
echo -e "library(Repitools)\n" >> $R_dir/$identifier"_figure.R"

## Define Variables ##
######################

echo -e "allmainchrs <- paste(\"chr\",1:22,sep=\"\")\n" >> $R_dir/$identifier"_figure.R"

## Read BAM File and Select Fragments with Insert Size of 167 bp (Insert Size Mode) ##
######################################################################################

echo -e "$identifier"_gr" <- GRanges(readGAlignmentPairs(\"$i\"))" >> $R_dir/$identifier"_figure.R"
echo -e "$identifier"_gr" <- $identifier"_gr"[seqnames($identifier"_gr") %in% allmainchrs]" >> $R_dir/$identifier"_figure.R"
echo -e "$identifier"_widths" <- width($identifier"_gr")" >> $R_dir/$identifier"_figure.R"
echo -e "$identifier"_widths" <- table($identifier"_widths")" >> $R_dir/$identifier"_figure.R"
echo -e "$identifier"_mode" <- $identifier"_widths"[$identifier"_widths" == max($identifier"_widths")]" >> $R_dir/$identifier"_figure.R"
echo -e "$identifier"_index" <- as.numeric($identifier"_mode")" >> $R_dir/$identifier"_figure.R"
echo -e "$identifier"_mode" <- as.numeric(names($identifier"_mode"))" >> $R_dir/$identifier"_figure.R"
echo -e "$identifier"_mode_gr" <- $identifier"_gr"[width($identifier"_gr") == $identifier"_mode"]" >> $R_dir/$identifier"_figure.R"

## Calculate the Proportion of Fragments in each Normal Patient Library with n CpGs ##
######################################################################################

echo -e "$identifier"_mode_gr"\$CpG_count <- cpgDensityCalc($identifier"_mode_gr", Hsapiens)\n" >> $R_dir/$identifier"_figure.R"

echo -e "$identifier"_mode_gr_cpgprop" <- table($identifier"_mode_gr"\$CpG_count) / sum(table($identifier"_mode_gr"\$CpG_count))" >> $R_dir/$identifier"_figure.R"
echo -e "save($identifier"_mode_gr_cpgprop", file=\"$output_dir/$identifier/$identifier"_mode_gr_cpgprop.RData"\")\n" >> $R_dir/$identifier"_figure.R"

## Annotate Windows based on Overlap with CpG Islands, Shores, Shelfs, and Open Sea ##
######################################################################################

echo -e "$identifier"_mode_gr"\$ID <- paste(seqnames($identifier"_mode_gr"), start($identifier"_mode_gr"), end($identifier"_mode_gr"), sep=\".\")" >> $R_dir/$identifier"_figure.R"
echo -e "$identifier"_mode_df" <- as.data.frame($identifier"_mode_gr")" >> $R_dir/$identifier"_figure.R"
echo -e "$identifier"_mode_df" <- Epigenome.hg19($identifier"_mode_df", is.CpG=T)" >> $R_dir/$identifier"_figure.R"
echo -e "$identifier"_mode_df_output" <- $identifier"_mode_df"[,c('seqnames','start','end','cgi')]" >> $R_dir/$identifier"_figure.R"
echo -e "colnames($identifier"_mode_df_output") <- c('chr','start','end','cpg_annot')" >> $R_dir/$identifier"_figure.R"
echo -e "save($identifier"_mode_df_output", file=\"$output_dir/$identifier/$identifier"_mode_df_R.RData"\")\n" >> $R_dir/$identifier"_figure.R"

## Randomly Subset Genome-Wide 167 bp Windows based on the Number of CpG Fragments in Sample ##
###############################################################################################

echo -e "chr_allmainchrs_167 <- genomeBlocks(seqlengths(Hsapiens)[allmainchrs],width=$identifier"_mode")\n" >> $R_dir/$identifier"_figure.R"

echo -e "set.seed(0)" >> $R_dir/$identifier"_figure.R"
echo -e "index <- sample(1:length(chr_allmainchrs_167), nrow($identifier"_mode_df_output"))" >> $R_dir/$identifier"_figure.R"
echo -e "chr_allmainchrs_167 <- chr_allmainchrs_167[index]" >> $R_dir/$identifier"_figure.R"
echo -e "chr_allmainchrs_167\$ID <- paste(seqnames(chr_allmainchrs_167), start(chr_allmainchrs_167), end(chr_allmainchrs_167), sep=\".\")" >> $R_dir/$identifier"_figure.R"
echo -e "chr_allmainchrs_df_167 <- as.data.frame(chr_allmainchrs_167)" >> $R_dir/$identifier"_figure.R"
echo -e "chr_allmainchrs_df_167 <- Epigenome.hg19(chr_allmainchrs_df_167, is.CpG=T)" >> $R_dir/$identifier"_figure.R"
echo -e "chr_allmainchrs_df_167_output <- chr_allmainchrs_df_167[,c('seqnames','start','end','cgi')]" >> $R_dir/$identifier"_figure.R"
echo -e "colnames(chr_allmainchrs_df_167_output) <- c('chr','start','end','cpg_annot')\n" >> $R_dir/$identifier"_figure.R"

## Calculate Observed over Expected Occurrence of CpG Annotations in Sample ##
##############################################################################

echo -e "$identifier"_mode_df_obsvsexp" <- table($identifier"_mode_df_output"\$cpg_annot) / table(chr_allmainchrs_df_167_output\$cpg_annot)" >> $R_dir/$identifier"_figure.R"
echo -e "save($identifier"_mode_df_obsvsexp", file=\"$output_dir/$identifier/$identifier"_mode_df_obsvsexp.RData"\")" >> $R_dir/$identifier"_figure.R"

###############################
## CREATE BASH SUBMIT SCRIPT ##
###############################

echo -e "#!/bin/bash\n#$ -S /bin/bash\n#$ -cwd\n" > $qsub_dir/$identifier"_R_submit.sh"
echo -e "module load R\n" >> $qsub_dir/$identifier"_R_submit.sh"
echo -e "R CMD BATCH $R_dir/$identifier"_figure.R" $R_dir/$identifier"_figure.Rout"" >> $qsub_dir/$identifier"_R_submit.sh"

done

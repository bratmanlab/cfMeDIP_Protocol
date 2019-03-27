###############
## Figure 4 ##
###############

# March 21st 2019
# 
# Data used was from 21 representative cfMeDIP-seq libraries generated from healthy donors
#
# BAM files were generated using bwa-mem and duplicates were removed by samtools
#
# Scripts were generated for each cfMeDIP-seq library in bash

###############################################################################################
# R Script that incorporates output from Bash generated R scripts for each cfMeDIP-seq library #
###############################################################################################

# Load R packages

library(GenomicAlignments)
library(GenomicRanges)
library(Repitools)
library(BSgenome.Hsapiens.UCSC.hg19)

## FIGURE 4A ##
###############

# Define variables

allmainchrs <- paste("chr",1:22,sep="")

# Load .RData files generated from Bash script

setwd("/cluster/projects/scottgroup/questionable_people/justin/Projects/cfMeDIP_Protocol")

files <- list.files("output",full.names = T) # cpgprop.RData files are located in their respective folders within the directory "output", output/sample_name/...

## List the full path of files from each of the 21 healthy donors, within their respective folders, load each file
for(i in files){
  tmp <- list.files(i, pattern="cpgprop", full.names = T)
  load(tmp)
}

files2 <- ls(pattern="cpgprop") # list R variable names for each object

# Generate a dataframe with each sample per row, and the number of CpGs per column (0 to 20)

df <- as.data.frame(matrix(nrow=0, ncol=21)) # empty dataframe to append each cpgprop.RData

for(i in files2){
  df <- rbind(df, get(i)[1:21])
}

rownames(df) <- files2
colnames(df) <- 0:20

# Divide genome into 167 bp bins, excluding chromosome X and Y

chr_gw_167 <- genomeBlocks(seqlengths(Hsapiens)[allmainchrs], width=167)

chr_gw_167$CpG_count <- cpgDensityCalc(chr_gw_167, Hsapiens)

# Randomly subset genome 167 bp windows to the average number of 167 fragments among the 21 healthy donors

set.seed(0)

sample_index <- sample(1:length(chr_gw_167),1567561) # 1567561 is the average number of 167 bp fragments among the 21 healthy donors

chr_gw_167 <- chr_gw_167[sample_index]

# Calculate the proportion of randomly subset 167 bp windows with n CpGs

chr_gw_167_cpgprop <- table(chr_gw_167$CpG_count) / sum(table(chr_gw_167$CpG_count))

# Generate plot for Figure 4a 

par(ps=7, lwd=1)

pdf("/cluster/projects/scottgroup/questionable_people/justin/Projects/cfMeDIP_Protocol/R/Figure_4a.pdf")
boxplot(df,
        ylim=c(0,0.45),
        xlim=c(1,21),
        ylab="Fraction of all windows",
        xlab="# of CpGs",
        names=0:20,
        las=2,
        cex=1,
        lwd=1,
        pch="")
lines(x=1:21, y=as.numeric(chr_gw_167_cpgprop)[1:21], lwd=1, lty=2)

title("# of CpGs in 167 bp Genome-Wide Windows and cfMeDIP Fragments", font.main=1)

legend("topright",lty=2, legend=c("Genome-wide"), bty='n')

dev.off()

## FIGURE 4B ##
###############

rm(list=ls())

# Define variables

allmainchrs <- paste("chr",1:22,sep="")

# Load obsvsexp.RData files in R/bam_file/... directories

files <- list.files("output", full.names = T)
for(i in files){
  tmp <- list.files(i, pattern="obs_vs_exp", full.names = T)
  load(tmp)
}


files2 <- ls(pattern="obs_vs_exp") # list R variable names for each object

# Create dataframe with each of the 21 healthy donors as rows, and each CpG annotation for columns (open sea, shelf, shore, island)

df <- as.data.frame(matrix(nrow=0,ncol=4))

for(i in files2){
  df <- rbind(df, get(i)[c(1,4,3,2)])
}

colnames(df) <- c("Island","Shore","Shelf","Open sea")
rownames(df) <- files2

# Generate plot for Figure 4b

par(ps=7)

pdf("/cluster/projects/scottgroup/questionable_people/justin/Projects/cfMeDIP_Protocol/R/Figure_4b.pdf")
boxplot(df,
        ylab="Observed/expected overlap of fragments",
        cex=1,
        lwd=1,
        pch="",
        xlab="Features")
abline(h=1, lwd=1, lty=2)

title("CpG annotation of 167 bp fragments", font.main=1)

dev.off()

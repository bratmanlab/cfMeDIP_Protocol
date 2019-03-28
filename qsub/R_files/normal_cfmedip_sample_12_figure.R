library(GenomicRanges)
library(GenomicAlignments)
library(BSgenome.Hsapiens.UCSC.hg19)
library(DeCarvalho)
library(dplyr)
library(Repitools)

allmainchrs <- paste("chr",1:22,sep="")

normal_cfmedip_sample_12_gr <- GRanges(readGAlignmentPairs("~/bam_directory/file.bam"))
normal_cfmedip_sample_12_gr <- normal_cfmedip_sample_12_gr[seqnames(normal_cfmedip_sample_12_gr) %in% allmainchrs]
normal_cfmedip_sample_12_widths <- width(normal_cfmedip_sample_12_gr)
normal_cfmedip_sample_12_widths <- table(normal_cfmedip_sample_12_widths)
normal_cfmedip_sample_12_mode <- normal_cfmedip_sample_12_widths[normal_cfmedip_sample_12_widths == max(normal_cfmedip_sample_12_widths)]
normal_cfmedip_sample_12_index <- as.numeric(normal_cfmedip_sample_12_mode)
normal_cfmedip_sample_12_mode <- as.numeric(names(normal_cfmedip_sample_12_mode))
normal_cfmedip_sample_12_mode_gr <- normal_cfmedip_sample_12_gr[width(normal_cfmedip_sample_12_gr) == normal_cfmedip_sample_12_mode]
normal_cfmedip_sample_12_mode_gr$CpG_count <- cpgDensityCalc(normal_cfmedip_sample_12_mode_gr, Hsapiens)

normal_cfmedip_sample_12_mode_gr_cpgprop <- table(normal_cfmedip_sample_12_mode_gr$CpG_count) / sum(table(normal_cfmedip_sample_12_mode_gr$CpG_count))
save(normal_cfmedip_sample_12_mode_gr_cpgprop, file="~/output/normal_cfmedip_sample_12/normal_cfmedip_sample_12_mode_gr_cpgprop.RData")

normal_cfmedip_sample_12_mode_gr$ID <- paste(seqnames(normal_cfmedip_sample_12_mode_gr), start(normal_cfmedip_sample_12_mode_gr), end(normal_cfmedip_sample_12_mode_gr), sep=".")
normal_cfmedip_sample_12_mode_df <- as.data.frame(normal_cfmedip_sample_12_mode_gr)
normal_cfmedip_sample_12_mode_df <- Epigenome.hg19(normal_cfmedip_sample_12_mode_df, is.CpG=T)
normal_cfmedip_sample_12_mode_df_output <- normal_cfmedip_sample_12_mode_df[,c('seqnames','start','end','cgi')]
colnames(normal_cfmedip_sample_12_mode_df_output) <- c('chr','start','end','cpg_annot')
save(normal_cfmedip_sample_12_mode_df_output, file="~/output/normal_cfmedip_sample_12/normal_cfmedip_sample_12_mode_df_R.RData")

chr_allmainchrs_167 <- genomeBlocks(seqlengths(Hsapiens)[allmainchrs],width=normal_cfmedip_sample_12_mode)

set.seed(0)
index <- sample(1:length(chr_allmainchrs_167), nrow(normal_cfmedip_sample_12_mode_df_output))
chr_allmainchrs_167 <- chr_allmainchrs_167[index]
chr_allmainchrs_167$ID <- paste(seqnames(chr_allmainchrs_167), start(chr_allmainchrs_167), end(chr_allmainchrs_167), sep=".")
chr_allmainchrs_df_167 <- as.data.frame(chr_allmainchrs_167)
chr_allmainchrs_df_167 <- Epigenome.hg19(chr_allmainchrs_df_167, is.CpG=T)
chr_allmainchrs_df_167_output <- chr_allmainchrs_df_167[,c('seqnames','start','end','cgi')]
colnames(chr_allmainchrs_df_167_output) <- c('chr','start','end','cpg_annot')

normal_cfmedip_sample_12_mode_df_obsvsexp <- table(normal_cfmedip_sample_12_mode_df_output$cpg_annot) / table(chr_allmainchrs_df_167_output$cpg_annot)
save(normal_cfmedip_sample_12_mode_df_obsvsexp, file="~/output/normal_cfmedip_sample_12/normal_cfmedip_sample_12_mode_df_obsvsexp.RData")

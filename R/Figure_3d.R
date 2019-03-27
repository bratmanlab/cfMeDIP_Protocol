###############
## Figure 3d ##
###############

# March 21st 2019
# 
# Data used was from a representative cfMeDIP-seq library generated from a healthy donor
#
# Insert size metrics were obtained from column 9 of the bam file via bash script, considering only the first mate
# "cat file.bam | awk '$9 > 0 {print $9}' > file_insertmetrics.txt"
#
# Insert size metrics were collected for both bam files aligned with bwa-mem and bowtie2

## Read insert size metric files...

bwa_is_file <- "path_to_bwa_file.txt/bwa_insertmetric.txt" 
bwa_is <- read.table(bwa_is_file)[,1]

bowtie_is_file <- "path_to_bowtie_file.txt/bowtie_insertmetric.txt"
bowtie_is <- read.table(bowtie_is_file)[,1]

## Count the instances of each insert size

bwa_is_tab <- table(bwa_is)
bowtie_is_tab <- table(bowtie_is)

## For each possible insert size, select the maximum number of reads from either bwa or bowtie2

max_frag_dis <- c() # create an empty vector to store values at each insert size
index <- 1 # counter for each insert size
i <- 1 # counter for next coordinate in empty vector

while(index <= 500){
  if(index %in% as.numeric(names(bwa_is_tab)) & index %in% as.numeric(names(bowtie_is_tab))){
    tmp <- max(c(bwa_is_tab[which(as.numeric(names(bwa_is_tab)) == index)], bowtie_is_tab[which(as.numeric(names(bowtie_is_tab)) == index)]))
    max_frag_dis <- c(max_frag_dis, tmp)
    names(max_frag_dis)[i] <- index
    i <- i + 1
  } else {
    if(index %in% as.numeric(names(bwa_is_tab))){
      tmp <- bwa_is_tab[which(as.numeric(names(bwa_is_tab)) == index)]
      max_frag_dis <- c(max_frag_dis, tmp)
      names(max_frag_dis)[i] <- index
      i <- i + 1
    } else {
      if(index %in% as.numeric(names(bowtie_is_tab))){
        tmp <- bowtie_is_tab[which(as.numeric(names(bowtie_is_tab)) == index)]
        max_frag_dis <- c(max_frag_dis, tmp)
        names(max_frag_dis)[i] <- index
        i <- i + 1
      }
    }
  }
  index <- index + 1
}

## Determine percentage of inserts that fall into boundaries of first (60 - 255 bp) and second peak (256 - 450 bp)

peak1_sum <- sum(max_frag_dis[which(names(max_frag_dis) == "60"):which(names(max_frag_dis) == "255")])
peak2_sum <- sum(max_frag_dis[which(names(max_frag_dis) == "256"):which(names(max_frag_dis) == "450")])

peak1_percent <- round(peak1_sum / (peak1_sum + peak2_sum), digits=3) * 100
peak2_percent <- round(peak2_sum / (peak1_sum + peak2_sum), digits=3) * 100

## Generate plot for Figure 3d

par(ps=7)

plot(x=as.numeric(names(max_frag_dis)), y=max_frag_dis, 
     type='l',
     xlab="Fragment size (bp)",
     ylab="Count",
     xaxt="n")  
axis(1,at=seq(0,500,50))

title("Fragment size distribution of cfMeDIP-seq library\nwith BWA-mem and Bowtie2 merged", font.main=1)    

abline(v=c(60,255,450), lty=2, lwd=0.5, col=rgb(0,0,0,1/4)) # visualize boundaries

lines(x=c(167,167), lwd=0.5, y=c(-10000, (57237)), lty=3) # define mode of distribution
# cfMeDIP-seq Protocol
*Author: Justin M. Burgener*

In-depth protocol for cfMeDIP-seq library prepartion, QC, and bioinformatic analysis.

## *About*

The code within the bratmanlab/cfMeDIP_Protocol repository is used to demonstrate the insert-size distribution, CpG enrichment, and mapping annotation of reads from representative cfMeDIP-seq libraries. This code uses proccessed data from BAM files (FASTQ to SAM by bwa-mem, SAM to BAM by samtools), requiring bash programming and pre-installed software. Access to a high-performance computing cluster is also recommended. Due to privacy concerns, the BAM files used are not available, however the script used to generate the processed .RData files in subsequent analysis is provided (generate_scripts_figure4.sh).

The output from the "generate_scripts_figure4.sh" script produces directories for each sample (within the directory qsub) with respective R scripts. The resulting .RData files generated for each sample are used to generatae the plots in Figure 4. PDF plots for Figure 4a and Figure 4b are also provided to show the expected output (assuming the same seed is set).

Below are the steps described to generate the plots shown (**NOTE: All scripts should be run in your project directory**):

## *Requirements:*
  * Computer running a Linux system (≥ 8 GB RAM)
    * Modules: bwa (version 0.7.15), bowtie2 (version 2.2.6), samtools (version 1.3.1), igenome-human/hg19
  * R/RStudio (version 3.5 or greater)
    * Packages: BSgenome.hsapiens.UCSC.hg19, GenomicRanges, AnnotationHub, Repitools

## *Procedure:*

  1. Align FASTQ files to reference genome (bwa/0.7.15, igenome-human/hg19)
     ```bash
     bwa mem -M -t4 $BWAINDEX $R1.fastq $R2.fastq > $file.sam
     ```
  2. Convert SAM file to BAM file (samtools/1.3.1)
     ```bash 
     samtools view -buS -f 2 -F 4 $file.sam | samtools sort -@4 -o $file.bam
     samtools index $file.bam
     ```
  3. Remove duplicates
     ```bash
     samtools rmdup $file.bam $bam_directory/$rmdup_file.bam
     samtools index $bam_directory/$rmdup_file.bam
     ```
  4. Create R scripts with generate_scripts_figure4.sh (R/3.5) (**NOTE: The output of this script will generate a qsub directory containing all of the R submission scripts and .R code [cfmedip_sample_n_figure.R]**)
     ```bash
     sh generate_scripts_figure4.sh $bam_directory $project_directory
     ```
  5. Generate QC metrics from each sample
     ```bash
     cd $project_directory/$qsub
     for i in *.sh
     do
     qsub $i
     done
     ```
  6. Performed grouped analysis of QC metrics from samples for Figure 4.
     ```bash
     R CMD BATCH Figure_4.R Figure_4.Rout
     ```
  7. From a representative sample, determine insert-size distribution of reads for Figure 3d. (**NOTE: Because bwa-mem removes reads outside of a specified interquartile range, bowtie2 (bowtie2/2.2.6) was also used to align the same sample, keeping reads with insert-sizes ranging from 50 to 750 bp**)
     ```
     cat file.bam | awk '$9 > 0 {print $9}' > file_insertmetrics.txt # repeat for bam files aligned with bowtie2
     
     R CMD BATCH Figure_3d.R Figure_3d.Rout
     ```
For further clarification of the bioinformatic analysis, please e-mail Justin Burgener at justin.burgener@uhnresearch.ca

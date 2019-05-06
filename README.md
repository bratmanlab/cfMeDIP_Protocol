# cfMeDIP-seq Protocol
In-depth protocol for cfMeDIP-seq library prepartion, QC, and bioinformatic analysis.

The code within the bratmanlab/cfMeDIP_Protocol repository is used to demonstrate the insert-size distribution, CpG enrichment, and mapping annotation of reads from representative cfMeDIP-seq libraries. This code uses proccessed data from BAM files (FASTQ to SAM by bwa-mem, SAM to BAM by samtools), requiring bash programming and pre-installed software. Access to a high-performance computing cluster is also recommended. Due to privacy concerns, the BAM files used are not available due to privacy concerns, however the script used to generate the processed .RData files in subsequent analysis is provided (generate_scripts_figure4.sh).

The output from the "generate_scripts_figure4.sh" script produces directories for each sample (within the directory qsub) with respective R scripts. The resulting .RData files generated for each sample were used to generatae the plots in Figure 3d and Figure 4. Due to file size constraints, insert-size metrics of each sample can be available upon request. PDF plots for Figure 4a and Figure 4b are also provided to show the expect output (assuming the same seed is set). 

## bash
Contains bash script used to generate independent R scripts, and their respective bash submission scripts, for each sample. Output of bash scripts are located in the qsub directory

## qsub
Contains the bash submission scripts and their independent R scripts for each sample. Jobs were submitted on a Unix cluster and performed in parallel. Outputs of each job were placed within the output directory

## output
Contains RData objects for each sample analyzed. 

### sample_cpgprop.RData
Proportion of n CpGs [0 - 20] within each 167 bp fragment.

### sample_obs_vs_exp.RData
Number of fragments with specified CpG annotation (i.e. Island, Open Sea, Shelf, Shore) over the expected number of CpG annotations based on random 167 bp genome-wide windows.

### sample_output.RData
Processed data used to calculate 'sample_obs_vs_exp.RData'. Columns include chromosome, start, end, and CpG annotation for each fragment.

## R 
Contains R scripts used to generate Figure 3d and Figure 4. Due to file size constraints, insert-size metrics of each sample is available upon request. PDF plots for Figure 4a and Figure 4b are also provided to show the expected output.

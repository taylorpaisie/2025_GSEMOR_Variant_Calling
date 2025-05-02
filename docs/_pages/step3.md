---
layout: default
title: "3. Trimming and Filtering"
nav_order: 3
parent: Tutorial
permalink: /step3/
---


# Step 3: Trimming and Filtering

1. Trimming the bad quality reads from our fastq files
    * In the previous episode, we took a high-level look at the quality of each of our samples using FastQC  
    * We visualized per-base quality graphs showing the distribution of read quality at each base across all reads in a sample and extracted information about which samples fail which quality checks  
    * It is very common to have some quality metrics fail, and this may or may not be a problem for your downstream application  
    * For our variant calling workflow, we will be removing some of the low quality sequences to reduce our false positive rate due to sequencing error  
    * We will use a program called [Trimmomatic](https://github.com/timflutre/trimmomatic) to filter poor quality reads and trim poor quality bases from our samples  

```
$ trimmomatic  
    Usage:  
       PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...  
   or:  
       SE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...  
   or:  
       -version 

```


   * This output shows us that we must first specify whether we have paired end (PE) or single end (SE) reads  
   * Next, we specify what flag we would like to run. For example, you can specify threads to indicate the number of processors on your computer that you want Trimmomatic to use  
   * In most cases using multiple threads (processors) can help to run the trimming faster  
   * These flags are not necessary, but they can give you more control over the command  
   * The flags are followed by positional arguments, meaning the order in which you specify them is important  
   * In paired end mode, Trimmomatic expects the two input files, and then the names of the output files  
   * While, in single end mode, Trimmomatic will expect 1 file as input, after which you can enter the optional settings and lastly the name of the output file  

   * Trimmomatic command options:  
     * ILLUMINACLIP: Cut adapter and other Illumina-specific sequences from the read  
     * SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average
  quality within the window falls below a threshold  
     * LEADING: Cut bases off the start of a read, if below a threshold quality  
     * TRAILING: Cut bases off the end of a read, if below a threshold quality  
     * CROP: Cut the read to a specified length  
     * HEADCROP: Cut the specified number of bases from the start of the read 4
     * MINLEN: Drop the read if it is below a specified length  
     * TOPHRED33: Convert quality scores to Phred-33  
     * TOPHRED64: Convert quality scores to Phred-64  

1. Running Trimmomatic:  
   
    * Move to the correct directory with untrimmed fastq files we downloaded:  
    `$ cd ~/variant_calling/data/untrimmed_fastq` 

    * Copy Illumina adapters from Trimmomatic into working directory:  
    `$ cp /home/gitpod/miniconda/envs/variant_calling/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa .`  

    * Run Trimmomatic:  

    ```
    trimmomatic PE \
    SRR1972917_1.fastq.gz SRR1972917_2.fastq.gz \
    SRR1972917_1.trim.fastq.gz SRR1972917_1un.trim.fastq.gz \
    SRR1972917_2.trim.fastq.gz SRR1972917_2un.trim.fastq.gz \
    SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
    ```


    * Output:  

    ```
    TrimmomaticPE: Started with arguments:
    SRR1972917_1.fastq.gz SRR1972917_2.fastq.gz SRR1972917_1.trim.fastq.gz SRR1972917_1un.trim.fastq.gz SRR1972917_2.trim.fastq.gz SRR1972917_2un.trim.fastq.gz SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
    Multiple cores found: Using 4 threads
    Using PrefixPair: 'AGATGTGTATAAGAGACAG' and 'AGATGTGTATAAGAGACAG'
    Using Long Clipping Sequence: 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
    Using Long Clipping Sequence: 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
    Using Long Clipping Sequence: 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
    Using Long Clipping Sequence: 'CTGTCTCTTATACACATCTGACGCTGCCGACGA'
    ILLUMINACLIP: Using 1 prefix pairs, 4 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
    Quality encoding detected as phred33
    Input Read Pairs: 4377867 Both Surviving: 1241328 (28.35%) Forward Only Surviving: 2133670 (48.74%) Reverse Only Surviving: 18680 (0.43%) Dropped: 984189 (22.48%)
    TrimmomaticPE: Completed successfully
    ```


   * List out files created by Trimmomatic:  
    `$ ls SRR1972917*`  

   * Trimmed files should be smaller in size than our untrimmed fastq files
  
2. Running a for loop on all fastq files  

    ```
    for infile in *_1.fastq.gz
        do
            base=$(basename ${infile} _1.fastq.gz)
            trimmomatic PE ${infile} ${base}_2.fastq.gz \
            ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \
            ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \
            SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
        done
    ```



3. Moving trimmed fastq files to a new directory:
    * We have now completed the trimming and filtering steps of our quality control process! Before we move on, let’s move our trimmed FASTQ files to a new subdirectory within our data/ directory  
  
    `$ cd ~/variant_calling/data/untrimmed_fastq`  
    `$ mkdir -p ~/variant_calling/data/trimmed_fastq`  
    `$ mv *trim* ~/variant_calling/data/trimmed_fastq`  
    `$ cd ~/variant_calling/data/trimmed_fastq`  
    `$ ls -al`  

4. Lets rerun FastQC on the trimmed fastq files  
    `$ fastqc *trim.fastq.gz`

---

👉 Ready? Continue to [Step 4 - Reference Based Mapping]({{ "/step4/" | relative_url }})
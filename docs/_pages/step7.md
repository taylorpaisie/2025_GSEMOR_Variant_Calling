---
layout: default
title: "7. Automating a Variant Calling Workflow"
nav_order: 7
parent: Tutorial
permalink: /step7/
---

# Step 7: Automating a Variant Calling Workflow

#### You wrote a simple shell script eariler in this tutorial
#### Shell scripts can be much more complicated than that and can be used to perform a large number of operations on one or many files  
#### This saves you the effort of having to type each of those commands over for each of your data files and makes your work less error-prone and more reproducible  
#### For example, the variant calling workflow we just carried out had about eight steps where we had to type a command into our terminal   
#### If we wanted to do this for all eight of our data files, that would be forty-eight steps  
#### If we had 50 samples, it would be 400 steps
#### We have also used for loops earlier to iterate one or two commands over multiple input files  
#### In these for loops, the filename was defined as a variable in the for statement, which enables you to run the loop on multiple files  


#### Running Trimmomatic on all samples:

```
$ for infile in *_1.fastq.gz
    do
        base=$(basename ${infile} _1.fastq.gz)  
        trimmomatic PE ${infile} ${base}_2.fastq.gz  
        ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz  
        ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz  
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15   
    done
``` 

* Within the Bash shell you can create variables at any time (as we did above, and during the `for` loop lesson)  
* Assign any name and the value using the assignment operator: `=`  
* You can check the current definition of your variable by typing into your script: echo `$variable_name`  

#### Notice that in this for loop, we used two variables, `infile`, which was defined in the for statement, and `base`, which was created from the filename during each iteration of the loop  

#### We can extend these principles to the entire variant calling workflow. 
#### To do this, we will take all of the individual commands that we wrote before, put them into a single file, add variables so that the script knows to iterate through our input files and write to the appropriate output files

#### Make a directory for automated variant calling script:

`$ mkdir scripts`  
`$ cd scripts/`  
`$ touch run_variant_calling.sh`  

Instead of using the `touch` command you can use `code` to open in VS code:  

`$ code run_variant_calling.sh`



#### Our variant calling workflow has the following steps:  
1. Index the reference genome for use by bwa and samtools  
2. Align reads to reference genome  
3. Convert the format of the alignment to sorted BAM, with some intermediate steps  
4. Calculate the read coverage of positions in the genome  
5. Detect the SNVs  
6. Filter and report the SNVs in VCF  


#### Let's write our script, which should look like this:  

```
#! bin/bash/

set -e
cd ~/variant_calling/results

genome=~/variant_calling/data/ref_genome/KJ660346.2.fasta

bwa index $genome

# makes directories (should already be made)
# mkdir -p sam bam bcf vcf

for fq1 in ~/variant_calling/data/trimmed_fastq/*_1.trim.fastq.gz
    do
    echo "working with file $fq1"

    base=$(basename $fq1 _1.trim.fastq.gz)
    echo "base name is $base"

    fq1=~/variant_calling/data/trimmed_fastq/${base}_1.trim.fastq.gz
    fq2=~/variant_calling/data/trimmed_fastq/${base}_2.trim.fastq.gz
    sam=~/variant_calling/results/sam/${base}.aligned.sam
    bam=~/variant_calling/results/bam/${base}.aligned.bam
    sorted_bam=~/variant_calling/results/bam/${base}.aligned.sorted.bam
    raw_bcf=~/variant_calling/results/bcf/${base}_raw.bcf
    variants=~/variant_calling/results/vcf/${base}_variants.vcf
    final_variants=~/variant_calling/results/vcf/${base}_final_variants.vcf 

    bwa mem $genome $fq1 $fq2 > $sam
    samtools view -S -b $sam > $bam
    samtools sort -o $sorted_bam $bam 
    samtools index $sorted_bam
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 
    vcfutils.pl varFilter $variants > $final_variants
   
    done
```
 

#### We change our working directory so that we can create new results subdirectories in the right location  

`$ cd ~/variant_calling/results`   

#### Next we tell our script where to find the reference genome by assigning the genome variable to the path to our reference genome:  

`$ genome=~/variant_calling/data/ref_genome/KJ660346.2.fasta`  

#### Now we want to index the reference genome for BWA:

`$ bwa index $genome`  

#### We already have all the directories in `results` that are needed


#### What is in our automated variant calling script?
* We will use a loop to run the variant calling workflow on each of our FASTQ files   
* The full list of commands within the loop will be executed once for each of the FASTQ files in the data/trimmed_fastq/ directory   
* We will include a few echo statements to give us status updates on our progress  
* The first thing we do is assign the name of the FASTQ file we are currently working with to a variable called fq1 and tell the script to echo the filename back to us so we can check which file we are on  

```
for fq1 in ~/variant_calling/data/trimmed_fastq/*_1.trim.fastq.gz
    do
    echo "working with file $fq1"
```  

* We then extract the base name of the file (excluding the path and .fastq.gz extension) and assign it to a new variable called base  

```
    base=$(basename $fq1 _1.trim.fastq.gz)
    echo "base name is $base"
```

* We can use the base variable to access both the base_1.fastq.gz and base_2.fastq.gz input files, and create variables to store the names of our output files   
* This makes the script easier to read because we do not need to type out the full name of each of the files: instead, we use the base variable, but add a different extension (e.g. .sam, .bam) for each file produced by our workflow  

```
    .# input fastq files
    fq1=~/variant_calling/data/trimmed_fastq/${base}_1.trim.fastq.gz
    fq2=~/variant_calling/data/trimmed_fastq/${base}_2.trim.fastq.gz
        
    # output files
    sam=~/variant_calling/results/sam/${base}.aligned.sam
    bam=~/variant_calling/results/bam/${base}.aligned.bam
    sorted_bam=~/variant_calling/results/bam/${base}.aligned.sorted.bam
    raw_bcf=~/variant_calling/results/bcf/${base}_raw.bcf
    variants=~/variant_calling/results/vcf/${base}_variants.vcf
    final_variants=~/variant_calling/results/vcf/${base}_final_variants.vcf   
```

#### Now to the actual steps of the workflow:

1. Align the reads to the reference genome and output a .sam file:

`$ bwa mem $genome $fq1 $fq2 > $sam`  

2. Convert the SAM file to BAM format:  

`$ samtools view -S -b $sam > $bam`  

3. Sort the BAM file:  

`$ samtools sort -o $sorted_bam $bam`  

4. Index the BAM file for display purposes:  

`$ samtools index $sorted_bam`  

5. Calculate the read coverage of positions in the genome:  

`$ bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam`  

6. Call SNVs with bcftools:  

`$ bcftools call --ploidy 1 -m -v -o $variants $raw_bcf`  

7. Filter and report the SNVs in a VCF:  

`$ vcfutils.pl varFilter $variants > $final_variants`  


#### Time to run the script!

`$ bash run_variant_calling.sh`  

#### Now your automated variant calling script should be running!!!
#### Tip:  using echo statements within your scripts is a great way to get an automated progress update
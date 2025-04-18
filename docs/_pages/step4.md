---
layout: default
title: "4. Reference Based Mapping "
nav_order: 4
parent: Tutorial
permalink: /step4/
---

# Step 4: Reference Based Mapping  

1. Aligning to a reference genome
    * We perform read alignment or mapping to determine where in the genome our reads originated from  
    * There are a number of tools to choose from and, while there is no gold standard, there are some tools that are better suited for particular NGS analyses  
    * We will be using the [Burrows Wheeler Aligner (BWA)](https://github.com/lh3/bwa), which is a software package for mapping low-divergent sequences against a large reference genome  
    * The alignment process consists of two steps:  
        * Indexing the reference genome  
        * Aligning the reads to the reference genome  


2. Downloading reference genome  
   * Navigate to NCBI and search for GenBank accession `KJ660346.2` and download fasta file  

    `$ mkdir -p ~/variant_calling/data/ref_genome`  
    `$ cd ~/variant_calling/data/ref_genome`  
    `$ wget -nv https://raw.githubusercontent.com/taylorpaisie/VEME_2024_NGS_Variant_Calling/main/KJ660346.2.fasta -O KJ660346.2.fasta`  
    

3. Create directories for the results that will be generated as part of this workflow    
    * We can do this in a single line of code, because mkdir can accept multiple new directory names as input  
    `$ cd ~/variant_calling/`  
    `$ mkdir -p results/sam results/bam results/bcf results/vcf`  

4. Index the reference genome  
    * Our first step is to index the reference genome for use by BWA  
    * Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment  
    * Indexing the reference only has to be run once  
    * The only reason you would want to create a new index is if you are working with a different reference genome or you are using a different tool for alignment  
    `$ bwa index data/ref_genome/KJ660346.2.fasta`  

5. Align reads to the reference genome  
    * The alignment process consists of choosing an appropriate reference genome to map our reads against and then deciding on an aligner  
    * We will use the BWA-MEM algorithm, which is the latest and is generally recommended for high-quality queries as it is faster and more accurate  
    `$ bwa mem data/ref_genome/KJ660346.2.fasta data/trimmed_fastq/SRR1972917_1.trim.fastq.gz data/trimmed_fastq/SRR1972917_2.trim.fastq.gz > results/sam/SRR1972917.aligned.sam`  

    * The SAM file, is a tab-delimited text file that contains information for each individual read and its alignment to the genome  
    * The compressed binary version of SAM is called a BAM file  
    * We use this version to reduce size and to allow for indexing, which enables efficient random access of the data contained within the file  
    * The file begins with a header, which is optional  
    * The header is used to describe the source of data, reference sequence, method of alignment, etc., this will change depending on the aligner being used  
    * Following the header is the alignment section  
    * Each line that follows corresponds to alignment information for a single read  
    * Each alignment line has 11 mandatory fields for essential mapping information and a variable number of other fields for aligner specific information  
    * An example entry from a SAM file is displayed below with the different fields highlighted  
    * We will use the program [samtools](https://www.htslib.org/) for these steps
  
6. Convert SAM file to BAM format  
    `$ samtools view -S -b results/sam/SRR1972917.aligned.sam > results/bam/SRR1972917.aligned.bam`  

7. Sort BAM file by coordinates  
    `$ samtools sort -o results/bam/SRR1972917.aligned.sorted.bam results/bam/SRR1972917.aligned.bam`   

    * Lets take a look at the statistics in our BAM file:  
    `$ samtools flagstat results/bam/SRR1972917.aligned.sorted.bam`  

Now lets view the header of our BAM file:  

`$ samtools view -H results/bam/SRR1972917.aligned.sorted.bam`

```
@HD	VN:1.5	SO:coordinate
@SQ	SN:KJ660346.2	LN:18959
@PG	ID:bwa	PN:bwa	VN:0.7.18-r1243-dirty	CL:bwa mem /Users/tpaisie/variant_calling/data/ref_genome/KJ660346.2.fasta /Users/tpaisie/variant_calling/data/trimmed_fastq/SRR1972917_1.trim.fastq.gz /Users/tpaisie/variant_calling/data/trimmed_fastq/SRR1972917_2.trim.fastq.gz
@PG	ID:samtools	PN:samtools	PP:bwa	VN:1.20	CL:samtools view -S -b /Users/tpaisie/variant_calling/results/sam/SRR1972917.aligned.sam
@PG	ID:samtools.1	PN:samtools	PP:samtools	VN:1.20	CL:samtools sort -o /Users/tpaisie/variant_calling/results/bam/SRR1972917.aligned.sorted.bam /Users/tpaisie/variant_calling/results/bam/SRR1972917.aligned.bam
@PG	ID:samtools.2	PN:samtools	PP:samtools.1	VN:1.20	CL:samtools view -H results/bam/SRR1972917.aligned.sorted.bam
```

    
Now lets view a small part of our BAM file:  

`$ samtools view results/bam/SRR1972917.aligned.sorted.bam | head -1` 


Mapping Quality (MAPQ) and Compact Idiosyncratic Gapped Alignment Representation (CIGAR):

`$ samtools view results/bam/SRR1972917.aligned.sorted.bam | cut -f 5,6 | grep -v "101" | head -5`

```
60	32M
60	33M
49	31M
60	2S98M
60	3S98M
```

The values in the MAPQ (Mapping Quality) column here are all 60. This column was designed to indicate the likelihood of the alignment being placed incorrectly. It is the same Phred score that we encountered in the FASTQ files. And we read it the same way, 60/10 = 6 so the chance of seeing this alignment being wrong is 10^-6 or 1/1,000,000 one in a million.

The numbers that an aligner puts into the MAPQ field are typically estimates. It is not possible to mathematically compute this value. What this field does is to inform us on a guess by the aligner’s algorithm. This guess is more of a hunch that should not be treated as a continuous, numerical value. Instead, it should be thought of as an ordered label, like “not so good” or “pretty good”. Aligner developers even generate unique MAPQ qualities to mark individual cases. For example, bwa will create a MAPQ=0 if a read maps equally well to more than one location.

The CIGAR string is a different beast altogether. It is meant to represent the alignment via numbers followed by letters:

- `M` match of mismatch
- `I` insertion
- `D` deletion
- `S` soft clip
- `H` hard clip
- `N` skipping

These are also meant to be “readable”; the 2S98M says that 2 bases are soft clipped and the next 98 are a match or mismatch.

The CIGAR representation is a neat concept, alas it was developed for short and well-matching reads. As soon as the reads show substantial differences the CIGAR representations are much more difficult to read.

There are also different variants of CIGAR. The “default” CIGAR encoding describes both the match and mismatch in the same way, with an M. There are some unfortunate rationales and explanations for having adopted this choice - one that complicates analysis even more.

The extended CIGAR encoding adopted by others uses the symbols of = and X to indicate matches and mismatches.


---

👉 Ready? Continue to [Step 5 - Variant Calling](./step5.md)
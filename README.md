# Identifying differentially expressed transposons across four life-cycle stages of *Fasciola hepatica*

## Introduction
Transposable elements (TEs) are highly repetitive mobile sequences, which play diverse roles in genome regulation. As well, it is expected that TEs participate in lncRNAs function. In trematodes, lncRNA might be involved in development processes and life cycle regulation. For future studies it is significant to explore connection between TEs expression and developmental stages.

## Mission
Our study is devoted to detecting transposons in transcriptomes of four different life-cycle stages of *F. hepatica* and analyzing their differential expression across stages.

## Goals
1. Create mapping of RNA-seq data onto a list of transposons sequences 
2. Quantify expression from mappings for each transposon
2. Normalize TE data using statistical approach
4. Discover differentially expressed TEs on each stage of *F. hepatica* life cycle

## Methods

### *F. hepatica* RNA-seq data
For our purposes we used public data which can be found in SRA NCBI archive.  
Acsessions: [ERX535560](https://www.ncbi.nlm.nih.gov/sra/?term=ERX535560), [ERX535563](https://www.ncbi.nlm.nih.gov/sra/?term=ERX535563), [ERX535561](https://www.ncbi.nlm.nih.gov/sra/?term=ERX535561), [ERX535360](https://www.ncbi.nlm.nih.gov/sra/?term=ERX535360), [ERX535363](https://www.ncbi.nlm.nih.gov/sra/?term=ERX535363), [ERX535362](https://www.ncbi.nlm.nih.gov/sra/?term=ERX535362), [ERX535371](https://www.ncbi.nlm.nih.gov/sra/?term=ERX535371), [ERX535370](https://www.ncbi.nlm.nih.gov/sra/?term=ERX535370).  
Info table: [SraRunInfo](https://github.com/NickPanyushev/IB_Fasciola/blob/master/SraRunInfo.csv)

### List of *F. hepatica* TEs
As well we applied fasta file with repeatitive elements of *F.hepatica* genome, created earlier with help of RepeatMasker tool.

### Project workflow 
Short outline of our work:
![workflow](https://github.com/NickPanyushev/IB_Fasciola/blob/master/some%20useful%20pictures/project_workflow.png)

First of all, row reads of 4 developmental stages were filtered using Trimmomatic and FastQC for quality control. Differential expression analysis was performed with two separate methods: kallisto + sleuth and TEtools + DeSeq2. To analyze any TEs that show change in expression across the different stages, likelihood ratio test (LRT) was implemented. In the end lists of significant TEs (FDR>0.05), which are differentially expressed across developmental stages, was obtained. 
All steps can be found in [lab notebook](https://github.com/NickPanyushev/IB_Fasciola/blob/master/lab_notebook.md)

### System requirements
* python v.3.6
* R v.3.4.4
* Trimmomatic v.0.36
* FastQC v.0.11.7
* kallisto v.0.44.0 
* sleuth R package v.0.29.0
* TEtools v.3 
* DeSeq2 R package v.1.18.1

## Results
Both approaches to DE analysis show similiar results.
All results can be found in repo:
[sleuth](https://github.com/NickPanyushev/IB_Fasciola/tree/master/sleuth_res)
[DeSeq2](https://github.com/NickPanyushev/IB_Fasciola/tree/master/DeSeq2_res)

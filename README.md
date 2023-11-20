![BSB_test_repo](https://github.com/CPalmer3200/bsb_test_repo/assets/145576128/4c8bf560-0113-4e4a-8619-e11b5b7b15f6)

# Results repository for BSB coding test

## 1. Repository architecture

This repository contains Nextflow and Python scripts designed to process DNA double strand breaks (DSBs) reads. These scripts are designed to be run on a Linux-based/WSL system through bash scripts. 

### 1.1 Nextflow pipeline

The Nextflow pipeline is comprised of three scripts designed to be run in the order: 1. index_pipeline.nf, 2. mapping_pipeline.nf, 3. intersect_pipeline.nf


- index_pipeline.nf: This script retrieves the chromosome.fasta file stored in the data directory and indexes it using Burrows-Wheeler Aligner (BWA). Please            install and activate BWA in the current environment before running this script with:

```
nextflow run index_pipeline.nf
```


- mapping_pipeline.nf: This script is designed to map the sample.fastq.gz files to the indexed chromosome with bwa mem and convert these into .bam files with           samtools. The resulting .bam files will then be converted to .bed files using bedtools and the reads will be trimmed to include just the break sites by               accessing the trim_reads.py script in this repo. Please install and activate BWA, samtools, Python3 and bedtools before running this script with:

```
nextflow run mapping_pipeline.nf
```


![BSB_test_repo](https://github.com/CPalmer3200/bsb_test_repo/assets/145576128/4c8bf560-0113-4e4a-8619-e11b5b7b15f6)

# Results repository for the BSB coding test (CPalmer3200)

## 1. Repository architecture

This repository contains Nextflow and Python scripts designed to process and analyse samples containing DNA double-strand breaks (DSBs). These scripts are designed to be run on a Linux-based/WSL system through bash scripts. These scripts were constructed with Nextflow version 23.10.0.5889 and Python version 3.12.0. Please make sure these are installed prior to running the scripts. 

### 1.1 Nextflow pipeline

The Nextflow pipeline is comprised of three scripts designed to be run in the following order: 1. index_pipeline.nf, 2. mapping_pipeline.nf, 3. intersect_pipeline.nf. Below you will find a brief explanation of the scripts and their dependencies:


- **index_pipeline.nf:** This script retrieves the chromosome.fasta file stored in the data directory and indexes it using Burrows-Wheeler Aligner (BWA). Please            install and activate BWA in the current environment before running this script with:

```
nextflow run index_pipeline.nf
```

> [!NOTE]  
> This script uses relative file paths so no internal modification is needed
<br>

- **mapping_pipeline.nf:** This script is designed to map the sample.fastq.gz files to the indexed chromosome with bwa mem and convert these into .bam files with           samtools. The resulting .bam files will then be converted to .bed files using bedtools and the reads will be trimmed to include just the break sites by               accessing the trim_reads.py script in this repo. Please install and activate BWA, samtools, Python3 and bedtools before running this script with:

```
nextflow run mapping_pipeline.nf
```


> [!WARNING]
> This version of mapping_pipeline.nf does not use relative file paths due to incompatibilities with external programs. Please modify the 'mapping_to_bam' parameters and trim_reads script (3 lines total) BEFORE running the script.
<br>

- **intersect_pipeline.nf:** This script uses bedtools to find the intersections between the trimmed .bed files and the AsiSI cut sites in the data folder. Please install and activate bedtools before running this script with:

```
nextflow run intersect_pipeline.nf
```


> [!WARNING]
> This version of intersect_pipeline.nf does not use relative file paths due to incompatibilities with external programs. Please modify the 'find_intersections' parameters (2 lines total) BEFORE running the script.
<br>

### 1.2 Python analysis

This repository also contains a Python file to analyse the results generated by the Nextflow pipeline. This file requires the packages 'pandas' and 'matplotlib' to be installed in the activated environment. The script can be run with:

```
python3 path/to/script path/to/trimmed/reads path/to/intersected/reads path/to/AsiSI/cut/sites
```

> [!NOTE]  
> This script uses relative file paths.
> The Python analysis script requires 4 arguments in the specified order and will raise a SystemExit exception if >4 or <4 arguments are given
<br>

## 2. Interpretation of results

![scatter_plot](https://github.com/CPalmer3200/bsb_test_repo/assets/145576128/791a0918-af95-4fbe-913d-53908646d959)

Figure 1: *Normalised DNA double-strand break count of Samples 1-16 at AsiSI cut sites. Samples 1, 2, 4, 5, 6, 7 and 8 are likely to be untreated while samples 9-16 are probably treated samples. The group that sample 3 belongs to is inconclusive.*
<br>

<br>

Table 1: *The output of the results_analysis.py script. This script calculates the total breaks, AsiSI cut site percentage and normalised break count for each sample.*
<br>
![image](https://github.com/CPalmer3200/bsb_test_repo/assets/145576128/94fc0c1b-2f75-4768-9232-986ab6f4ea23)


### 2.1 Control and treated samples

The results from the Python data analysis strongly indicate that samples 1, 2, 4, 5, 6, 7 and 8 were untreated DIvA control cells as no evidence of DSBs at AsiSI cut sites were found. Samples 9-16 were likely DIvA cells that have been treated with 4-hydroxytamoxifen (4OHT) as they all exhibited DSBs at AsiSI cut sites.
<br>

### 2.2 The identity of sample 3 is uncertain

I am uncertain of sample 3 as a single DSB was found at an AsiSI cut site - however, this DSB may have already been present in the sample and might not have been induced by treatment with 4OHT.
<br>

### 2.3 Samples 9 and 15 had the greatest percentage of possible AsiSI cut sites
Of the 71 possible AsiSI break sites, samples 9 and 15 had the greatest percentage (5.63%, n=4) of AsiSI cut sites. This was calculated as part of the results_analysis.py script.


## 3. Mapping of questions to answers

*This repository was created in response to the Broken String Biosciences Bioinformatician pre-interview coding test:*

- Step 1 is covered by section 1.1 of this repository

- Figure 1 illustrates the plot that step 2.4 instructs to be generated

- Table 1 shows the output of the python data analysis that is requested in step 2

- Question 1 is answered in section 2.1 of this repository

- Questions 2 and 3 are answered in section 2.2 of this repository

- Question 4 is answered in section 2.3 of this repository


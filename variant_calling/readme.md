# Synopsis
This pipeline aligns **paired-end** FASTQ reads against a reference genome using **BWA** then calls variants using **GATK**, **FreeBayes** and **BCFtools**.

This pipeline consists of 3 main steps:

1. Read mapping and variant calling
2. SNP analysis
3. CNV analysis

> **Notes:** SNP and CNV analyses are based on mutant/WT comparison.

# Organisation
This directory contains 3 sub-directories corresponding to 3 main steps of this pipeline, and a readme.md file.
```bash
variant_calling/
├── CNV_analysis
│   ├── add_gene_information.py
│   ├── compute_log2_ratio.py
│   ├── count_read_by_window.py
│   └── plot_log2_ratio.py
├── read_mapping_variant_calling
│   ├── combine_vcf_columns.py
│   ├── GenerateScripts_AlignAndCall.sh
│   └── run_bwa_gatk_freebayes_bcftools.template.sh
├── readme.md
└── SNP_analysis
    ├── add_allele_frequency.py
    ├── add_gene_description.py
    ├── aggregate_variant.py
    ├── build_presence_absence_matrix.py
    ├── calculate_allele_frequency.py
    ├── combine_variant.py
    ├── compare_variant.py
    ├── remove_WT_variant.py
    ├── remove_WT_variant.sh
    └── run_snp_analysis.sh
```

# Preparation
## 1. Prepare your reference genome
* Sequence file: fasta format 
* Annotation file: gff format
## 2. Prepare your FASTQ directory
The FASTQ directory should be organized so that each sample has its own folder containing the corresponding R1 and R2 read files. 

R1 read file must be named as `SampleName_R1_paired.fastq.gz` and R2 read file must be named as `SampleName_R2_paired.fastq.gz`, where `SampleName` matches the name of the sample.

> **Notes:** The FASTQ directory structure and FASTQ file naming pattern can be modified by updating the `R1_path` and `R2_path` variables in the script `run_bwa_gatk_freebayes_bcftools.template.sh`

# Pipeline
## I. Read mapping and variant calling
### Usage
#### 0. Preparation
* Create your working directory for this step and navigate into it
```bash
mkdir read_mapping
cd read_mapping
```

* Copy the scripts from the `read_mapping_variant_calling` subdirectory into your working directory
```bash
cp /project/def-mouellet/Scripts_MOU/PNP/alliancecan/variant_calling/read_mapping_variant_calling/* .
```

* If it is your first time using GATK, run this command on a login node:
```bash
module load StdEnv/2023 gatk/4.4.0.0
```
Then you will see this message:
```
gatk/4.4.0.0:
============================================================================================
Using this software requires you to accept a license on the software website.
Please confirm that you registered on the website below (yes/no).

Utiliser ce logiciel nécessite que vous acceptiez une licence sur le site de l'auteur.
Veuillez confirmer que vous vous êtes enregistrés sur le site web ci-dessous (oui/non).
============================================================================================

https://software.broadinstitute.org/gatk/download/licensing.php
```
Answer yes to continue.

#### 1. Create a sample list file
This file must list the name of the samples to analyse, one sample per line.

#### 2. Generate scripts
To display usage instructions for the script:
```bash
bash GenerateScripts_AlignAndCall.sh
```
Then you will see:
```
Usage: bash GenerateScripts_AlignAndCall.sh <template.sh> <sample_list> <scheduler> <ref_fasta_path> <fastq_dir_path>
```
In which:
* `template.sh`: the template script, in this case `run_bwa_gatk_freebayes_bcftools.template.sh`
* `sample_list`: your sample list file (e.g. `sample_list.txt`)
* `scheduler`: sbatch
* `ref_fasta_path`: path to your reference fasta file
* `fastq_dir_path`: path to your FASTQ directory

Example: 
```bash
bash GenerateScripts_AlignAndCall.sh run_bwa_gatk_freebayes_bcftools.template.sh sample_list.txt sbatch ../ref/TriTrypDB-68_LinfantumJPCM5_Genome.fasta ../fastqFiles/
```
A launch script to submit simultaneously all the scripts generated from the template script will be created. This script is named after the template script with a "-Launch.sh" suffix. In this case, it will be `run_bwa_gatk_freebayes_bcftools.template.sh-Launch.sh`.

#### 3. Run the launch script
```bash
bash run_bwa_gatk_freebayes_bcftools_filter.template.sh-Launch.sh
```

### Output
Output files for each sample will be saved in a folder named after the sample.

Extension | Description
----------|--------------
*_sorted_dedup_reads.bam | Alignment in BAM format. BAM file is sorted by coordinate and duplicates are marked
*_sorted_dedup_reads.bam.bai | Index file accompanying the BAM file
*_sorted_dedup_reads.bam.sbi | Bam splitting index file used by Spark to split the job into n-parts running in parallel
*_alignment_metrics.txt | Summary statistics that describe the alignment of sequencing reads in the BAM file
*_insert_size_metrics.txt | Statistics about the insert sizes (aka fragment lengths) of paired-end reads
*_insert_size_histogram.pdf | Histogram of the insert size distribution
*_raw_variants.gatk.vcf | Raw variants called by GATK HaplotypeCaller
*_raw_variants.gatk.vcf.idx | Index file accompanying the vcf file
*_raw_variants.freebayes.vcf | Raw variants called by FreeBayes
*_raw_variants.bcftools.vcf | Raw variants called by Bcftools mpileup
*_filtered_variants.gatk.vcf | GATK HaplotypeCaller variants filtered according to https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
*_filtered_variants.freebayes.vcf | FreeBayes variants filtered by quality score (QUAL>=20, i.e. 99% confidence or 1% error probability) and depth (DP>=10)
*_filtered_variants.bcftools.vcf | Bcftools mpileup variants filtered by quality score (QUAL>=20, i.e. 99% confidence or 1% error probability) and depth (DP>=10)

> **Notes:** For debugging purpose, intermediate files are placed in your scratch in a folder named `SampleName_bwa_gatk_freebayes_bcftools`, where SampleName is the name of the sample. If the script completes successfully, you can safely delete these files.



## II. SNP analysis
### Usage
#### 0. Preparation
* Create your working directory for this step and navigate into it
```bash
mkdir SNP_analysis
cd SNP_analysis
```
* Copy the scripts from the `SNP_analysis` subdirectory to your working directory
```bash
cp /project/def-mouellet/Scripts_MOU/PNP/alliancecan/variant_calling/SNP_analysis/* .
```
* Copy or create symbolic links of the snpEff container `snpeff_v5.2.sif` to your working directory
```bash
ln -s /project/def-mouellet/Scripts_MOU/PNP/containers/snpeff_v5.2.sif
```

#### 1. Remove WT variants & N-stretch-containing variants
For **each lineage** of mutants:
##### a. Create a lineage mutant list file
This file must list the name of the mutants of this lineage, one sample per line.

##### b. Run the script `remove_WT_variant.sh`

```bash
bash remove_WT_variant.sh <calling_folder> <wt_sample_name> <mutant_list>
```
* calling_folder: path to directory where you performed the first step "Read mapping and variant calling" for the WT and its mutants

#### 2. Run SNP analysis
##### a. Create a mutant list file
This file must list the name of the mutants to analyse, one sample per line.

##### b. Run the script run_snp_analysis.sh
```bash
bash run_snp_analysis.sh <path_to_snpEff_container> <genome_name> <path_to_ref_fasta> <path_to_ref_gff> <mutant_list>
```
E.g.
```bash
bash run_snp_analysis.sh snpeff_v5.2.sif TriTrypDB-68_LinfantumJPCM5 ../ref/TriTrypDB-68_LinfantumJPCM5_Genome.fasta ../ref/TriTrypDB-68_LinfantumJPCM5.gff mutant_list.txt
```

### Output
Two main output files are: all_variant_nosyn.tsv and all_variant_nosyn_matrix.tsv
#### 1. all_variant_nosyn.tsv
This file list all filtered, non synonymous variants detected by one of the three algorithms: gatk (GATK HaplotypeCaller), freebayes (FreeBayes) and bcftools (Bcftools mpileup). The columns are:
* #CHROM
* POS
* GENOTYPE
* ALLELE_FREQ
* GENE_ID
* GENE_NAME
* DESCRIPTION
* FEATURE_ID
* VARIANT_TYPE
* NU_CHANGE
* AA_CHANGE
* CDS_POS/CDS_LENGTH
* AA_POS/AA_LENGTH
* SAMPLE
* CALLED_BY

#### 2. all_variant_nosyn_matrix.tsv
This file is a matrix of presence/absence of SNP in genes. The columns are:
* #CHROM
* GENE_ID
* GENOTYPE
* GENE_NAME
* DESCRIPTION
* VARIANT_TYPE
* CALLED_BY

Samples start from the 8th column, one sample per column.

## III. CNV analysis 



# Dependencies
## Modules
* gcc/12.3
* bwa/0.7.18
* samtools/1.20
* gatk/4.4.0.0
* freebayes/1.3.7
* bcftools/1.19
* r/4.4.0
## Containers
* snpeff_v5.2.sif
* python_3.11.5.sif

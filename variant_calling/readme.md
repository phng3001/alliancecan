# Synopsis
This pipeline aligns paired-end FASTQ reads against a reference genome using BWA then calls variants using GATK, FreeBayes and BCFtools.

# Pipeline
## Read mapping and variant calling
### Usage
#### 0. Preparation
* Copy the scripts in the **read_mapping_variant_calling** directory to your working directory
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
This file must list the name of the samples to analyse, one sample per line

#### 2. Generate scripts
```
bash GenerateScripts_AlignAndCall.sh <template.sh> <sample_list> <scheduler> <ref_fasta_path> <fastq_dir_path>
```
* template script: run_bwa_gatk_freebayes_bcftools.template.sh
* scheduler: sbatch
* ref_fasta_path: path to the reference fasta file
* fastq_dir_path: path to directory containing FASTQ reads



## SNP analysis



## CNV analysis 



# Dependencies
* gcc/12.3
* bwa/0.7.18
* samtools/1.20
* gatk/4.4.0.0
* freebayes/1.3.7
* bcftools/1.19
* r/4.4.0



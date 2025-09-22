# Synopsis
This pipeline aligns **paired-end** FASTQ reads against a reference genome using **BWA** then calls variants using **GATK**, **FreeBayes** and **BCFtools**.

This pipeline consists of 3 main parts:

1. Read mapping and variant calling
2. SNP analysis
3. CNV analysis

> **Notes:** SNP and CNV analyses are based on mutant/WT comparison.

# Organisation
This directory contains:
* 3 sub-directories corresponding to 3 main parts of this pipeline: `read_mapping_variant_calling`, `snp_analysis` and `cnv_analysis`.
* A container subdirectory `containers` containing the necessary containers to run this pipeline.
* A `readme.md` file.
```
variant_calling/
├── cnv_analysis
│   ├── add_gene_information.py
│   ├── compute_log2_ratio.py
│   ├── count_read_by_window.py
│   └── plot_log2_ratio.py
├── containers
│   ├── python_3.11.5.def
│   ├── python_3.11.5.sha256
│   ├── python_3.11.5.sif
│   ├── snpeff_v5.2.sha256
│   └── snpeff_v5.2.sif
├── read_mapping_variant_calling
│   ├── combine_vcf_columns.py
│   ├── GenerateScripts_AlignAndCall.sh
│   └── run_bwa_gatk_freebayes_bcftools.template.sh
├── readme.md
└── snp_analysis
    ├── ncbi
    │   ├── add_allele_frequency.py
    │   ├── add_gene_description.py
    │   ├── aggregate_variant.py
    │   ├── build_presence_absence_matrix.py
    │   ├── calculate_allele_frequency.py
    │   ├── combine_variant.py
    │   ├── compare_variant.py
    │   ├── filter_out_tsv.py
    │   ├── remove_WT_variant.py
    │   ├── remove_WT_variant.sh
    │   └── run_snp_analysis.sh
    └── tritrypdb
        ├── add_allele_frequency.py
        ├── add_gene_description.py
        ├── aggregate_variant.py
        ├── build_presence_absence_matrix.py
        ├── calculate_allele_frequency.py
        ├── combine_variant.py
        ├── compare_variant.py
        ├── filter_out_tsv.py
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
## Part I. Read mapping and variant calling
### Protocol
#### Step 0. Preparation
##### a. Create your working directory for this part and navigate into it
```bash
mkdir read_mapping && cd read_mapping
```

##### b. Copy the scripts from the subdirectory `read_mapping_variant_calling` into your working directory
```bash
cp /project/def-mouellet/Scripts_MOU/PNP/alliancecan/variant_calling/read_mapping_variant_calling/* .
```

##### c. If it is your first time using GATK, run this command on a login node:
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

#### Step 1. Create a sample list file
This file must list the name of the samples to analyse, one sample per line.

#### Step 2. Generate scripts
To display usage instructions for the script `GenerateScripts_AlignAndCall.sh`:
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
A launch script to submit simultaneously all the scripts generated from the template script will be created. This script is named after the template script with a "-Launch.sh" suffix. In this case, the launch script is `run_bwa_gatk_freebayes_bcftools.template.sh-Launch.sh`.

#### Step 3. Run the launch script
```bash
bash run_bwa_gatk_freebayes_bcftools.template.sh-Launch.sh
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



## Part II. SNP analysis
### Step 0. Preparation
#### a. Put all output folders from Part I - `Read mapping and variant calling` into a single directory

#### b. Create your working directory for this part and navigate into it
```bash
mkdir SNP_analysis && cd SNP_analysis
```

#### c. Copy the appropriate scripts from the subdirectory `snp_analysis` to your working directory
* If your reference genome annotation gff file is from TriTrypDB:
```bash
cp /project/def-mouellet/Scripts_MOU/PNP/alliancecan/variant_calling/snp_analysis/tritrypdb/* .
```
* If your reference genome annotation gff file is from NCBI or Prokka: 
```bash
cp /project/def-mouellet/Scripts_MOU/PNP/alliancecan/variant_calling/snp_analysis/ncbi/* .
```

#### d. Copy or create symbolic links of the snpEff container `snpeff_v5.2.sif` to your working directory
```bash
ln -s /project/def-mouellet/Scripts_MOU/PNP/alliancecan/variant_calling/containers/snpeff_v5.2.sif
```

### Step 1. Remove WT variants & N-stretch-containing variants
#### a. Create mutant list file(s)
For **each lineage** of mutants, create a mutant list file. This file must list the name of the mutants of this lineage, one sample per line.

#### b. Run the script `remove_WT_variant.sh`
##### Usage
To display usage instructions for this script:
```bash
bash remove_WT_variant.sh
```
Then you will see:
```
Usage: bash remove_WT_variant.sh <calling_folder> <wt_sample_name> <mutant_list>
```
In which: 
* `calling_folder`: path to directory containing the outputs of Part I - `Read mapping and variant calling` of the WT and its mutants
* `wt_sample_name`: the name of the WT of this lineage
* `mutant_list`: the lineage mutant list file

Example:
```bash
bash remove_WT_variant.sh ../read_mapping/ ldi263WT_cl1_FF mutant_list.txt
```

##### Output
For each sample: 
Extension | Description
----------|--------------
*_filtered_variants.noWT.gatk.vcf | Filtered GATK HaplotypeCaller variants, excluding variants found in WT
*_filtered_variants.noWT.freebayes.vcf | Filtered FreeBayes variants, excluding variants found in WT
*_filtered_variants.noWT.bcftools.vcf | Filtered Bcftools mpileup variants, excluding variants found in WT

### Step 2. Run SNP analysis
#### a. Create a mutant list file
This file must list the name of the mutants to analyse, one sample per line.

#### b. Run the script `run_snp_analysis.sh`
##### Usage
To display usage instructions for this script:
```bash
bash run_snp_analysis.sh
```
Then you will see:
```
Usage: bash run_snp_analysis.sh <snpEff_container_path> <genome_name> <ref_fasta_path> <ref_gff_path> <mutant_list>
```
In which:
* `snpEff_container_path`: path to the snpEff container, in this case `snpeff_v5.2.sif`
* `genome_name`: name of your reference genome
* `ref_fasta_path`: path to your reference fasta file
* `ref_gff_path`: path to your reference gff file
* `mutant_list`: your mutant list file

Example
```bash
bash run_snp_analysis.sh snpeff_v5.2.sif TriTrypDB-68_LinfantumJPCM5 ../ref/TriTrypDB-68_LinfantumJPCM5_Genome.fasta ../ref/TriTrypDB-68_LinfantumJPCM5.gff mutant_list.txt
```

##### Output
###### Main outputs
Two main output files are `all_variant.tsv` and `all_variant_matrix.tsv`

1. `all_variant.tsv`

This file lists all filtered variants detected by one of the three algorithms: gatk (GATK HaplotypeCaller), freebayes (FreeBayes) and bcftools (Bcftools mpileup). The columns are:
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

2. `all_variant_matrix.tsv`

This file is a matrix of presence/absence of SNP in genes. The columns are:
* #CHROM
* GENE_ID
* GENOTYPE
* GENE_NAME
* DESCRIPTION
* VARIANT_TYPE
* CALLED_BY

Samples start from the 8th column, one sample per column.

###### Other outputs
* For each sample:

Extension | Description
----------|--------------
*_filtered_variants.noWT.snpEff.gatk.vcf | Filtered GATK HaplotypeCaller variants, excluding variants found in WT, annotated by SnpEff
*_filtered_variants.noWT.snpEff.freebayes.vcf | Filtered FreeBayes variants, excluding variants found in WT, annotated by SnpEff
*_filtered_variants.noWT.snpEff.bcftools.vcf | Filtered Bcftools mpileup variants, excluding variants found in WT, annotated by SnpEff

* `all_variant.gatk.tsv`: lists all filtered variants detected by gatk
* `all_variant.freebayes.tsv`: lists all filtered variants detected by freebayes
* `all_variant.bcftools.tsv`: lists all filtered variants detected by bcftools

### Step 3. Filter out variants with the script `filter_out_tsv.py` (optional)
In this example we will use the script `filter_out_tsv.py` to filter out some variant types from the output variant list `all_variant.tsv`.

To display usage instructions for this script:
```python
python filter_out_tsv.py -h
```
Then you will see:
```
usage: filter_out_tsv.py [-h] -i INPUT -c COLUMN -f FILTER -o OUTPUT

Remove rows from a TSV file where the specified column matches a given value.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input TSV file
  -c COLUMN, --column COLUMN
                        Column name to filter on
  -f FILTER, --filter FILTER
                        Value to filter out (i.e., exclude rows with this value)
  -o OUTPUT, --output OUTPUT
                        Output TSV file
```
The column for filtering in this example is `VARIANT_TYPE`. Pass this to the `-c` argument of the script.

#### Filter out synonymous variants
SnpEff annotates synonymous variants with the value `synonymous_variant`. Pass this to the `-f` argument of the script.

```python
python filter_out_tsv.py \
  -i all_variant.tsv \
  -c VARIANT_TYPE \
  -f synonymous_variant \
  -o all_variant_nosyn.tsv
```

#### Filter out intergenic variants
SnpEff annotates intergenic variants with the value `intergenic_region`. Pass this to the `-f` argument of the script.

```python
python filter_out_tsv.py \
  -i all_variant_nosyn.tsv \
  -c VARIANT_TYPE \
  -f intergenic_region \
  -o all_variant_nosyn_noig.tsv
```



## Part III. CNV analysis 
### Step 0. Preparation
#### a. Create your working directory for this part and navigate into it
```bash
mkdir CNV_analysis && cd CNV_analysis
```

#### b. Copy the scripts from the subdirectory `cnv_analysis` to your working directory
```bash
cp /project/def-mouellet/Scripts_MOU/PNP/alliancecan/variant_calling/cnv_analysis/* .
```

#### c. Copy or create symbolic links of the python container `python_3.11.5.sif` to your working directory
```bash
ln -s /project/def-mouellet/Scripts_MOU/PNP/alliancecan/variant_calling/containers/python_3.11.5.sif
```

#### d. Load Apptainer
```bash
module load apptainer
```

### Step 1. Create a sample list file
This file must list the name of the samples to analyse, one sample per line.

### Step 2. Make symbolic links of BAM files to your working directory
>**Notes:** BAM files were generated at Part I of this pipeline "Read mapping and variant calling"

Example:
```bash
for X in $(cat sample_list.txt); 
do 
  ln -s ../read_mapping/$X/${X}_sorted_dedup_reads.bam ${X}.bam; 
done
```

### Step 3. Copy or make symbolic links of the reference fasta and gff files to your working directory

### Step 4. Count mapped reads by genomic windows for each sample
#### a. Load Samtools
```bash
module load samtools
```
#### b. Run the script `count_read_by_window.py`
##### Usage
To display usage instructions for this script:
```python
python count_read_by_window.py -h
```
Then you will see:
```
usage: count_read_by_window.py [-h] -i BAM_FILE -r FASTA_FILE -w WINDOW_SIZE -o OUTPUT_FILE

Count reads in BAM file by genomic windows.

options:
  -h, --help            show this help message and exit
  -i BAM_FILE, --bam_file BAM_FILE
                        Input BAM file path
  -r FASTA_FILE, --fasta_file FASTA_FILE
                        Reference FASTA file path
  -w WINDOW_SIZE, --window_size WINDOW_SIZE
                        Size of each genomic window (e.g., 5000).
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output TSV file path
```
Example:
```bash
for X in $(cat sample_list.txt); 
do 
  python count_read_by_window.py \
  -i ${X}.bam \
  -r TriTrypDB-68_LinfantumJPCM5_Genome.fasta \
  -w 5000 \
  -o ${X}_cov.tsv;
done
```

##### Output
This script outputs TSV file showing the number of reads mapped to each genomic window. The columns are:
* #CHROM
* START
* END
* NB_READS

### Step 5. Compute the normalized log2 mutant/WT read ratio
#### a. Create mutant list file(s)
For **each lineage** of mutants, create a mutant list file. This file must list the name of the mutants of this lineage, one sample per line.

#### b. Run the script `compute_log2_ratio.py`
##### Usage
To display usage instructions for this script:
```python
python compute_log2_ratio.py -h
```
Then you will see:
```
usage: compute_log2_ratio.py [-h] -i MUTANT_FILE -r WT_FILE -o OUTPUT_FILE

Compute the mutant/wild type read ratio in log2 for each corresponding genomic window. Read counts are normalized by the total read count per sample.

options:
  -h, --help            show this help message and exit
  -i MUTANT_FILE, --mutant_file MUTANT_FILE
                        Mutant read counts file path
  -r WT_FILE, --wt_file WT_FILE
                        Wild type read counts file path
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output file path
```
The `MUTANT_FILE` and `WT_FILE` are outputs generated by the script `count_read_by_window.py` in step 4. Pass these to the `-i` and `-r` arguments of the script, respectively.

Example:
```bash
for X in $(cat mutant_list.txt);
do
  apptainer run python_3.11.5.sif \
  python compute_log2_ratio.py \
  -i ${X}_cov.tsv \
  -r ldi263WT_cl1_FF_cov.tsv \
  -o ${X}_vs_ldi263WT_cl1_FF.tsv;
  echo "$X done";
done
```

##### Output
This script outputs a TSV file showing the normalized log2 mutant/WT read ratio for each genomic window. The columns are:
* #CHROM
* START
* END
* NB_READS_MUTANT
* NB_READS_WT
* NB_READS_RATIO_MUTANT/WT
* TOTAL_READS_MUTANT
* TOTAL_READS_WT
* TOTAL_READS_RATIO_MUTANT/WT
* NORMALIZED_LOG2_READS_RATIO_MUTANT/WT

### Step 6. Add gene information with the script `add_gene_information`
#### Usage
To display usage instructions for this script:
```python
python add_gene_information.py -h
```
Then you will see:
```
usage: add_gene_information.py [-h] -i INPUT_FILE -r GFF_FILE [-f FEATURE] -o OUTPUT_FILE

Annotate genomic regions with gene information retrieved from a GFF file. Input file must contain these columns: '#CHROM', 'START', and 'END'. Gene information will be written to a new column named 'GENE'.

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file INPUT_FILE
                        Input TSV file path
  -r GFF_FILE, --gff_file GFF_FILE
                        Input GFF file path
  -f FEATURE, --feature FEATURE
                        Feature type to annotate (e.g., 'gene', 'CDS'). Default is 'gene'
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output TSV file path
```
In which:
* `INPUT_FILE`: Output generated by the script `compute_log2_ratio.py` in step 5. Pass this to the `-i` argument of the script.
* `GFF_FILE`: The reference annotation gff file. Pass this to the `-r` argument of the script.
* `FEATURE`: Feature type in gff file to be used for annotation. Pass this to the `-f` argument of the script.
>**Notes:**  Check your GFF file to select the appropriate feature type in the 3rd column to pass. The script will extract from the 9th column the following tags: `ID` or `locus_tag`, `Name` or `gene`, `product` or `description`.

Example:
```bash
for X in $(cat mutant_list.txt);
do
  apptainer run python_3.11.5.sif \
  python add_gene_information.py \
  -i ${X}_vs_ldi263WT_cl1_FF.tsv \
  -r TriTrypDB-68_LinfantumJPCM5.gff \
  -f protein_coding_gene \
  -o ${X}_vs_ldi263WT_cl1_FF.annot.tsv;
  echo "$X done";
done
```

#### Output
This script adds a new column named `GENE` containing information about gene(s) found in the corresponding genomic window in the following format: 
* Unique gene: (gene_locus_tag,gene_name,gene_description)
* Multiple genes: (gene1_locus_tag,gene1_name,gene1_description); (gene2_locus_tag,gene2_name,gene2_description)

### Step 7. Make plots
#### a. Create a list of files to plot
These files are outputs generated by the script `add_gene_information` in step 6. List the names of the files that you want to plot, one per line.

#### b. Run the script `plot_log2_ratio.py`
##### Usage
To display usage instructions for this script:
```python
python plot_log2_ratio.py -h
```
Then you will see:
```
usage: plot_log2_ratio.py [-h] -i INPUT [INPUT ...] -o OUTDIR [--ylim YLIM YLIM]
                          [--gridstep GRIDSTEP] [--dpi DPI] [--width WIDTH] [--height HEIGHT]

Generate interactive plots of the normalized log2 mutant/wild type read ratio for each chromosome.

options:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Input TSV files (one per sample)
  -o OUTDIR, --outdir OUTDIR
                        Output directory to save plots
  --ylim YLIM YLIM      y-axis limits as min max (e.g., --ylim -2 2)
  --gridstep GRIDSTEP   Step size for y-axis gridlines. Default is 1
  --dpi DPI             Plot resolution in dpi (dots per inch). Default is 96
  --width WIDTH         Figure width in inches. Default is 12
  --height HEIGHT       Figure height in inches. Default is 3
```

Example:
```bash
apptainer run python_3.11.5.sif \
python plot_log2_ratio.py \
-i $(paste -s -d ' ' mutant_file.txt) \
-o interactive_plots \
--ylim -4 4 \
--width 16 \
--height 4
```

##### Output
This script generates an interactive plot for each chromosome, showing the normalized log2 ratio of mutant to WT read coverage across genomic coordinates.



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

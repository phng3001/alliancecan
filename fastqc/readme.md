# Usage
## Multiple samples
### 1. Check the FASTQ directory structure and FASTQ file naming pattern 
Modify R1_path and R2_path variables in the template script run_fastqc.template.sh if necessary.

### 2. Create a directory to save fastqc output files 
```bash
mkdir fastqc_output
```
### 3. Generate scripts
```bash
bash GenerateScripts_fastqc.sh <template.sh> <sample_list> <scheduler> <fastq_dir_path> <output_dir_path>
```
* `template.sh`: the template script `run_fastqc.template.sh`
* `sample_list`: list of sample, one per line
* `scheduler`: sbatch
* `fastq_dir_path`: path to directory containing fastq reads (e.g. fastq in the example below)
```
fastq/
├── ERR065289
│   ├── ERR065289_R1_paired.fastq.gz
│   └── ERR065289_R2_paired.fastq.gz
└── ERR065290
    ├── ERR065290_R1_paired.fastq.gz
    └── ERR065290_R2_paired.fastq.gz
```
* `output_dir_path`: path to the directory to save fastqc output files (e.g. `fastqc_output` in the example at the 2nd step)

### 4. Launch scripts
```bash
bash run_fastqc.template.sh-Launch.sh
```

# More info
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

https://github.com/s-andrews/FastQC
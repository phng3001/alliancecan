# Usage
## Multiple samples
### 1. Create a directory to save fastqc output files 
```bash
mkdir fastqc_output
```
### 2. Generate scripts
```bash
bash GenerateScripts_fastqc.sh <template.sh> <sample_list> <scheduler> <fastq_dir_path> <output_dir_path>
```
* template script: run_fastqc.template.sh
* fastq_dir_path: path to directory containing fastq reads (e.g. fastq in the example below)
```
fastq/
├── ERR065289
│   ├── ERR065289_R1_paired.fastq.gz
│   └── ERR065289_R2_paired.fastq.gz
└── ERR065290
    ├── ERR065290_R1_paired.fastq.gz
    └── ERR065290_R2_paired.fastq.gz
```
* output_dir_path: path to the directory to save fastqc output files (e.g. fastqc_output in the example at the 1st step)

### 3. Launch scripts
```bash
bash run_fastqc.template.sh-Launch.sh
```

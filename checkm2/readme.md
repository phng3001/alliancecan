# Source
https://github.com/chklovski/CheckM2

# Function
Rapid assessment of genome bin quality (completeness and contamination) using machine learning

# Container
checkm2_v1.1.0.sif

# Script
## `run_apptainer_checkm2.sh`
> Notes: Check all genomes at once

### Usage
#### 1. Check if the database path set in the variable `CHECKM2DB` is valid

#### 2. Run the script
```bash
bash run_apptainer_checkm2.sh <container> <fasta_dir_path> <fasta_extension> <output_dir_name>
```
or submit a slurm job:
```bash
sbatch run_apptainer_checkm2.sh <container> <fasta_dir_path> <fasta_extension> <output_dir_name>
```

* `container`: the CheckM2 container `checkm2_v1.1.0.sif`
* `fasta_dir_path`: path to directory containing genome fasta files
* `fasta_extension`: extension of the genome fasta files
* `output_dir_name`: output directory name

### Output
The main CheckM2 output written in the output directory is a TSV file named `quality_report.tsv` 

## `run_apptainer_checkm2.template.sh`
> Notes: Split the list of genomes into batches and check batches independently

### Usage
#### 1. Check if the database path set in the variable `CHECKM2DB` in the template script is valid

#### 2. Generate scripts
```bash
bash GenerateScriptsBySubsamples_apptainer_checkm2.sh <template.sh> <sample_list> <batch_size> <scheduler> <container> <fasta_dir_path> <fasta_extension>
```
* `template.sh`: the template script `run_apptainer_checkm2.template.sh`
* `sample_list`: list of sample, one per line
* `batch_size`: the sample list will be splitted into sub-list of this size
* `scheduler`: sbatch
* `container`: the CheckM2 container, in this case `checkm2_v1.1.0.sif`
* `fasta_dir_path`: path to directory containing genome fasta files
* `fasta_extension`: extension of the genome fasta files

#### 3. Launch scripts
```bash
bash run_apptainer_checkm2.template.sh-Launch.sh
```

### Output
#### The input sample list file (e.g. `sample_list.txt`) will be splitted into sub-lists (batches) named after its name with `batch_*` suffixes (e.g. `sample_list.txt_batch_000`, `sample_list.txt_batch_001` etc)
#### For each batch: 
* The input genome fasta files will be copied to a directory named after the corresponding sub-list file with a `_dir` suffix (e.g. `sample_list.txt_batch_000_dir`, `sample_list.txt_batch_001_dir` etc)
* The CheckM2 output files will be written to a directory named after the corresponding sub-list file with a `_checkm2_result` suffix (e.g. `sample_list.txt_batch_000_checkm2_result`, `sample_list.txt_batch_001_checkm2_result` etc)

### Postprocessing
Example
```bash
# Get the header
head -n 1 sample_list.txt_batch_000_checkm2_result/quality_report.tsv > quality_report.tsv
# Merge CheckM2 outputs
tail -n +2 -q *_checkm2_result/quality_report.tsv >> quality_report.tsv
```

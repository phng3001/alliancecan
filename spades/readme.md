# `run_spades_auto.template.sh`
> Notes: k-mer sizes are automatically chosen by SPAdes

## Usage
### 1. Set variables in the template script
#### a. Check the FASTQ directory structure and FASTQ file naming pattern, modify the `R1_path` and `R2_path` variables if necessary
#### b. Modify the variable `min_contig_lgth` if necessary to change the minimum contig length to keep

### 2. Generate scripts
```bash
bash GenerateScripts_spades.sh <template.sh> <sample_list> <scheduler> <fastq_dir_path>
```
* `template.sh`: the template script
* `sample_list`: list of sample, one per line 
* `scheduler`: sbatch
* `fastq_dir_path`: path to directory containing fastq reads

### 3. Launch scripts
```bash
bash run_spades_auto.template.sh-Launch.sh
``` 

# `run_spades.template.sh`

## Usage
Similar to `run_spades_auto.template.sh` but k-mer sizes could be modified in the template script


# More info
https://github.com/ablab/spades
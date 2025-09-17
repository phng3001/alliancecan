# Usage
## Multiple samples
### 1. Check the FASTQ directory structure and FASTQ file naming pattern
Modify `R1_path` and `R2_path` variables in the template script `run_trimmomatic.template.sh` if necessary.

### 2. Generate scripts
```bash
bash GenerateScripts_trimmomatic.sh <template.sh> <sample_list> <scheduler> <fastq_dir_path> <adapter_file_path>
```
* `template.sh`: the template script `run_trimmomatic.template.sh`
* `sample_list`: list of sample, one per line 
* `scheduler`: sbatch
* `fastq_dir_path`: path to directory containing fastq reads
* `adapter_file_path`: path to the adapter fasta file (e.g. `TruSeq3-PE-2.fa`)

### 3. Launch scripts
```bash
bash run_trimmomatic.template.sh-Launch.sh
``` 

# More info
https://github.com/usadellab/Trimmomatic
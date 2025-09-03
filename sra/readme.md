# Usage
## Download SRA reads
### 1 sample
```bash
bash run_download_sra.sh <run_id>
```
#### Example
```bash
bash run_download_sra.sh ERR065289
```
### Multiple samples
#### Generate scripts
```bash
bash GenerateScripts.sh <template.sh> <sample_list> <scheduler>
```
* template script: run_download_sra.template.sh
#### Launch scripts
```bash
bash run_download_sra.template.sh-Launch.sh
```
### Multiple subsamples
#### Generate scripts
```bash
bash GenerateScriptsBySubsamples.sh <template.sh> <sample_list> <batch_size> <scheduler>
```
* template script: run_download_sra.subsample.template.sh
#### Launch scripts
```bash
bash run_download_sra.subsample.template.sh-Launch.sh
```

## Rename SRA reads
Paired reads downloaded from SRA have different names, e.g. "ERR065289.1.1" and "ERR065289.1.2", which could cause errors in subsequent analyses (e.g. BWA-MEM)

This could be fixed by the python script **rename_sra_read.py**, which will remove the 2 last characters in the read names, so that paired reads have the same name, e.g. "ERR065289.1"

### 1 sample
```python
python rename_read.py <sra_prefix> <input_fastq> <output_fastq>
```
#### Example
```python
python rename_read.py ERR ERR065289_1.fastq ERR065289_R1_paired.fastq
python rename_read.py ERR ERR065289_2.fastq ERR065289_R2_paired.fastq
```
### Multiple samples
#### Generate scripts
```bash
bash GenerateScripts_rename_sra_read.sh <template.sh> <sample_list> <scheduler> <sra_prefix> <fastq_dir_path>
```
* **Note:** The script rename_read.py should be present in the same directory
* **template script:** run_rename_sra_read.template.sh
* **fastq_dir_path:** path to directory where SRA raw reads were downloaded (e.g. sra in the example below)
```
sra/
├── ERR065289
│   ├── ERR065289_1.fastq.gz
│   └── ERR065289_2.fastq.gz
└── ERR065290
    ├── ERR065290_1.fastq.gz
    └── ERR065290_2.fastq.gz
```

#### Launch scripts
```bash
bash run_rename_sra_read.template.sh-Launch.sh
```

### Multiple subsamples
#### Generate scripts
```bash
bash GenerateScriptsBySubsamples.sh <template.sh> <sample_list> <batch_size> <scheduler>
```
* template script: run_rename_sra_read.subsample.template.sh
#### Launch scripts
```bash
bash run_rename_sra_read.subsample.template.sh-Launch.sh
```

# Note
> 2025-08-29
* Nibi: worked on login node and compute node
* Fir: worked on login node and compute node
* Rorqual: worked only on login node
* Narval: worked only on login node

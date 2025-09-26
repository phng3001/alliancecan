# Usage
## 1 sample
```bash
# The first two arguments are mandatory
bash run_prokka.sh <sample_name> <fasta_file> [kingdom] [genus] [species] [strain]
```

## Multiple subsamples
### 1. Check the template script `run_prokka.template.sh`, modify the fasta file path pattern if necessary

### 2. Generate scripts
```bash
bash GenerateScriptsBySubsamples_prokka.sh <template.sh> <sample_list> <batch_size> <scheduler> <fasta_dir_path> <fasta_extension> [kingdom] [genus] [species]
```

### 3. Launch scripts
```bash
bash run_prokka.template.sh-Launch.sh
```
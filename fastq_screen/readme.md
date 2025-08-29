# Source
https://github.com/StevenWingett/FastQ-Screen

# Function
FastQ Screen allows us to screen a library of sequences in FastQ format against a set of sequence databases so we can see if the composition of the library matches with what we expect

# Container 
fastq-screen_v0.16.0.sif
> **Note:** This container already has bowtie2 installed inside at /usr/local/bin/bowtie2

# Setup databases
- Build Bowtie2 indexes for the genomes we want to screen
> **Note:** FastQ Screen can process up to a maximum of 32 reference genomes
- Editing the configuration file to reflect the actual set up databases
> **Note:** The current configuration file (fastq_screen.conf) has been set up to scan for Human, Mouse, *E. coli*, PhiX, common adapters, common vectors, *Leishmania* (*L. major*, *L. infantum* and *L. tarentolae*) and *S. pneumoniae* genomes

# Usage
## 1 sample
```bash
bash run_apptainer_fastq_screen.sh <container> <conf_file> <fastq_file>
```
> **Note:** By default sample is subsampled at 100000 reads.

## Multiple samples
### Generate scripts
```bash
bash GenerateScripts_fastq_screen.sh <template.sh> <sample_list> <scheduler> <container> <conf_file> <fastq_dir_path>
```
> **Template script:** run_apptainer_fastq_screen.template.sh
### Launch scripts
```bash
bash run_apptainer_fastq_screen.template.sh-Launch.sh
```

# Note
These scripts were tested on Nibi using *Leishmania* sequencing run 20240506_LEPP048. It took less than 5 minutes per sample

# More info
https://stevenwingett.github.io/FastQ-Screen/

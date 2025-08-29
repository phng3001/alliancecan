# Source
https://github.com/kblin/ncbi-genome-download/

# Function
Download genomes from the NCBI FTP servers

# Container
ncbi-genome-download_v0.3.3.sif

# Usage
## Shell
### Options
```bash
apptainer exec ncbi-genome-download_v0.3.3.sif ncbi-genome-download --help
```
### Example
#### Download all reference genomes of the genus Clostridium
> Check which genomes to download
```bash
apptainer exec \
--bind /etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem:/etc/pki/tls/certs/ca-bundle.crt \
ncbi-genome-download_v0.3.3.sif ncbi-genome-download bacteria \
--genera Clostridium \
--refseq-categories reference \
--dry-run
```
> Download these genomes in the fasta format right into the output folder and save genome metadata
```bash
apptainer exec \
--bind /etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem:/etc/pki/tls/certs/ca-bundle.crt \
ncbi-genome-download_v0.3.3.sif ncbi-genome-download bacteria \
--genera Clostridium \
--refseq-categories reference \
--formats fasta \
--flat-output \
-p 4 \
-o clostridium_genomes \
--metadata-table clostridium_metadata \
--progress-bar
```

## Scripting
```bash
bash run_apptainer_ncbi-genome-download_refseq_bacteria.sh <container> <species_taxid> <format> <output_prefix>
```
or
```bash
bash run_apptainer_ncbi-genome-download_genbank.sh <container> <species_taxid> <format> <output_prefix>
```
### Example
#### Download fasta format of Streptococcus acidominimus assemblies from RefSeq
```bash
bash run_apptainer_ncbi-genome-download_bacteria.sh ncbi-genome-download_v0.3.3.sif 1326 fasta streptococcus_acidominimus
```
#### Download all available formats of Leishmania guyanensis assemblies from GenBank
```bash
bash run_apptainer_ncbi-genome-download_genbank.sh ncbi-genome-download_v0.3.3.sif 5670 all leishmania_guyanensis
```

# GFF to GTF
## Alliance
### 1. Load gffread
> **Notes:** Until now (2027-07-23), gffread is only available on `StdEnv/2020` environment.
```bash
module load StdEnv/2020 gffread/0.12.3
```
### 2. Convert GFF to GTF
```bash
# -F: keep all GFF attributes, -T: output GTF
gffread input.gff -F -T -o output.gtf
```

## Container
gffread_v0.12.7.sif
### 1. Load apptainer
```bash
module load apptainer
```
### 2. Convert GFF to GTF
```bash
apptainer run gffread_v0.12.7.sif gffread input.gff -F -T -o output.gtf
```

# Scripts
## gff_to_tsv_tritrypdb.py
### Usage
```bash
python gff_to_tsv_tritrypdb.py
```
### Example
```bash
python gff_to_tsv_tritrypdb.py \
--gff_file TriTrypDB-68_LinfantumJPCM5.gff \
--output_tsv TriTrypDB-68_LinfantumJPCM5.tsv
```

## gff_to_tsv_ncbi.py
### Usage
```bash
python gff_to_tsv_ncbi.py
```

### Examples
```bash
python gff_to_tsv_ncbi.py \
--gff_file SpneumoniaeD39V.gff \
--feature_types CDS tRNA rRNA ncRNA \
--output_tsv SpneumoniaeD39V.tsv
```

```bash
python gff_to_tsv_ncbi.py \
--gff_file LguyanensisM4147.gff \
--feature_types exon \
--output_tsv LguyanensisM4147.tsv
```

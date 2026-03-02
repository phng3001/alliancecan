# gff
## gff to gtf
### Binary
```bash
# -F: keep all GFF attributes, -T: output GTF
gffread input.gff -F -T -o output.gtf
```
### Container
```bash
apptainer run gffread_v0.12.7.sif gffread input.gff -F -T -o output.gtf
```

# Miscellaneous
## Replace string in each filename
```bash
# Replace 'genes' with 'features' in filenames
for f in *genes*; do mv -- "$f" "${f//genes/features}"; done
```
> Preview before renaming:
```bash
for f in *genes*; do echo mv -- "$f" "${f//genes/features}"; done
```

## Remove the first N lines
### tail
```bash
tail -n +$((N+1)) file.txt
```
### sed
```bash
sed '1,Nd' file.txt
```
### awk
```bash
awk 'NR>N' file.txt
```

## Extract FASTA entries longer than 1000 bases
```bash
grep -A1 '^>' sequences.fasta | awk 'NR%2==0 {if(length($0)>1000) print prev"\n"$0} {prev=$0}'
```

## Find files not belonging to specific users
```bash
find /path/to/directory -type f ! -user user1 ! -user user2
```

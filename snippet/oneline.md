# awk
## Print selected columns from input file
* By default, `awk` treats any whitespace (spaces or tabs) as the field separator (FS)
```bash
# 1st, 3rd and 5th columns
awk '{print $1,$3,$5}' input.txt
```
* Use `-F` flag  to specify the input field separator
```bash
# csv file
awk -F',' '{print $1,$3,$5}' input.csv
```
* Use `-v` flag and the built-in variable `OFS` to specify the output field separator
```bash
awk -v OFS="\t" '{print $1,$3,$5}' input.txt
```
> **Notes:** Use double quotes `""` to escape sequences like `\t` (tab), `\n` (new line), etc

### Example: Convert a CSV to a TSV by selecting and reordering columns
* Input: `data.csv`
```
gene_id,sample_id,expression_value,status
GENE_1,SAMPLE_A,45.6,OK
GENE_2,SAMPLE_B,88.8,FAIL
```
* Goal: Extract columns 2, 1, and 3 in that order, converting from CSV (comma-separated) to TSV (tab-separated)
* `awk` solution: 
```bash
awk -F',' -v OFS="\t" '{print $2, $1, $3}' data.csv
```
* Output:
```
sample_id   gene_id expression_value
SAMPLE_A	GENE_1	45.6
SAMPLE_B	GENE_2	88.8
```

## Filter by column value
### String
```bash
# Print lines where the 2nd column value is "folA"
awk -F"\t" '$2 == "folA"' genes.tsv
```
```bash
# case-insensitive
awk -F"\t" 'tolower($2) == "fola"' genes.tsv
```
```bash
# Print lines where the 2nd column value is not "folA"
awk -F"\t" '$2 != "folA"' genes.tsv
```
```bash
# Print lines where the 3rd column value starts with "spr"
awk '$3 ~ /^spr/' genes.tsv
# or
awk -v str="spr" '$3 ~ "^"str {print}' genes.tsv
```
```bash
# Print lines where the 3rd column value does not ends with "_gene"
awk '$3 !~ /_gene$/' genes.tsv
# or
awk -v str="_gene" '$3 !~ str"$"' genes.tsv
```

### Numeric value
```bash
# Print lines where the 3rd column value > 50
awk -F"\t" '$3 > 50' genes.tsv
```
```bash
# Print lines where the 3rd column value > the 5th column value
awk '$3 > $5' genes.tsv
```

### Combine multiple conditions
```bash
# Print lines where the 3rd column > 30 AND the 5th column < 80
awk -F"\t" '$3 > 30 && $5 < 80' genes.tsv
```
```bash
# Print lines where the 3rd column > 30 OR the 2nd column is "folA"
awk -F"\t" '$3 > 30 || $2 == "folA"' genes.tsv
```
```bash
# Skip header (#) AND print lines where the 3rd column < 50
awk -F"\t" '!/^#/ && $3 < 50' genes.tsv
```

### Regular expression
```bash
# Print lines where the 5th column matches the regular expression
# This regular expression means "start with a lowercase letter a–f"
awk '$5  ~ /^[a-f]/' file.txt
```
```bash
# Print lines where the 5th column does not match the regular expression
awk '$5  !~ /^[a-f]/' file.txt
```

## Duplicate handling 
```bash
# Removes exact duplicate lines (based on all columns)
awk '!seen[$0]++' file.txt
```
```bash
# Print duplicate lines (based on all columns)
awk 'seen[$0]++' file.txt
```
```bash
# Print unique value of the 1st column
awk '!arr[$1]++ {print $1}' file.txt
```
```bash
# Keep first line per unique value in the 2nd column
awk '!arr[$2]++' file.txt
```
```bash
# Keep last occurence per unique value in the 2nd column
tac file.txt | awk '!arr[$2]++' | tac
```
```bash
# Keep first line per unique combination of columns 1 and 2
awk '!seen[$1 FS $2]++' file.txt
```

## Compute sum & mean
* By default `awk` converts non-numeric strings to 0
```bash
# Print sum of the 1st column 
awk '{sum+=$1} END {print sum}' file.txt
```
```bash
# Print mean of the 2nd column
awk '{x+=$2} END {print x/NR}' file.txt # no header
# or
awk 'NR>1 {x+=$2; n++} END {if(n>0) print x/n}' file.txt # skip the first line
```

## Calculate the length of each sequence in a fasta file
```bash
awk 'BEGIN {OFS="\t"} /^>/ {if (seqlen) print seqid, seqlen; seqid = substr($1, 2); seqlen = 0; next} {seqlen += length($0)} END {if (seqlen) print seqid, seqlen}' file.fasta
```



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


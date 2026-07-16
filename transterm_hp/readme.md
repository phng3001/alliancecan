# Installation
```bash
unzip transterm_hp_v2.09.zip
cd transterm_hp_v2.09
sed -i "s/if(i == EOF) return false/if(i == EOF) return NULL/g" seq.cc
make clean transterm
make no_obj
```

# Usage
0. Prepare fasta and annotation files

Annotation file (extension`.coords` or `.crd`) is tab-separated and each line has the following format:
```
gene_name   start   end chrom_id
```
Cf. `test_data` folder for examples
1. Make symbolic link of `transterm` executable file to the working directory
2. Make `transterm` executable
```bash
chmod u+x transterm
```
3. Make symbolic link of `expterm.dat` file to the working directory
4. Run transtem
```bash
transterm -p expterm.dat seq.fasta annotation.coords > output.tt
```

# Source
Source file `transterm_hp_v2.09.zip` downloaded from https://transterm.cbcb.umd.edu/

# More info
Cf. `USAGE.txt` in the installation folder

ARNold: Another online tool for identification of Rho–independent transcription terminators (http://rssf.i2bc.paris-saclay.fr/toolbox/arnold/index.php)

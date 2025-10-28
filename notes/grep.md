# grep basic

```bash
grep [options] <search_string> <path/to/file>
```

## Print lines in file containing a string 
```bash 
grep folA file.txt
```
* Use single quotes `''` for literal patterns
```bash
# Literally search for the string $USER
grep '$USER' file.txt 
```
* Use double quotes `""` for variable and command expansion 
```bash
# Search for your username, as $USER is a Bash shell's default environment variable
grep "$USER" file.txt
```

## grep's options
### `-i`: Case insensitive
```bash
grep -i fola file.txt
```
### `-w`: Search only for the whole word
```bash
grep -w 100 file.txt
```
### `-v`: Return lines that don't match

### `-c`: Count the number of matching lines

### `-o`: Only print the matching part of a line

### `-l`: Print the name of the file(s) containing matches

### `-n`: Show the line numbers of the matches

### `-m`: Stop after max count is reached
```bash
grep -m 1 folA file.txt
```
### `-r`: Recursive search (within subdirectories)
```bash
grep -r <search_string> <path/to/directory>
```
### `-e`: Specify a pattern explicitly
Useful cases: 
* Pass multiple patterns
```bash
# Search for lines containing either "folA" or "folP"
grep -e "folA" -e "folP" file.txt
``` 
* When the pattern starts with a dash (-) that could be mistaken for an option
```bash
grep -e "-abc" file.txt
```
### `-E`: Extended Regular Expressions
So `|`, `+`, `?`, `{}`, `()` can be used **without escaping (backslashes)**.
```bash
grep -E "folA|folP" file.txt
```
instead of 
```bash
grep "folA\|folP" file.txt
```
### `-A`: Print lines after match
```bash
# Display 3 lines after matching string
grep -A 3 'Exception' error.log
```
### `-B`: Print lines before match
```bash
# Display 4 lines before matching string
grep -B 4 'Exception' error.log
```
### `-C`: Print lines around match
```bash
# Display 5 lines before and after matching string
grep -C 5 'Exception' error.log
```
### `-F`: Fixed string match (disables regex)
```bash
grep -F "LR1234.1" annotations.gff
```

### `-q`: Quiet mode, useful for scripting
This does not print any output to the terminal, only sets its **exit status**:
* `0`: Pattern found
* `1`: Pattern not found
* `>1`: Error (e.g. file not found)

Example:
```bash
grep -q "PASS" variants.vcf && echo "High-quality variants present" || echo "No PASS variants"
```

## Remove empty lines
```bash
grep . file.txt
```
or 
```bash
grep -v '^$' file.txt
```

## Remove empty or whitespace-only lines
```bash
grep -v '^[[:space:]]*$' file.txt
```

## Count empty lines
```bash
grep -c '^$' file.txt
```

## Count empty or whitespace-only lines
```bash
grep -c '^[[:space:]]*$' file.txt
```

## Remove lines starting with a string
```bash
# Remove lines which start with "#", like this line
grep -v '^#' file.txt
```

## Count occurences
```bash
grep -o folA file.txt | wc -l
```

## Search all files in directory
### Prints only the matching lines
```bash
grep -rh <search_string> <path/to/directory>
```
### Prints only the filenames with matches
```bash
grep -rl <search_string> <path/to/directory>
```

## Skip directories
```bash
grep folA -d skip 'dir_name' .
```
To avoid this error:
```bash
grep: dir_name: Is a directory
```

## Search patterns from a file (one per line)
```bash
grep -F -w -f gene_list.txt annotation.gff
```

## Detect lines containing ambiguous bases in FASTA file
```bash
grep -E '[^ATGCNatgcn]' sequences.fasta
```

## Extract FASTA entries longer than 1000 bases
```bash
grep -A1 '^>' sequences.fasta | awk 'NR%2==0 {if(length($0)>1000) print prev"\n"$0} {prev=$0}'
```

## More info
https://ryanstutorials.net/linuxtutorial/cheatsheetgrep.php

https://quickref.me/grep.html

https://devhints.io/grep
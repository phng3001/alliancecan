# tar basic

## Create an archive
```bash
tar -cf archive.tar file_1 file_2 file_n
```

## Compress an archive
```bash
# Using gzip (faster compression/decompression)
tar -cvzf archive.tar.gz <path_to_directory_or_file>
```

```bash
# Using bzip2 (better compression ratio i.e. smaller archive size)
tar -cvjf archive.tar.bz2 <path_to_directory_or_file>
```

```bash
# Using xz (even better compression ratio but slower)
tar -cvJf archive.tar.xz <path_to_directory_or_file>
```

## Extract an archive
```bash
# *.tar file
tar -xvf archive.tar
```

```bash
# *.tar.gz file
tar -xvzf archive.tar.gz
```

```bash
# *.tar.bz2 file
tar -xvjf archive.tar.bz2
```

```bash
# *.tar.xz file
tar -xvJf archive.tar.xz
```

## Extract an archive to a defined directory
```bash
tar -xvzf archive.tar.gz -C <path_to_directory>
```

## List archive content without extracting it
```bash
# Add -v for more details (like ls -l)
tar -tzf archive.tar.gz
```

## Extract a single file or directory from an archive
```bash
tar -xzf archive.tar.gz "file_name"
```

## Decompress an archive
```bash
# *.tar.gz file
gunzip archive.tar.gz
gunzip -c archive.tar.gz > archive.tar # to keep the original
```

```bash
# *.tar.bz2 file
# Add -k to keep the original
bzip2 -d archive.tar.bz2
```

```bash
# *.tar.xz file
# Add -k to keep the original
xz -d archive.tar.xz
```

## Append files to an existing archive
> **Notes:** This only works with uncompressed `.tar` archives
```bash
# Only update newer files
tar -uf archive.tar file1.txt dir2/ script3.sh
```
```bash
# Append files regardless of whether they exist
tar -rf archive.tar file1.txt dir2/ script3.sh
```

## Append an archive to another one
> **Notes:** This only works with uncompressed `.tar` archives
```bash
# This appends archive2.tar content into archive.tar
# Only one archive can be appended at a time
tar -Af archive.tar archive2.tar
```

## More info
https://linux-audit.com/cheat-sheets/tar/

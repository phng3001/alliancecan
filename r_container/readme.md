# R images
Container | Source | Description
----------|--------|------------
rocker_verse.sif | docker://rocker/verse | tidyverse, devtools, tex & publishing-related packages
rocker_verse_DESeq2_v0.0.3.sif | PNP | docker://rocker/verse + DESeq2 & related packages

# Launch RStudio in an Apptainer image

## Check if rserver is available in the image (required)
```bash
apptainer exec r_image.sif which rserver
```

## HPC

### In the working directory, create local directories for RStudio Server to write logs, databases, session info
```bash
mkdir -p rstudio/var/lib
mkdir -p rstudio/var/run
mkdir -p rstudio/tmp
```

### Start an interactive session on a compute node with salloc
```bash
salloc --time=1:0:0 --ntasks=1 --cpus-per-task=2 --mem-per-cpu=4Gb --account=def-professor
```

### Launch rserver on the compute node
```bash
apptainer exec \
  --bind rstudio/var/lib:/var/lib/rstudio-server \
  --bind rstudio/var/run:/var/run/rstudio-server \
  --bind rstudio/tmp:/tmp \
  r_image.sif \
  rserver \
    --www-port=8888 \
    --server-user=$(whoami) \
    --server-daemonize=0 \
    --auth-none=1
```

### SSH tunneling
```bash
ssh -L 8888:somecomputenode:8888 someuser@someserver.alliancecan.ca
```

### Open localhost on a web browser
```
http://localhost:8888
```
RStudio will be launched automatically

## WSL

### In the working directory, create a directory to bind to RStudio Server 
```bash
mkdir -p rstudio/project
```

### Launch rserver
```bash
apptainer exec \
  --writable-tmpfs \
  --bind rstudio/project/:/home/rstudio/project \
  r_image.sif \
  rserver \
    --www-address=0.0.0.0 \
    --www-port=8888 \
    --server-user=$(whoami) \
    --server-daemonize=0 \
    --auth-none=1
```
### Open localhost on a web browser
```
http://localhost:8888
```
RStudio will be launched automatically

> **Notes:** 
> - Only the bound directories and files inside them are persistent, others will be erased once the session closed 
> - In R scripts don't forget to set the working directory in a bound directory (E.g. `setwd("/home/rstudio/project")` in this example)

# Ressources
https://rocker-project.org/images/

# R images
Container | Source | Description
----------|--------|------------
rocker_verse.sif | docker://rocker/verse | tidyverse, devtools, tex & publishing-related package

# Launch RStudio in Apptainer image

## Check if rserver is available in the image (required)
```bash
apptainer exec r_image.sif which rserver
```

## In the working directory, create local directories for RStudio Server to write logs, databases, session info
```bash
mkdir -p rstudio/var/lib
mkdir -p rstudio/var/run
mkdir -p rstudio/tmp
```

## Start an interactive session on a compute node with salloc
```bash
salloc --time=1:0:0 --ntasks=1 --cpus-per-task=2 --mem-per-cpu=4Gb --account=def-someuser
```

## Launch rserver on the compute node
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

## SSH tunneling
```bash
ssh -L 8888:rc12601:8888 someuser@rorqual.alliancecan.ca
```

## Open localhost on web browser
```
http://localhost:8888
```

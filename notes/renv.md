---
title:  "renv"
excerpt: "renv basic"
collection: cheatsheets
---

## Install / Load R
```bash
# On alliancecan clusters
module load StdEnv/2023 r/4.5.0
```

## Open R and install `renv` package
```R
> install.packages("renv")
```
> **Notes:** This only needs to be done once per server per R version

## Open R and navigate to the project directory
```R
> setwd("/path/to/directory")
```

## Initialize renv for the project
```R
> renv::init()
```

## Restart the R session
```R
> q()
```

## Work in the project directory
### Activate `renv`
> **Notes:** This should be done everytime

```R
> setwd("/path/to/directory")
> renv::activate()
```

or

```R
> setwd("/path/to/directory")
> source("renv/activate.R")
```

### Example: Install ggplot2
```R
> install.packages("ggplot2") 
```

## Save the state of the project dependencies to the lockfile (`renv.lock`)
> **Notes:** Only packages being used inside R scripts in the project will be saved

```R
> renv::snapshot()
```

## Restore the project dependencies from the lockfile (`renv.lock`)
```R
> renv::restore()
```

## Other useful commands
### Check `renv` status
```R
> renv::status()
```
### Check dependencies in a project
```R
> renv::dependencies()
```

### Deactivate `renv`
```R
> renv::deactivate()
```

## More info
https://rstudio.github.io/renv/articles/renv.html

https://docs.posit.co/ide/user/ide/guide/environments/r/renv.html



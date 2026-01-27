# Load R
```bash
module load StdEnv/2023 r/4.5.0
```

# Initialize renv
```R
> renv::init()
The following package(s) will be updated in the lockfile:

# CRAN -----------------------------------------------------------------------
- renv   [* -> 1.1.6]

The version of R recorded in the lockfile will be updated:
- R      [* -> 4.5.0]

- Lockfile written to "/project/6002430/Scripts_MOU/PNP/alliancecan/deseq2/renv.lock".
- renv activated -- please restart the R session.
> q()
Save workspace image? [y/n/c]: n
```

# Install packages
```bash
Rscript packages.R
```

# Save the project's state
```bash
Rscript renv_snapshot.R
```

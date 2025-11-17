# Create a shared directory

```bash
mkdir shared_dir_name
chgrp def-mouellet shared_dir_name # def-professor
chmod -R 2775 shared_dir_name # owner = rwx, group = rwx, others = r-x
setfacl -R -m g::rwx shared_dir_name
setfacl -d -m g::rwx shared_dir_name
setfacl -d -m o::r-x shared_dir_name
```

# Data stream
Operator | Function | Syntax
---------|----------|-------
`>` | Redirect stdout | `command > file.txt` or `command 1> file.txt`
`2>` | Redirect stderr | `command 2> error.txt`
`&>` | Redirect all | `command &> log.txt`
`>>` | Append stdout | `command >> files.txt`
`2>>` | Append stderr | `command 2>> errors.txt`
`&>>` | Append all | `command &> log.txt`

> **Notes:** Redirect to `/dev/null` to hide the corresponding stream (E.g. `command &> /dev/null`)


# Google Colab
## Upload
```python
from google.colab import files
upload = files.upload()
```

## Download
```python
from google.colab import files
download = files.download('output.csv')
```
# Create a shared directory

```bash
mkdir shared_dir_name
chgrp def-mouellet shared_dir_name # def-professor
chmod -R 2775 shared_dir_name # owner = rwx, group = rwx, others = r-x
setfacl -R -m g::rwx shared_dir_name
setfacl -d -m g::rwx shared_dir_name
setfacl -d -m o::r-x shared_dir_name
```

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
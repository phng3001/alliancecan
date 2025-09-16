# Install Aspera Connect
> **Notes:** This only needs to be done once per server.
## Run the setup script
```bash
bash ibm-aspera-connect_4.2.13.820_linux_x86_64.sh
```
## Make the library files executable
```bash
chmod u+x ~/.aspera/connect/plugins/*/*.so ~/.aspera/connect/lib/*
```
## Run the `setrpaths` script
```bash
setrpaths.sh --path $HOME/.aspera
```
## Add the ascp binaries to PATH (optional)
```bash
export PATH=~/.aspera/connect/bin:$PATH
```

# Upload files via Aspera
```bash
bash run_aspera_4+.sh <aspera_key_path> <upload_dir_path>
```
* `aspera_key_path`: path to aspera ssh key
> **Notes:** aspera ssh key is regularly updated by NCBI. **Be sure to have the lastest version of the key**.
* `upload_dir_path`: path to the local folder that contains all of the files to upload
> **Notes:** Sequence files for one submission should be in one directory. Do not create complex directory structures *i.e.* directories within directories.

# More info
https://docs.alliancecan.ca/wiki/Transferring_files_with_Aspera_Connect/ascp
https://www.ncbi.nlm.nih.gov/sra/docs/submitfiles/#aspera-connect
https://submit.ncbi.nlm.nih.gov/subs/sra/#files_upload_aspera_cmd

# Notes
> 2025-09-16
- [ ] Nibi
```
Session Stop  (Error: Client unable to connect to server (check UDP port and firewall))
```
- [x] Fir
- [x] Rorqual
- [x] Narval

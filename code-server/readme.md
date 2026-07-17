# Containers
Container | Source | Description
----------|--------|------------
code-server.sif|docker://codercom/code-server:latest|Version 4.129.0-39

# Launch VSCode from a code-server image
## HPC
### Launch VSCode
```bash
apptainer run code-server.sif --auth none
```
### SSH tunneling
```bash
ssh -L 8888:somenode:8080 someuser@someserver.alliancecan.ca
```
### Open localhost on a web browser
```
http://localhost:8888
```
> **Notes:** Server and localhost ports can be customized. By default code-server uses port `8080`. To customize server port, for example use port `8888` instead of the default port:
```bash
apptainer run code-server.sif --auth none --bind-addr 0.0.0.0:8888
```

## WSL
### Launch VSCode using a different port than `8080`
```bash
apptainer run code-server.sif --auth none --bind-addr 0.0.0.0:8889
```
### Open localhost on a web browser
```
http://localhost:8889
```

# Ressources
https://hub.docker.com/r/codercom/code-server

https://github.com/coder/code-server

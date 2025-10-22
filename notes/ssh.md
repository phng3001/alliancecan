# SSH basic

## `ssh`: Connect to a remote server
```bash
ssh username@servername
```

## `ssh-keygen`: Generating SSH keys
```bash
# Generate an RSA 4096-bit key with "your_email@example.com" as a comment
ssh-keygen -t rsa -b 4096 -C "your_email@example.com"
```
* `-t`: Type of key (`rsa`, `ed25519`, `dsa`, `ecdsa`)
* `-b`: Number of bits in the key
* `-C`: Comment
> **Notes:** `ed25519` is often the best choice (small, fast, and secure). Its key size is fixed at 256 bits

### GitHub
```bash
ssh-keygen -t ed25519 -C "your_email@example.com"
```
https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent#generating-a-new-ssh-key

## `scp`: Copy files securely between servers
### Remote to local
```bash
scp user@server:/folder/file.ext dest/
```
### Local to remote
```bash
scp dest/file.ext user@server:/folder
```
### Server to server
```bash
scp user1@server1:/folder/file.ext user2@server2:/folder
```

## `rsync`: Copy files securely between servers, preserving file attributes
```bash
rsync -av -e ssh user@source:/path/to/dir/ user@dest:/path/to/dir/
```

## More info
https://www.stationx.net/ssh-commands-cheat-sheet/
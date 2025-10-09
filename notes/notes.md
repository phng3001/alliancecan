# SSH

## Basic

### `ssh`: Connect to a remote server
```bash
ssh username@servername
```

### `ssh-keygen`: Generating SSH keys
```bash
# Generate an RSA 4096-bit key with "your_email@example.com" as a comment
ssh-keygen -t rsa -b 4096 -C "your_email@example.com"
```
* `-t`: Type of key (`rsa`, `ed25519`, `dsa`, `ecdsa`)
* `-b`: Number of bits in the key
* `-C`: Comment
> **Notes:** `ed25519` is often the best choice (small, fast, and secure). Its key size is fixed at 256 bits

#### GitHub
```bash
ssh-keygen -t ed25519 -C "your_email@example.com"
```
https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent#generating-a-new-ssh-key

### `scp`: Copy files securely between servers
#### Remote to local
```bash
scp user@server:/folder/file.ext dest/
```
#### Local to remote
```bash
scp dest/file.ext user@server:/folder
```
#### Server to server
```bash
scp user1@server1:/folder/file.ext user2@server2:/folder
```

### `rsync`: Copy files securely between servers, preserving file attributes
```bash
rsync -av -e ssh user@source:/path/to/dir/ user@dest:/path/to/dir/
```

## More info
https://www.stationx.net/ssh-commands-cheat-sheet/



# git

## Setup

* Set the name that will be attached to your commits and tags
```bash
git config --global user.name "your_username"
```
* Set the e-mail address that will be attached to your commits and tags
```bash
git config --global user.email "your_email@example.com"
```
* Set automatic command line coloring for Git
```bash
git config --global color.ui auto
```
* Set git's default editor (e.g. vim, nano, notepad etc)
```bash
git config --global core.editor "nano"
```
* Local directory, single project (default)
```bash
git config --local user.email "your_email@example.com"
```
* List all global configurations set
```bash
git config --global --list
```
* Check a specific global setting (e.g. user.name)
```bash
git config --global <key>
```
* View git global config file
```bash
cat ~/.gitconfig
```
* Check a specific configuration key (e.g. user.email)
```bash
git config --get <key>
```

## Init

* Initializes a new Git repository in the current directory
```bash
git init
```
* Creates a new Git repository in the specified directory
```bash
git init <directory>
```
* Retrieve an entire repository from a hosted location via URL
```bash
git clone <repository_url>
```
* Clones a specific branch from a remote repository
```bash
git clone --branch <branch_name> <repository_url>
```
* Show your current local branch
```bash
git branch
```

## Stage
* Show modified files in working directory
```bash
git status
```
* Adds a specific file to the staging area
```bash
git add <file>
```
* Adds all modified and new files to the staging area
```bash
git add .
```
* Commit your staged content
```bash
git commit -m "your_message"
```
* Correct a git commit message if the commit has not been pushed yet
```bash
git commit --amend
```
* Show difference of what is changed but not staged
```bash
git diff
```
* Show difference of what is staged but not yet committed
```bash
git diff --staged
```

## Remote repositories
* Lists the names of the remotes
```bash
git remote
```
* Lists remote names + URLs for fetch and push
```bash
git remote -v
```
* Pushes local commits to the remote repository
```bash
git push
```
* Pushes local commits to the specified remote repository
```bash
git push <remote>
```
* Pushes local commits to the specified branch of the remote repository
```bash
git push <remote> <branch>
```
E.g. 
```bash
git push origin main
```
* Retrieves change from a remote repository
```bash
git fetch
```
* See commit messages from a remote repository
```bash
git fetch
# then
git log HEAD..origin/main
# or
git log HEAD..origin/main --oneline
# or
git log HEAD..origin/main --pretty=format:"%h %an %s"
```
* See difference between your branch and the remote repository
```bash
git fetch
# then
# single combined diff
git diff HEAD..origin/main
# full commit history
git log HEAD..origin/main -p
```
* Fetches changes from the remote repository and merges them into the current branch
```bash
git pull
```

## More info
https://www.geeksforgeeks.org/git/git-cheat-sheet/



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
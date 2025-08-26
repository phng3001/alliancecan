# SSH key

## Generating a new SSH key
ssh-keygen -t ed25519 -C "your_email@example.com"

## More info
https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent#generating-a-new-ssh-key



# git

## Setup

* Set the name that will be attached to your commits and tags<br>
git config --global user.name "your_username"
* Set the e-mail address that will be attached to your commits and tags<br>
git config --global user.email "your_email"
* Set automatic command line coloring for Git<br>
git config --global color.ui auto
* Set git's default editor (e.g. vim, nano, notepad etc)<br>
git config --global core.editor "nano"
* Local directory, single project (default)<br>
git config --local user.email "your_email"
* List all global configurations set<br>
git config --global --list
* Check a specific global setting<br>
git config --global user.name
* View git global config file<br>
cat ~/.gitconfig

## Init

* Initializes a new Git repository in the current directory<br>
git init
* Creates a new Git repository in the specified directory<br>
git init [directory]
* Retrieve an entire repository from a hosted location via URL<br>
git clone [repository_url]
* Clones a specific branch from a remote repository<br>
git clone --branch [branch_name] [repository_url]

## Stage
* Show modified files in working directory<br>
git status
* Adds a specific file to the staging area<br>
git add [file]
* Adds all modified and new files to the staging area.<br>
git add .
* Commit your staged content<br>
git commit -m "your_message"
* Show difference of what is changed but not staged<br>
git diff
* Show difference of what is staged but not yet committed<br>
git diff --staged

## Remote repositories
* Lists the names of the remotes<br>
git remote
* Lists remote names + URLs for fetch and push<br>
git remote -v
* Pushes local commits to the remote repository<br>
git push
* Pushes local commits to the specified remote repository<br>
git push [remote]
* Pushes local commits to the specified branch of the remote repository<br>
git push [remote] [branch]<br>
E.g. git push origin main
* Retrieves change from a remote repository<br>
git fetch








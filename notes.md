# SSH key

## Generating a new SSH key
ssh-keygen -t ed25519 -C "your_email@example.com"

## More info
https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent#generating-a-new-ssh-key



# git

## Setup

* Set the name that will be attached to your commits and tags
git config --global user.name "[your_username]"
* Set the e-mail address that will be attached to your commits and tags
git config --global user.email "[your_email]"
* Set automatic command line coloring for Git
git config --global color.ui auto
* Set git's default editor (e.g. vim, nano, notepad etc)
git config --global core.editor "nano"
* Local directory, single project (default)
git config --local user.email "[your_email]"
* List all global configurations set
git config --global --list
* Check a specific global setting
git config --global user.name
* View git global config file
cat ~/.gitconfig





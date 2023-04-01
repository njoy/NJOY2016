---
layout: page
---
# Jekyll-website
This repository contains the required stuff to make the NJOY website and keep it consistent across all the projects.

## Using jekyll-website
For all projects under the NJOY organization, this repository should be created as a git 'subtree' in the `docs` directory. This will ensure consistency across all the projects as well as make the documentation for each accessible from the website. 

In order to prevent pushing to this repository one should set the push url for this remote in git like this:

```bash
git remote set-url --push <name> no_push
```
Here `<name>` is the name of the remote (probably `jekyll`). What this tells git is whenever you try to push to this remote, simply push to the url `no_push` which doesn't exist. 

Alternatively, you can set the push url to the url of the remote of the real repository like:

```bash
git remote set-url --push jekyll git@github.com:njoy/njoy.github.io.git
```
as is done for the [njoy.github.io](https://github.com/njoy/njoy.github.io) repository.

## License
This repository---and all repositories under the NJOY organization---are licensed according to the [LICENSE](LICENSE.md) file

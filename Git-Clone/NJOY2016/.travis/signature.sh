#! /bin/bash

if [ "$TRAVIS_PULL_REQUEST" = false ] && 
   [ "$TRAVIS_BRANCH" = "master" ] &&
   [ "$build_type" = "coverage" ]; then
  git config --global user.email "jlconlin@lanl.gov"
  git config --global user.name "Jeremy Lloyd Conlin"

  timestamp=`date +%F_%T`

  git clone https://github.com/njoy/signatures.git

  project=$(basename $TRAVIS_REPO_SLUG)
  DIR="$PWD/signatures/$project"
  mkdir -p $DIR

  filename=$DIR/$timestamp

  ./metaconfigure/signature.py $filename

  cd $DIR
  git checkout master
  git add $filename.json
  git commit -m "Adding signature file for $project."

  git remote add origin-travis https://user:${GH_TOKEN}@github.com/njoy/signatures.git > /dev/null 2>&1

git push --quiet --set-upstream origin-travis master
fi

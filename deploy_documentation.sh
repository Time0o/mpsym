#!/bin/bash

USER_NAME="Time0o"
REPO_NAME="TUD_computational_group_theory"
REPO_URL="github.com/$USERNAME/$REPO_NAME.git"

echo "Cloning gh-pages branch..."
git clone -b gh-pages "https://$REPO_TOKEN@$REPO_URL"
cd "$REPO_NAME"

echo "Configuring git..."
git config push.default simple
git config user.name "Timo Nicolai"
git config user.name "timonicolai@arcor.de"

if [ -d ../doxygen/html ] && [ -f ../doxygen/html/index.html ]; then
  echo "Pushing documentation..."
  cp -r ../doxygen/html .
  git add .
  git commit --amend -m "Deploy documentation, current commit is $TRAVIS_COMMIT"
  git push --force "https://$REPO_TOKEN@$REPO_URL" > /dev/null 2>&1
else
  echo "Failed to find generated documentation" >&2
  exit 1
fi

#!/bin/bash

USER_NAME="Time0o"
REPO_NAME="mpsym"
REPO_URL="github.com/$USER_NAME/$REPO_NAME.git"

if [ "${MATRIX_EVAL}" == "CC=gcc-7 && CXX=g++-7" ]; then
  # upload coverage data
  echo "Uploading coverage data..."
  mkdir gcov
  cd gcov
  gcov ../source/CMakeFiles/**/*
  bash <(curl -s https://codecov.io/bash)
  cd ..

  # deploy documentation
  echo "Deploying documentation..."
  cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_DOC=ON ..
  make doxygen

  echo "Cloning gh-pages branch..."
  git clone -b gh-pages "https://$REPO_URL"
  cd "$REPO_NAME"

  echo "Configuring git..."
  git config push.default simple
  git config user.name "Timo Nicolai"
  git config user.name "timonicolai@arcor.de"

  if [ -d ../doxygen/html ] && [ -f ../doxygen/html/index.html ]; then
    echo "Pushing documentation..."
    cp -r ../doxygen/html/* .
    git add .
    git commit --amend -m "Deploy documentation, current commit is $TRAVIS_COMMIT"
    git push --force "https://$REPO_TOKEN@$REPO_URL" > /dev/null 2>&1
  else
    echo "Failed to find generated documentation" >&2
    exit 1
  fi
fi

#!/bin/bash

AUTHOR="Timo Nicolai"
AUTHOR_EMAIL="timonicolai@arcor.de"

USER_NAME="Time0o"
REPO_NAME="mpsym"
REPO_URL="https://github.com/$USER_NAME/$REPO_NAME.git"

# upload coverage data

if [ -z "$CODECOV_TOKEN" ]; then
  echo >&2 "CODECOV_TOKEN must be set"
  exit 1
fi

mkdir -p build_coverage
(
  cd build_coverage
  cmake .. -DCMAKE_BUILD_TYPE=Debug -DPYTHON_BINDINGS=ON -DCOVERAGE=ON
  make -j $(nproc)
  make test

  mkdir -p gcov
  (
    cd gcov
    echo "Generating coverage data"
    gcov ../source/CMakeFiles/**/* > gcov.log 2> gcov.err
    echo "Uploading coverage data"
    bash <(curl -s https://codecov.io/bash) > codecov.log 2> codecov.err
  )
)

# upload documentation
mkdir -p build_documentation
(
  cd build_documentation
  cmake .. -DCMAKE_BUILD_TYPE=Release -DNO_INSTALL=ON -DBUILD_DOC=ON
  make -j $(nproc)
  make doc
  echo "Cloning gh-pages branch..."
  git clone -b gh-pages "$REPO_URL"
  cd "$REPO_NAME"
  echo "Pushing documentation..."
  cp -r ../doxygen/html/* .
  git add .
  git commit --amend -m "Deploy documentation, current commit is $(git rev-parse HEAD)"
  git push --force "$REPO_URL" > /dev/null 2>&1
)

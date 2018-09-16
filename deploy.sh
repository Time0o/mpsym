#!/bin/bash

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
  cmake -DCMAKE_BUILD_TYPE=Release ..
  make doxygen
  bash "${TRAVIS_BUILD_DIR}"/deploy_documentation.sh
fi

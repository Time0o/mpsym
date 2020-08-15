cmake_minimum_required(VERSION 3.6)

project(nauty-download NONE)

include(ExternalProject)

set(EXTRA_CFLAGS "'-O3 -fPIC'")
set(NCPU 4)

ExternalProject_Add(nauty_traces
  URL               http://pallini.di.uniroma1.it/nauty26r10.tar.gz
  SOURCE_DIR        "${NAUTY_WORK_DIR}"
  BINARY_DIR        "${NAUTY_WORK_DIR}"
  CONFIGURE_COMMAND "${NAUTY_WORK_DIR}/configure" "CFLAGS=${EXTRA_CFLAGS}"
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
)

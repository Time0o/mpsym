cmake_minimum_required(VERSION 3.6)

project(nauty-download NONE)

include(ExternalProject)

set(NCPU 4)

ExternalProject_Add(nauty_traces
  URL               http://pallini.di.uniroma1.it/nauty26r10.tar.gz
  SOURCE_DIR        "${NAUTY_WORK_DIR}"
  BINARY_DIR        "${NAUTY_WORK_DIR}"
  CONFIGURE_COMMAND "${NAUTY_WORK_DIR}/configure"
  BUILD_COMMAND     make -C "${NAUTY_WORK_DIR}" nauty.a -j "${NCPU}" &&
                    mv "${NAUTY_WORK_DIR}/nauty.a" "${NAUTY_LIB}"
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
)

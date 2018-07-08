cmake_minimum_required(VERSION 3.6)

project(nauty-download NONE)

include(ExternalProject)

set(NCPU 12)
ExternalProject_Add(nauty
  URL               http://pallini.di.uniroma1.it/nauty26r10.tar.gz
  SOURCE_DIR        "${NAUTY_WORKDIR}"
  BINARY_DIR        "${NAUTY_WORKDIR}"
  CONFIGURE_COMMAND "${NAUTY_WORKDIR}/configure"
  BUILD_COMMAND     make -C "${NAUTY_WORKDIR}" nauty.a -j "${NCPU}"
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
)

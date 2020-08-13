cmake_minimum_required(VERSION 3.6)

project(pybind11-download NONE)

include(ExternalProject)

ExternalProject_Add(pybind11
  GIT_REPOSITORY    "https://github.com/pybind/pybind11.git"
  GIT_TAG           "master"
  SOURCE_DIR        "${PYBIND11_SRC_DIR}"
  BINARY_DIR        "${PYBIND11_BIN_DIR}"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
)

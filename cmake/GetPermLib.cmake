cmake_minimum_required(VERSION 3.6)

project(permlib-download NONE)

include(ExternalProject)

ExternalProject_Add(permlib
  GIT_REPOSITORY    "https://github.com/tremlin/PermLib.git"
  GIT_TAG           "v0.2.9"
  SOURCE_DIR        "${PERMLIB_WORK_DIR}"
  BINARY_DIR        "${PERMLIB_WORK_DIR}"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
)

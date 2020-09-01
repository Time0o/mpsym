cmake_minimum_required(VERSION 3.6)

project(nlohmann_json-download NONE)

include(ExternalProject)

ExternalProject_Add(nlohmann_json
  GIT_REPOSITORY    "https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent"
  GIT_TAG           "master"
  SOURCE_DIR        "${NLOHMANN_JSON_SRC_DIR}"
  BINARY_DIR        "${NLOHMANN_JSON_BIN_DIR}"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
)

cmake_minimum_required(VERSION 3.6)

project(CGTL)

### Download Google Test #######################################################

string(CONCAT GTEST_CMAKE
  "cmake_minimum_required(VERSION 3.0)\n"
  "\n"
  "project(googletest-download NONE)\n"
  "include(ExternalProject)\n"
  "\n"
  "ExternalProject_Add(googletest\n"
  "  GIT_REPOSITORY    https://github.com/google/googletest.git\n"
  "  GIT_TAG           master\n"
  "  SOURCE_DIR        \"\${CMAKE_BINARY_DIR}/googletest-src\"\n"
  "  BINARY_DIR        \"\${CMAKE_BINARY_DIR}/googletest-build\"\n"
  "  CONFIGURE_COMMAND \"\"\n"
  "  BUILD_COMMAND     \"\"\n"
  "  INSTALL_COMMAND   \"\"\n"
  "  TEST_COMMAND      \"\"\n"
  ")")

set (GTEST_CMAKE_TMPFILE "GTestCMakeList.txt")

file(WRITE "${GTEST_CMAKE_TMPFILE}" "${GTEST_CMAKE}")
configure_file("${GTEST_CMAKE_TMPFILE}" googletest-download/CMakeLists.txt)
file(REMOVE "${GTEST_CMAKE_TMPFILE}")

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )

if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )

if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
                 ${CMAKE_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

### Add subdirectories #########################################################

add_subdirectory(source)
add_subdirectory(test)

enable_testing()

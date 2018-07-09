cmake_minimum_required(VERSION 3.6)


################################################################################
# Project Structure
################################################################################

project(CGTL "CXX")
set(MAIN_LIB "cgtl")

set(CMAKE_SCRIPT_DIR "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH "${CMAKE_SCRIPT_DIR}")

set(NAUTY_CONFIG "${CMAKE_SCRIPT_DIR}/GetNauty.cmake")
set(GTEST_CONFIG "${CMAKE_SCRIPT_DIR}/GetGTest.cmake")

file(DOWNLOAD
     https://raw.githubusercontent.com/bilke/cmake-modules/master/CodeCoverage.cmake
     "${CMAKE_SCRIPT_DIR}/CodeCoverage.cmake")

set(SRC_DIR "${PROJECT_SOURCE_DIR}/source")
set(INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

set(TEST_SRC_DIR "${PROJECT_SOURCE_DIR}/test/source")
set(TEST_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/test/include")

set(NAUTY_DIR "${CMAKE_BINARY_DIR}/nauty")
set(NAUTY_DOWNLOAD_DIR "${NAUTY_DIR}/download")
set(NAUTY_WORKDIR "${NAUTY_DIR}/workdir")

set(GTEST_DIR "${CMAKE_BINARY_DIR}/gtest")
set(GTEST_DOWNLOAD_DIR "${GTEST_DIR}/download")
set(GTEST_SRC_DIR "${GTEST_DIR}/source")
set(GTEST_BIN_DIR "${GTEST_DIR}/build")


################################################################################
# External Dependencies
################################################################################

function(add_external_project)
  set(options NONE)
  set(args NAME CONFIGFILE WORKDIR)
  cmake_parse_arguments(Arg "${options}" "${args}" "" "${ARGN}")

  configure_file("${Arg_CONFIGFILE}" "${Arg_WORKDIR}/CMakeLists.txt")

  execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY "${Arg_WORKDIR}")

  if(result)
    message(FATAL_ERROR "CMake step for ${Arg_NAME} failed: ${result}")
  endif()

  execute_process(COMMAND "${CMAKE_COMMAND}" --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY "${Arg_WORKDIR}")

  if(result)
    message(FATAL_ERROR "Build step for ${Arg_NAME} failed: ${result}")
  endif()
endfunction()

# Boost
find_package(Boost 1.40 REQUIRED COMPONENTS graph)

# Nauty
add_external_project(NAME "nauty/traces"
                     CONFIGFILE "${NAUTY_CONFIG}"
                     WORKDIR "${NAUTY_DOWNLOAD_DIR}")

# Google Test
add_external_project(NAME "googletest"
                     CONFIGFILE "${GTEST_CONFIG}"
                     WORKDIR "${GTEST_DOWNLOAD_DIR}")

add_subdirectory("${GTEST_SRC_DIR}" "${GTEST_BIN_DIR}" EXCLUDE_FROM_ALL)


################################################################################
# Testing
################################################################################

enable_testing()

# Code coverage
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  include(CodeCoverage)

  set(COVERAGE_TARGET "coverage")

  # TODO

  #file(GLOB_RECURSE COVERAGE_STDLIB_BOOST_EXCLUDES
  #     "${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES}/*")

  #set(COVERAGE_EXCLUDES "${COVERAGE_STDLIB_BOOST_EXCLUDES}")

  setup_target_for_coverage(NAME "${COVERAGE_TARGET}"
                            EXECUTABLE "${CMAKE_CTEST_COMMAND}")
endif()


################################################################################
# Profiling
################################################################################

# TODO


################################################################################
# Add Subdirectories
################################################################################

add_subdirectory("${SRC_DIR}")
add_subdirectory("${TEST_SRC_DIR}")

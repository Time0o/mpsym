cmake_minimum_required(VERSION 3.6)

project(CGTL "CXX")
set (CMAKE_CXX_STANDARD 11)


################################################################################
# Project Structure
################################################################################

set(MAIN_LIB "cgtl")

set(CMAKE_SCRIPT_DIR "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH "${CMAKE_SCRIPT_DIR}")

set(NAUTY_CONFIG "${CMAKE_SCRIPT_DIR}/GetNauty.cmake")
set(GTEST_CONFIG "${CMAKE_SCRIPT_DIR}/GetGTest.cmake")

set(SRC_DIR "${PROJECT_SOURCE_DIR}/source")
set(INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

set(TEST_SRC_DIR "${PROJECT_SOURCE_DIR}/test/source")
set(TEST_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/test/include")
set(TEST_RESOURCE_DIR "${PROJECT_SOURCE_DIR}/test/resources")

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

#Lua
find_package(Lua 5.3 REQUIRED)

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
# Code Coverage
################################################################################

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  find_program(GCOV_PATH gcov)
  find_program(LCOV_PATH  NAMES lcov lcov.bat lcov.exe lcov.perl)
  find_program(GENHTML_PATH NAMES genhtml genhtml.perl genhtml.bat)

  if(NOT GCOV_PATH)
    message(FATAL_ERROR "Failed to find gcov")
  endif()

  if(NOT LCOV_PATH)
    message(FATAL_ERROR "Failed to find lcov")
  endif()

  if(NOT GENHTML_PATH)
    message(FATAL_ERROR "Failed to find genhtml")
  endif()

  set(Coverage_TARGET "coverage")
  set(Coverage_OUTDIR "coverage")

  add_custom_target(${Coverage_TARGET}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}

    # Cleanup lcov
    COMMAND ${LCOV_PATH} -z -d . --include "'${PROJECT_SOURCE_DIR}/*'"

    # Run tests
    COMMAND ${CMAKE_CTEST_COMMAND}

    # Capturing lcov counters and generating report
    COMMAND ${LCOV_PATH} -c -d . -o "lcov.info" --include "'${PROJECT_SOURCE_DIR}/*'"

    # generate html output
    COMMAND ${GENHTML_PATH} "lcov.info" -o "${Coverage_OUTDIR}"

    # remove temporary files
    COMMAND ${CMAKE_COMMAND} -E remove "lcov.base" "lcov.info" "lcov.total"
  )
endif()


################################################################################
# Profiling
################################################################################

# TODO


################################################################################
# Add Subdirectories
################################################################################

add_subdirectory("${SRC_DIR}")

enable_testing()
add_subdirectory("${TEST_SRC_DIR}")

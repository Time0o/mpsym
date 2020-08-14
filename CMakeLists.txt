################################################################################
# Project Settings
################################################################################

cmake_minimum_required(VERSION 3.6)

set(PROJECT_AUTHOR "Timo Nicolai")
set(PROJECT_COPYRIGHT "2020")
set(PROJECT_LICENSE "MIT")

set(MPSYM "mpsym")
set(PYMPSYM "pympsym")
set(DESCRIPTION "MPSoC Symmetry Reduction")
set(VERSION "1.0")

add_definitions(-DMPSYM="${MPSYM}" -DPYMPSYM=${PYMPSYM}
                -DDESCRIPTION="${DESCRIPTION}" -DVERSION="${VERSION}")

project("${MPSYM}"
        DESCRIPTION "${DESCRIPTION}"
        VERSION "${VERSION}"
        LANGUAGES "CXX")

set (CMAKE_CXX_STANDARD 11)


################################################################################
# Project Structure
################################################################################

set(CMAKE_SCRIPT_DIR "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH "${CMAKE_SCRIPT_DIR}")

set(NAUTY_CONFIG "${CMAKE_SCRIPT_DIR}/GetNauty.cmake")
set(GTEST_CONFIG "${CMAKE_SCRIPT_DIR}/GetGTest.cmake")
set(PERMLIB_CONFIG "${CMAKE_SCRIPT_DIR}/GetPermLib.cmake")
set(PYBIND11_CONFIG "${CMAKE_SCRIPT_DIR}/GetPybind11.cmake")

set(SRC_DIR "${PROJECT_SOURCE_DIR}/source")
set(INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
set(MPSYM_LIB "${MPSYM}")

set(TEST_SRC_DIR "${PROJECT_SOURCE_DIR}/test/source")
set(TEST_TESTS_DIR "${PROJECT_SOURCE_DIR}/test/tests")
set(TEST_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/test/include")

set(PROFILE_SRC_DIR "${PROJECT_SOURCE_DIR}/profile/source")
set(PROFILE_COMMON_DIR "${PROJECT_SOURCE_DIR}/profile/common")
set(PROFILE_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/profile/include")

set(NAUTY_DIR "${CMAKE_BINARY_DIR}/nauty")
set(NAUTY_DOWNLOAD_DIR "${NAUTY_DIR}/download")
set(NAUTY_WORKDIR "${NAUTY_DIR}/workdir")
set(NAUTY_LIB "${NAUTY_WORKDIR}/libnauty.a")

set(GTEST_DIR "${CMAKE_BINARY_DIR}/gtest")
set(GTEST_DOWNLOAD_DIR "${GTEST_DIR}/download")
set(GTEST_SRC_DIR "${GTEST_DIR}/source")
set(GTEST_BIN_DIR "${GTEST_DIR}/build")

set(PERMLIB_DIR "${CMAKE_BINARY_DIR}/permlib")
set(PERMLIB_DOWNLOAD_DIR "${PERMLIB_DIR}/download")
set(PERMLIB_WORKDIR "${PERMLIB_DIR}/workdir")
set(PERMLIB_INCLUDE_DIR "${PERMLIB_WORKDIR}/include")

set(PYBIND11_DIR "${CMAKE_BINARY_DIR}/pybind11")
set(PYBIND11_DOWNLOAD_DIR "${PYBIND11_DIR}/download")
set(PYBIND11_SRC_DIR "${PYBIND11_DIR}/source")
set(PYBIND11_BIN_DIR "${PYBIND11_DIR}/build")

set(PYTHON_SRC_DIR "${PROJECT_SOURCE_DIR}/python")
set(PYTHON_MODULE "${PYMPSYM}")

set(INSTALL_LIB_DIR "lib")
set(INSTALL_INCLUDE_DIR "include/${MPSYM}")

set(DOC_DOXYGEN_IN_DIR "${PROJECT_SOURCE_DIR}/doxygen")
set(DOC_DOXYGEN_OUT_DIR "${CMAKE_BINARY_DIR}/doxygen")


################################################################################
# Common Compilation Settings
################################################################################

set(MAIN_COMPILER_FLAGS "-Wall -Wextra -pedantic")
set(MAIN_LINK_FLAGS "")

function(set_default_target_properties)
  set(options NONE)
  set(args TARGET)
  cmake_parse_arguments(Arg "${options}" "${args}" "" "${ARGN}")

  set_target_properties("${Arg_TARGET}" PROPERTIES
    LINKER_LANGUAGE CXX
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    COMPILE_FLAGS "${MAIN_COMPILER_FLAGS}"
    LINK_FLAGS "${MAIN_LINK_FLAGS}"
  )
endfunction()


################################################################################
# Convencience Function for Adding External Dependencies
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


################################################################################
# Common Dependencies
################################################################################

# Boost
find_package(Boost 1.40 REQUIRED COMPONENTS graph)

# Lua
find_package(Lua 5.2 REQUIRED)

# Nauty
add_external_project(NAME "nauty_traces"
                     CONFIGFILE "${NAUTY_CONFIG}"
                     WORKDIR "${NAUTY_DOWNLOAD_DIR}")


################################################################################
# Debug Build
################################################################################

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  # Google Test
  add_external_project(NAME "googletest"
                       CONFIGFILE "${GTEST_CONFIG}"
                       WORKDIR "${GTEST_DOWNLOAD_DIR}")

  add_subdirectory("${GTEST_SRC_DIR}" "${GTEST_BIN_DIR}" EXCLUDE_FROM_ALL)

  # Coverage
  set(DO_COVERAGE TRUE)

  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(DO_COVERAGE FALSE)
    message(WARNING "Coverage generation currently only available with g++")
  endif()

  find_program(GCOV_PATH gcov)
  if(NOT GCOV_PATH)
    set(DO_COVERAGE FALSE)
    message(WARNING "Failed to find gcov")
  endif()

  if(DO_COVERAGE)
    set(MAIN_COMPILER_FLAGS "${MAIN_COMPILER_FLAGS} -g -O0 --coverage")
    set(MAIN_LINK_FLAGS "${MAIN_LINK_FLAGS} --coverage")

    find_program(LCOV_PATH lcov)
    find_program(GENHTML_PATH genhtml)

    if(NOT LCOV_PATH)
      set(LCOV_COVERAGE FALSE)
      message(WARNING "Failed to find lcov")
    endif()

    if(NOT GENHTML_PATH)
      set(LCOV_COVERAGE FALSE)
      message(WARNING "Failed to find genhtml")
    endif()

    if(LCOV_COVERAGE)
      set(Coverage_TARGET "lcov")
      set(Coverage_OUTDIR "lcov")

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
  endif()

  # Enable and build tests
  enable_testing()

  add_subdirectory("${TEST_TESTS_DIR}")
endif()


################################################################################
# Profile Build
################################################################################

if(CMAKE_BUILD_TYPE STREQUAL "Profile")
  # Permlib
  add_external_project(NAME "permlib"
                       CONFIGFILE "${PERMLIB_CONFIG}"
                       WORKDIR "${PERMLIB_DOWNLOAD_DIR}")

  set(MAIN_COMPILER_FLAGS "${MAIN_COMPILER_FLAGS} -g -O2 -DNDEBUG")

  # Build profiling programs
  add_subdirectory("${PROFILE_SRC_DIR}")
endif()


################################################################################
# Common Sources
################################################################################

add_subdirectory("${SRC_DIR}")


################################################################################
# Bindings
################################################################################

if(PYTHON_BINDINGS STREQUAL "Yes")
  # pybind11
  add_external_project(NAME "pybind11"
                       CONFIGFILE "${PYBIND11_CONFIG}"
                       WORKDIR "${PYBIND11_DOWNLOAD_DIR}")

  add_subdirectory("${PYBIND11_SRC_DIR}" "${PYBIND11_BIN_DIR}" EXCLUDE_FROM_ALL)

  add_subdirectory("${PYTHON_SRC_DIR}")
endif()


################################################################################
# Release Build
################################################################################

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  # Installation
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(FATAL_ERROR "CMAKE_INSTALL_PREFIX must be specified explicitly")
  endif()

  file(GLOB PUBLIC_HEADERS "${INCLUDE_DIR}/*.hpp")

  set_target_properties("${MPSYM_LIB}" PROPERTIES
    PUBLIC_HEADER "${PUBLIC_HEADERS}"
  )

  install(FILES "${NAUTY_LIB}"
    DESTINATION "${INSTALL_LIB_DIR}"
  )

  install(TARGETS "${MPSYM_LIB}"
    LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
    PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}"
  )

  # Documentation
  message(STATUS "Configuring documentation build")

  set(BUILD_DOC TRUE)

  find_package(Doxygen)
  if(NOT DOXYGEN_FOUND)
    set(BUILD_DOC FALSE)
    message(WARNING "Failed to find Doxygen")
  else()
    message(STATUS "Found Doxygen")
  endif()

  if (NOT BUILD_DOC)
    message(WARNING "Documentation will not be buildable")
  else()
    # Doxygen
    set(DOXYFILE_IN "${DOC_DOXYGEN_IN_DIR}/Doxyfile.in")
    set(DOXYFILE_OUT "${DOC_DOXYGEN_OUT_DIR}/Doxyfile")

    configure_file("${DOXYFILE_IN}" "${DOXYFILE_OUT}" @ONLY)

    add_custom_target(doxygen
      COMMAND "${DOXYGEN_EXECUTABLE}" "${DOXYFILE_OUT}"
      VERBATIM
    )

    add_custom_target(doc DEPENDS doxygen)
  endif()
endif()

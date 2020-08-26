################################################################################
# Project Settings
################################################################################

cmake_minimum_required(VERSION 3.6)

set(CMAKE_PROJECT_AUTHOR "Timo Nicolai")
set(CMAKE_PROJECT_COPYRIGHT "2020")
set(CMAKE_PROJECT_URL "https://github.com/Time0o/mpsym")
set(CMAKE_PROJECT_DOC_URL "https://time0o.github.io/mpsym/")
set(CMAKE_PROJECT_LICENSE "MIT")

set(MPSYM "mpsym")
set(PYMPSYM "pympsym")

project(
  "${MPSYM}"
  VERSION "0.1"
  DESCRIPTION "MPSoC Symmetry Reduction"
  HOMEPAGE_URL "https://github.com/Time0o/mpsym"
  LANGUAGES "CXX" "C"
)

add_definitions(
  -DMPSYM="${MPSYM}"
  -DPYMPSYM=${PYMPSYM}
  -DDESCRIPTION="${CMAKE_PROJECT_DESCRIPTION}"
  -DVERSION="${CMAKE_PROJECT_VERSION}"
)

set (CMAKE_CXX_STANDARD 11)


################################################################################
# Project Structure
################################################################################

# CMake build types
set(CMAKE_BUILD_TYPE_DEBUG "Debug")
set(CMAKE_BUILD_TYPE_PROFILE "Profile")
set(CMAKE_BUILD_TYPE_RELEASE "Release")

# CMake scripts
set(CMAKE_SCRIPT_DIR "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH "${CMAKE_SCRIPT_DIR}")

set(CMAKE_ADD_EXTERNAL_PROJECT "${CMAKE_SCRIPT_DIR}/AddExternalProject.cmake")
set(CMAKE_RUN_DOXYGEN "${CMAKE_SCRIPT_DIR}/RunDoxygen.cmake")

set(CMAKE_NAUTY_CONFIG "${CMAKE_SCRIPT_DIR}/GetNauty.cmake")
set(CMAKE_GTEST_CONFIG "${CMAKE_SCRIPT_DIR}/GetGTest.cmake")
set(CMAKE_PERMLIB_CONFIG "${CMAKE_SCRIPT_DIR}/GetPermLib.cmake")
set(CMAKE_PYBIND11_CONFIG "${CMAKE_SCRIPT_DIR}/GetPybind11.cmake")

# mpsym
set(MPSYM_SRC_DIR "${PROJECT_SOURCE_DIR}/source")
set(MPSYM_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
set(MPSYM_LIB "${MPSYM}")
set(MPSYM_INSTALL_LIB_DIR "lib")
set(MPSYM_INSTALL_INCLUDE_DIR "include/${MPSYM}")

# mpsym tests
set(MPSYM_TEST_SRC_DIR "${PROJECT_SOURCE_DIR}/test/source")
set(MPSYM_TEST_TESTS_DIR "${PROJECT_SOURCE_DIR}/test/tests")
set(MPSYM_TEST_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/test/include")

# mpsym profiling
set(MPSYM_PROFILE_SRC_DIR "${PROJECT_SOURCE_DIR}/profile/source")
set(MPSYM_PROFILE_COMMON_DIR "${PROJECT_SOURCE_DIR}/profile/common")
set(MPSYM_PROFILE_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/profile/include")

# nauty
set(NAUTY_SRC_DIR "${PROJECT_SOURCE_DIR}/nauty")
set(NAUTY_BIN_DIR "${CMAKE_BINARY_DIR}/nauty")
set(NAUTY_DOWNLOAD_DIR "${NAUTY_BIN_DIR}/download")
set(NAUTY_WORK_DIR "${NAUTY_BIN_DIR}/workdir")
set(NAUTY_LIB "nauty")
set(NAUTY_INSTALL_LIB_DIR "lib")

# googletest
set(GTEST_BIN_DIR "${CMAKE_BINARY_DIR}/gtest")
set(GTEST_DOWNLOAD_DIR "${GTEST_BIN_DIR}/download")
set(GTEST_SRC_DIR "${GTEST_BIN_DIR}/source")
set(GTEST_BIN_DIR "${GTEST_BIN_DIR}/build")
set(GTEST_GMOCK_MAIN "gmock_main")

# gcov
set(GCOV_LIB "gcov")

# permlib
set(PERMLIB_BIN_DIR "${CMAKE_BINARY_DIR}/permlib")
set(PERMLIB_DOWNLOAD_DIR "${PERMLIB_BIN_DIR}/download")
set(PERMLIB_WORK_DIR "${PERMLIB_BIN_DIR}/workdir")
set(PERMLIB_INCLUDE_DIR "${PERMLIB_WORK_DIR}/include")

# pybind11
set(PYBIND11_BIN_DIR "${CMAKE_BINARY_DIR}/pybind11")
set(PYBIND11_DOWNLOAD_DIR "${PYBIND11_BIN_DIR}/download")
set(PYBIND11_SRC_DIR "${PYBIND11_BIN_DIR}/source")
set(PYBIND11_BIN_DIR "${PYBIND11_BIN_DIR}/build")

set(PYTHON_SRC_DIR "${PROJECT_SOURCE_DIR}/python")
set(PYTHON_BIN_DIR "${CMAKE_BINARY_DIR}/python")
set(PYTHON_SETUP_PY_IN "${PYTHON_SRC_DIR}/setup.py.in")
set(PYTHON_SETUP_PY_OUT "${PYTHON_BIN_DIR}/setup.py")
set(PYTHON_MODULE "${PYMPSYM}")

# Doxygen
set(DOC_DOXYGEN_IN_DIR "${PROJECT_SOURCE_DIR}/doxygen")
set(DOC_DOXYGEN_OUT_DIR "${CMAKE_BINARY_DIR}/doxygen")
set(DOC_DOXYGEN_INDEX_HTML "${DOC_DOXYGEN_OUT_DIR}/html/index.html")


################################################################################
# CMake Functions
################################################################################

include("${CMAKE_ADD_EXTERNAL_PROJECT}")


################################################################################
# External Dependencies
################################################################################

# Boost
message(STATUS "Finding boost...")
find_package(Boost 1.40 REQUIRED COMPONENTS graph)

# Lua
message(STATUS "Finding Lua...")
find_package(Lua 5.2 REQUIRED)

# Nauty
message(STATUS "Acquiring Nauty...")
add_external_project(NAME "nauty_traces"
                     CONFIG_FILE "${CMAKE_NAUTY_CONFIG}"
                     WORK_DIR "${NAUTY_DOWNLOAD_DIR}")

add_subdirectory("${NAUTY_SRC_DIR}")

# googletest
if(CMAKE_BUILD_TYPE STREQUAL "${CMAKE_BUILD_TYPE_DEBUG}")
  message(STATUS "Acquiring googletest...")
  add_external_project(NAME "googletest"
                       CONFIG_FILE "${CMAKE_GTEST_CONFIG}"
                       WORK_DIR "${GTEST_DOWNLOAD_DIR}"
                       SRC_DIR "${GTEST_SRC_DIR}"
                       BIN_DIR "${GTEST_BIN_DIR}")
endif()


# Permlib
if(CMAKE_BUILD_TYPE STREQUAL "${CMAKE_BUILD_TYPE_PROFILE}")
  message(STATUS "Acquiring PermLib...")
  add_external_project(NAME "permlib"
                       CONFIG_FILE "${CMAKE_PERMLIB_CONFIG}"
                       WORK_DIR "${PERMLIB_DOWNLOAD_DIR}")
endif()

# pybind11
if(PYTHON_BINDINGS)
  message(STATUS "Acquiring pybind11...")
  add_external_project(NAME "pybind11"
                       CONFIG_FILE "${CMAKE_PYBIND11_CONFIG}"
                       WORK_DIR "${PYBIND11_DOWNLOAD_DIR}"
                       SRC_DIR "${PYBIND11_SRC_DIR}"
                       BIN_DIR "${PYBIND11_BIN_DIR}")
endif()


################################################################################
# Debug Build
################################################################################

if(CMAKE_BUILD_TYPE STREQUAL "${CMAKE_BUILD_TYPE_DEBUG}")
  # Coverage
  message(STATUS "Configuring coverage data generation...")

  set(DO_COVERAGE TRUE)

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    find_program(GCOV_PATH gcov)
    if(NOT GCOV_PATH)
      set(DO_COVERAGE FALSE)
      message(STATUS "Failed to find gcov")
    endif()
  else()
    set(DO_COVERAGE FALSE)
    message(STATUS "Coverage data generation currently only available with g++")
  endif()

  if(DO_COVERAGE)
    message(STATUS "Coverage data will be generated")

    set(MAIN_COMPILER_FLAGS "${MAIN_COMPILER_FLAGS} -g -O0 --coverage")
    set(MAIN_LINK_FLAGS "${MAIN_LINK_FLAGS} --coverage")
  else()
    message(STATUS "Coverage data will be generated")
  endif()

  # Enable and build tests
  enable_testing()

  add_subdirectory("${MPSYM_TEST_TESTS_DIR}")
endif()


################################################################################
# Profile Build
################################################################################

if(CMAKE_BUILD_TYPE STREQUAL "${CMAKE_BUILD_TYPE_PROFILE}")
  add_subdirectory("${MPSYM_PROFILE_SRC_DIR}")
endif()


################################################################################
# Release Build
################################################################################

if(CMAKE_BUILD_TYPE STREQUAL "${CMAKE_BUILD_TYPE_RELEASE}")
  # Prefix check
  message(STATUS "Checking install prefix...")

  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(FATAL_ERROR "CMAKE_INSTALL_PREFIX must be specified explicitly")
  endif()

  # Documentation
  message(STATUS "Configuring documentation...")

  set(BUILD_DOC TRUE)

  find_package(Doxygen)
  if(NOT DOXYGEN_FOUND)
    set(BUILD_DOC FALSE)
    message(STATUS "Failed to find Doxygen")
  else()
    message(STATUS "Found Doxygen")
  endif()

  if(BUILD_DOC)
    message(STATUS "Documentation will be generated")

    # Doxygen
    set(DOXYFILE_IN "${DOC_DOXYGEN_IN_DIR}/Doxyfile.in")
    set(DOXYFILE_OUT "${DOC_DOXYGEN_OUT_DIR}/Doxyfile")

    configure_file("${DOXYFILE_IN}" "${DOXYFILE_OUT}" @ONLY)

    add_custom_command(
      OUTPUT "${DOC_DOXYGEN_INDEX_HTML}"
      COMMAND "${DOXYGEN_EXECUTABLE}" "${DOXYFILE_OUT}"
      DEPENDS "${MPSYM_LIB}"
    )

    add_custom_target(doc DEPENDS "${DOC_DOXYGEN_INDEX_HTML}")
  else()
    message(STATUS "Documentation will not be generated")
  endif()
endif()


################################################################################
# Bindings
################################################################################

# Python
if(PYTHON_BINDINGS)
  add_subdirectory("${PYTHON_SRC_DIR}")
endif()


################################################################################
# Source
################################################################################

add_subdirectory("${MPSYM_SRC_DIR}")

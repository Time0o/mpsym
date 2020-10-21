cmake_minimum_required(VERSION 3.6)


################################################################################
# Project Settings
################################################################################

set(CMAKE_PROJECT_NAME "mpsym")
set(CMAKE_PROJECT_VERSION "0.2")
set(CMAKE_PROJECT_DESCRIPTION "MPSoC Symmetry Reduction")
set(CMAKE_PROJECT_HOMEPAGE_URL "https://github.com/Time0o/mpsym")
set(CMAKE_PROJECT_DOC_URL "https://time0o.github.io/mpsym/")
set(CMAKE_PROJECT_AUTHOR "Timo Nicolai")
set(CMAKE_PROJECT_AUTHOR_EMAIL "timonicolai@arcor.de")
set(CMAKE_PROJECT_COPYRIGHT "2020")
set(CMAKE_PROJECT_LICENSE "MIT")

project(
  "${CMAKE_PROJECT_NAME}"
  VERSION "${CMAKE_PROJECT_VERSION}"
  DESCRIPTION "${CMAKE_PROJECT_DESCRIPTION}"
  HOMEPAGE_URL "${CMAKE_PROJECT_HOMEPAGE_URL}"
  LANGUAGES "CXX" "C"
)

set(CMAKE_CXX_STANDARD 11)

set(DEBUG_FLAGS " -g3 -fno-omit-frame-pointer")
string(APPEND CMAKE_C_FLAGS_DEBUG "${DEBUG_FLAGS}")
string(APPEND CMAKE_CXX_FLAGS_DEBUG "${DEBUG_FLAGS}")
string(APPEND CMAKE_C_FLAGS_RELWITHDEBINFO "${DEBUG_FLAGS}")
string(APPEND CMAKE_CXX_FLAGS_RELWITHDEBINFO "${DEBUG_FLAGS}")

################################################################################
# Project Structure
################################################################################

# CMake build types
set(CMAKE_BUILD_TYPE_DEBUG "Debug")
set(CMAKE_BUILD_TYPE_PROFILE "Profile")
set(CMAKE_BUILD_TYPE_RELEASE "Release")

# CMake scripts
set(CMAKE_SCRIPT_DIR "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH "${CMAKE_SCRIPT_DIR}")

set(CMAKE_ADD_EXTERNAL_PROJECT "${CMAKE_SCRIPT_DIR}/AddExternalProject.cmake")
set(CMAKE_RUN_DOXYGEN "${CMAKE_SCRIPT_DIR}/RunDoxygen.cmake")

set(CMAKE_NAUTY_CONFIG "${CMAKE_SCRIPT_DIR}/GetNauty.cmake")
set(CMAKE_NLOHMANN_JSON_CONFIG "${CMAKE_SCRIPT_DIR}/GetNlohmannJson.cmake")
set(CMAKE_GTEST_CONFIG "${CMAKE_SCRIPT_DIR}/GetGTest.cmake")
set(CMAKE_PERMLIB_CONFIG "${CMAKE_SCRIPT_DIR}/GetPermLib.cmake")
set(CMAKE_PYBIND11_CONFIG "${CMAKE_SCRIPT_DIR}/GetPybind11.cmake")

# Other scripts
set(EMBED_FILE "${CMAKE_SOURCE_DIR}/embedfile.sh")

# mpsym
set(MPSYM_SRC_DIR "${CMAKE_SOURCE_DIR}/source")
set(MPSYM_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
set(MPSYM_LIB "${CMAKE_PROJECT_NAME}")
set(MPSYM_INSTALL_LIB_DIR "lib")
set(MPSYM_INSTALL_INCLUDE_DIR "include/${MPSYM}")

# mpsym tests
set(MPSYM_TEST_SRC_DIR "${CMAKE_SOURCE_DIR}/test/source")
set(MPSYM_TEST_TESTS_DIR "${CMAKE_SOURCE_DIR}/test/tests")
set(MPSYM_TEST_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/test/include")

# mpsym profiling
set(MPSYM_PROFILE_SRC_DIR "${CMAKE_SOURCE_DIR}/profile/source")
set(MPSYM_PROFILE_COMMON_DIR "${CMAKE_SOURCE_DIR}/profile/common")
set(MPSYM_PROFILE_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/profile/include")

# nlohmann_json
set(NLOHMANN_JSON_BIN_DIR "${CMAKE_BINARY_DIR}/nlohmann_json")
set(NLOHMANN_JSON_DOWNLOAD_DIR "${NLOHMANN_JSON_BIN_DIR}/download")
set(NLOHMANN_JSON_SRC_DIR "${NLOHMANN_JSON_BIN_DIR}/source")

# nauty
set(NAUTY_SRC_DIR "${CMAKE_SOURCE_DIR}/nauty")
set(NAUTY_BIN_DIR "${CMAKE_BINARY_DIR}/nauty")
set(NAUTY_DOWNLOAD_DIR "${NAUTY_BIN_DIR}/download")
set(NAUTY_WORK_DIR "${NAUTY_BIN_DIR}/workdir")
set(NAUTY_LIB "nauty")
set(NAUTY_INSTALL_LIB_DIR "lib")

# googletest
set(GTEST_BIN_DIR "${CMAKE_BINARY_DIR}/gtest")
set(GTEST_DOWNLOAD_DIR "${GTEST_BIN_DIR}/download")
set(GTEST_SRC_DIR "${GTEST_BIN_DIR}/source")
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

# Lua
set(LUA_SRC_DIR "${CMAKE_SOURCE_DIR}/lua")
set(LUA_BIN_DIR "${CMAKE_BINARY_DIR}/lua")
set(LUA_MODULE "${CMAKE_PROJECT_NAME}")
set(LUA_MODULE_PATH "${LUA_SRC_DIR}/${LUA_MODULE}.lua")
set(LUA_MODULE_EMBED "${LUA_BIN_DIR}/${LUA_MODULE}_lua_embed.c")
set(LUA_ROCKSPEC_IN "${LUA_SRC_DIR}/${LUA_MODULE}.rockspec.in")
set(LUA_ROCKSPEC_OUT "${LUA_BIN_DIR}/${LUA_MODULE}-${CMAKE_PROJECT_VERSION}-1.rockspec")

# Python
set(PYTHON_SRC_DIR "${CMAKE_SOURCE_DIR}/python")
set(PYTHON_BIN_DIR "${CMAKE_BINARY_DIR}/python")
set(PYTHON_SETUP_PY_IN "${PYTHON_SRC_DIR}/setup.py.in")
set(PYTHON_SETUP_PY_OUT "${PYTHON_BIN_DIR}/setup.py")
set(PYTHON_MODULE "py${CMAKE_PROJECT_NAME}")

# Doxygen
set(DOC_DOXYGEN_IN_DIR "${CMAKE_SOURCE_DIR}/doxygen")
set(DOC_DOXYGEN_OUT_DIR "${CMAKE_BINARY_DIR}/doxygen")
set(DOC_DOXYGEN_INDEX_HTML "${DOC_DOXYGEN_OUT_DIR}/html/index.html")


################################################################################
# Definitions
################################################################################

add_definitions(
  -DPYTHON_MODULE=${PYTHON_MODULE}
  -DDESCRIPTION="${CMAKE_PROJECT_DESCRIPTION}"
  -DVERSION="${CMAKE_PROJECT_VERSION}"
)


################################################################################
# External Dependencies
################################################################################

include("${CMAKE_ADD_EXTERNAL_PROJECT}")

# Boost
if(LINK_STATIC)
  set(Boost_USE_STATIC_LIBS ON)
endif()

message(STATUS "Finding boost...")
find_package(Boost 1.40 REQUIRED COMPONENTS graph)

# Lua
message(STATUS "Finding Lua...")
find_package(Lua 5.2 REQUIRED)

# LuaRocks
if(CMAKE_BUILD_TYPE STREQUAL "${CMAKE_BUILD_TYPE_RELEASE}" AND
   NOT NO_INSTALL AND NOT LUA_NO_ROCK AND NOT LUA_EMBED)
  message(STATUS "Finding LuaRocks...")
  find_program(LUAROCKS_EXECUTABLE luarocks)
  if(NOT LUAROCKS_EXECUTABLE)
    message(FATAL_ERROR "LuaRocks not found")
  endif()
endif()

# Nauty
message(STATUS "Acquiring Nauty...")
add_external_project(NAME "nauty_traces"
                     CONFIG_FILE "${CMAKE_NAUTY_CONFIG}"
                     WORK_DIR "${NAUTY_DOWNLOAD_DIR}")

add_subdirectory("${NAUTY_SRC_DIR}")

# nlohmann_json
message(STATUS "Acquiring nlohmann_json...")
add_external_project(NAME "nlohmann_json"
                     CONFIG_FILE "${CMAKE_NLOHMANN_JSON_CONFIG}"
                     WORK_DIR "${NLOHMANN_JSON_DOWNLOAD_DIR}"
                     SRC_DIR "${NLOHMANN_JSON_SRC_DIR}")

set(JSON_BuildTests OFF CACHE INTERNAL "")

add_subdirectory("${NLOHMANN_JSON_SRC_DIR}")

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

if(PYTHON_BINDINGS)
  # Python
  message(STATUS "Finding Python...")
  if(NOT PYTHON_EXECUTABLE)
    find_program(PYTHON_EXECUTABLE python3)
    if(NOT PYTHON_EXECUTABLE)
      message(FATAL_ERROR "Python not found")
    endif()
  endif()

  # pybind11
  message(STATUS "Acquiring pybind11...")
  add_external_project(NAME "pybind11"
                       CONFIG_FILE "${CMAKE_PYBIND11_CONFIG}"
                       WORK_DIR "${PYBIND11_DOWNLOAD_DIR}"
                       SRC_DIR "${PYBIND11_SRC_DIR}"
                       BIN_DIR "${PYBIND11_BIN_DIR}")

  # pip
  if(NOT PYTHON_NO_SETUP_PY AND CMAKE_BUILD_TYPE STREQUAL "${CMAKE_BUILD_TYPE_RELEASE}")
    find_program(PIP_EXECUTABLE pip3)
    if(NOT PIP_EXECUTABLE)
      message(FATAL_ERROR "pip not found")
    endif()
  endif()
endif()


################################################################################
# Debug Build
################################################################################

if(CMAKE_BUILD_TYPE STREQUAL "${CMAKE_BUILD_TYPE_DEBUG}")
  if (NOT NO_COVERAGE)
    # Coverage
    message(STATUS "Configuring coverage data generation...")

    set(DO_COVERAGE TRUE)

    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      find_program(GCOV_EXECUTABLE gcov)
      if(NOT GCOV_EXECUTABLE)
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
  if(NOT NO_INSTALL AND NOT IGNORE_INSTALL_PREFIX AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(FATAL_ERROR "CMAKE_INSTALL_PREFIX must be specified explicitly")
  endif()

  if(BUILD_DOC)
    # Prefix check
    message(STATUS "Checking install prefix...")

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
endif()


################################################################################
# Bindings
################################################################################

# Lua
add_subdirectory("${LUA_SRC_DIR}")

# Python
if(PYTHON_BINDINGS)
  add_subdirectory("${PYTHON_SRC_DIR}")
endif()


################################################################################
# Source
################################################################################

add_subdirectory("${MPSYM_SRC_DIR}")

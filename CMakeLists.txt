cmake_minimum_required(VERSION 3.6)

project("MPSoC Symmetry Reduction Library"
        VERSION 0.1
        LANGUAGES "CXX")

set (CMAKE_CXX_STANDARD 11)

set(PROJECT_AUTHOR "Timo Nicolai")
set(PROJECT_COPYRIGHT "2018")


################################################################################
# Project Structure
################################################################################

set(MAIN_LIB "cgtl")

set(CMAKE_SCRIPT_DIR "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH "${CMAKE_SCRIPT_DIR}")

set(NAUTY_CONFIG "${CMAKE_SCRIPT_DIR}/GetNauty.cmake")
set(GTEST_CONFIG "${CMAKE_SCRIPT_DIR}/GetGTest.cmake")
set(PERMLIB_CONFIG "${CMAKE_SCRIPT_DIR}/GetPermLib.cmake")

set(SRC_DIR "${PROJECT_SOURCE_DIR}/source")
set(INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
set(RESOURCE_DIR "${CMAKE_BINARY_DIR}/resources")

set(TEST_SRC_DIR "${PROJECT_SOURCE_DIR}/test/source")
set(TEST_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/test/include")
set(TEST_RESOURCE_DIR "${PROJECT_SOURCE_DIR}/test/resources")

set(PROFILE_SRC_DIR "${PROJECT_SOURCE_DIR}/profile/source")
set(PROFILE_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/profile/include")
set(PROFILE_RESOURCE_DIR "${PROJECT_SOURCE_DIR}/profile/resources")
set(PROFILE_BIN_DIR "${CMAKE_BINARY_DIR}/profile/source")

set(NAUTY_DIR "${CMAKE_BINARY_DIR}/nauty")
set(NAUTY_DOWNLOAD_DIR "${NAUTY_DIR}/download")
set(NAUTY_WORKDIR "${NAUTY_DIR}/workdir")

set(GTEST_DIR "${CMAKE_BINARY_DIR}/gtest")
set(GTEST_DOWNLOAD_DIR "${GTEST_DIR}/download")
set(GTEST_SRC_DIR "${GTEST_DIR}/source")
set(GTEST_BIN_DIR "${GTEST_DIR}/build")

set(PERMLIB_DIR "${CMAKE_BINARY_DIR}/permlib")
set(PERMLIB_DOWNLOAD_DIR "${PERMLIB_DIR}/download")
set(PERMLIB_WORKDIR "${PERMLIB_DIR}/workdir")
set(PERMLIB_INCLUDE_DIR "${PERMLIB_WORKDIR}/include")

set(DOC_DOXYGEN_IN_DIR "${PROJECT_SOURCE_DIR}/doxygen")
set(DOC_DOXYGEN_OUT_DIR "${CMAKE_BINARY_DIR}/doxygen")
set(DOC_SPHINX_IN_DIR "${PROJECT_SOURCE_DIR}/sphinx")
set(DOC_SPHINX_OUT_DIR "${CMAKE_BINARY_DIR}/sphinx")


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
find_package(Lua 5.2 REQUIRED)

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
# Documentation
################################################################################

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  message(STATUS "Configuring documentation build")

  set(BUILD_DOC TRUE)

  find_package(Doxygen)
  if(NOT DOXYGEN_FOUND)
    set(BUILD_DOC FALSE)
    message(WARNING "Failed to find Doxygen")
  else()
    message(STATUS "Found Doxygen")
  endif()

  find_program(SPHINX_EXECUTABLE "sphinx-build")
  if(SPHINX_EXECUTABLE STREQUAL "sphinx-build-NOTFOUND")
    set(BUILD_DOC FALSE)
    message(WARNING "Failed to find Sphinx")
  else()
    message(STATUS "Found Sphinx")
  endif()

  if (NOT BUILD_DOC)
    message(WARNING "Documentation will not be buildable")
  else()
    set(DOXYFILE_IN "${DOC_DOXYGEN_IN_DIR}/Doxyfile")
    set(DOXYFILE_OUT "${DOC_DOXYGEN_OUT_DIR}/Doxyfile")

    configure_file("${DOXYFILE_IN}" "${DOXYFILE_OUT}" @ONLY)

    string(CONCAT SPHINX_CONF
           "# AUTOMATICALLY GENERATED BY CMAKE\n\n"
           "project = '${CMAKE_PROJECT_NAME}'\n"
           "author = '${PROJECT_AUTHOR}'\n"
           "copyright = '${PROJECT_COPYRIGHT}, ${PROJECT_AUTHOR}'\n"
           "version = '${PROJECT_VERSION_MAJOR}'\n"
           "release = '${PROJECT_VERSION}'\n"
           "master_doc = 'index'\n"
           "source_suffix = '.rst'\n"
           "language = 'cpp'\n"
           "highlight_language = 'cpp'\n"
           "primary_domain = 'cpp'\n"
           "pygments_style = 'sphinx'\n"
           "extensions = [\n"
           "  'sphinx.ext.githubpages',\n"
           "  'sphinx.ext.imgmath',\n"
           "  'breathe',\n"
           "  'exhale'\n"
           "]\n"
           "htmlhelp_basename = 'MPSoCSymmetryReductionLibrarydoc'\n"
           "html_theme = 'sphinxdoc'\n"
           "imgmath_image_format = 'svg'\n"
           "breathe_projects = {\n"
           "  '${PROJECT_NAME}' : '${DOC_DOXYGEN_OUT_DIR}/xml'\n"
           "}\n"
           "breathe_default_project = '${PROJECT_NAME}'\n"
           "breathe_domain_by_extension = {'h' : 'cpp'}\n"
           "exhale_args = {\n"
           "  'containmentFolder'   : './api',\n"
           "  'rootFileName'        : 'library_root.rst',\n"
           "  'rootFileTitle'       : 'Library API',\n"
           "  'doxygenStripFromPath': '..',\n"
           "  'createTreeView'      : True\n"
           "}")

    file(WRITE "${DOC_SPHINX_IN_DIR}/conf.py" "${SPHINX_CONF}")

    add_custom_target(doxygen
      COMMAND "${DOXYGEN_EXECUTABLE}" "${DOXYFILE_OUT}"
      VERBATIM
    )

    add_custom_target(sphinx
      COMMAND "${SPHINX_EXECUTABLE}" -b html
              "${DOC_SPHINX_IN_DIR}" "${DOC_SPHINX_OUT_DIR}"
      VERBATIM
      DEPENDS doxygen
    )

    message("Build the documentation with 'make doc'")
  endif()
endif()


################################################################################
# Code Coverage
################################################################################

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  find_program(GCOV_PATH gcov)
  find_program(LCOV_PATH  NAMES lcov lcov.bat lcov.exe lcov.perl)
  find_program(GENHTML_PATH NAMES genhtml genhtml.perl genhtml.bat)

  set(DO_COVERAGE TRUE)

  if(NOT GCOV_PATH)
    set(DO_COVERAGE FALSE)
    message(WARNING "Failed to find gcov")
  endif()

  if(DO_COVERAGE)
    set(MAIN_COMPILER_FLAGS "${MAIN_COMPILER_FLAGS} -g -O0 --coverage")
    set(MAIN_LINK_FLAGS "${MAIN_LINK_FLAGS} --coverage")
  endif()

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


################################################################################
# Profiling
################################################################################

if(CMAKE_BUILD_TYPE STREQUAL "Profile")
  find_program(PERF_PATH perf)

  if(NOT PERF_PATH)
    message(FATAL_ERROR "Failed to find perf")
  endif()

  set(MAIN_COMPILER_FLAGS "${MAIN_COMPILER_FLAGS} -g -O2 -DNDEBUG")

  set(PROFILE_TARGET "profile")
  set(PROFILE_OUTDIR "${CMAKE_BINARY_DIR}/perf")

  file(MAKE_DIRECTORY "${PROFILE_OUTDIR}")

  file(GLOB PROFILE_SOURCES "${PROFILE_SRC_DIR}/*_profile.cc")
  foreach(PROFILE_SOURCE "${PROFILE_SOURCES}")
    get_filename_component(PROFILE_PROG "${PROFILE_SOURCE}" NAME)
    string(REPLACE "_profile.cc" "" PROFILE_PROG "${PROFILE_PROG}")

    add_custom_target("${PROFILE_TARGET}_${PROFILE_PROG}"
      COMMAND ${PERF_PATH} record
              -e cpu-clock,page-faults "${PROFILE_BIN_DIR}/${PROFILE_PROG}"
      COMMAND ${PERF_PATH} annotate "--dsos=${PROFILE_PROG}"
              --stdio > "${PROFILE_OUTDIR}/${PROFILE_PROG}.annotate"
      COMMAND ${CMAKE_COMMAND} -E remove "perf.data"
    )
  endforeach()
endif()


################################################################################
# Add Subdirectories
################################################################################

add_subdirectory("${SRC_DIR}")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  enable_testing()
  add_subdirectory("${TEST_SRC_DIR}")
elseif(CMAKE_BUILD_TYPE STREQUAL "Profile")
  add_external_project(NAME "permlib"
                       CONFIGFILE "${PERMLIB_CONFIG}"
                       WORKDIR "${PERMLIB_DOWNLOAD_DIR}")

  add_subdirectory("${PROFILE_SRC_DIR}")
endif()

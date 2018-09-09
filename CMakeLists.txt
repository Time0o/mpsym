cmake_minimum_required(VERSION 3.6)

project("MPSoC Symmetry Reduction Library" "CXX")
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

set(DOC_IN_DIR "${PROJECT_SOURCE_DIR}/doc")
set(DOC_OUT_DIR "${CMAKE_BINARY_DIR}/doc")
set(DOC_PUBLISH_DIR "${PROJECT_SOURCE_DIR}/docs")
set(DOC_DOXYGEN_IN_DIR "${DOC_IN_DIR}/doxygen")
set(DOC_DOXYGEN_OUT_DIR "${DOC_OUT_DIR}/doxygen")
set(DOC_IN_SPHINX_DIR "${DOC_IN_DIR}/sphinx")
set(DOC_OUT_SPHINX_DIR "${DOC_OUT_DIR}/sphinx")


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
# Documentation
################################################################################

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  find_package(Doxygen)

  if(DOXYGEN_FOUND)
    set(DOXYFILE_IN "${DOC_DOXYGEN_IN_DIR}/Doxyfile")
    set(DOXYFILE_OUT "${DOC_DOXYGEN_OUT_DIR}/Doxyfile")

    configure_file("${DOXYFILE_IN}" "${DOXYFILE_OUT}" @ONLY)

    message(STATUS "Building documentation")

    add_custom_target(doxygen ALL
      COMMAND "${DOXYGEN_EXECUTABLE}" "${DOXYFILE_OUT}"
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
      VERBATIM
    )

    find_program(SPHINX_EXECUTABLE "sphinx-build")
    if(SPHINX_EXECUTABLE STREQUAL "sphinx-build-NOTFOUND")
      message(WARNING
              "Failed to find sphinx, sphinx documentation will not be generated")
    else()
      file(COPY "${DOC_IN_SPHINX_DIR}/conf/conf.py"
           DESTINATION "${DOC_IN_SPHINX_DIR}")

      string(CONCAT BREATHE_EXT
             "\n\n# Set up breathe\n"
             "import sys\n"
             "sys.path.append('${BREATHE_SRC_DIR}')\n"
             "extensions.append('breathe')\n"
             "breathe_projects = {\n"
             "  '${PROJECT_NAME}' : '${DOC_DOXYGEN_OUT_DIR}/xml'\n"
             "}\n"
             "breathe_default_project = '${PROJECT_NAME}'\n"
             "breathe_domain_by_extension = {'h' : 'cpp'}")

      file(APPEND "${DOC_IN_SPHINX_DIR}/conf.py" "${BREATHE_EXT}")

      add_custom_target(sphinx
        COMMAND "${SPHINX_EXECUTABLE}" -b html
                "${DOC_IN_SPHINX_DIR}" "${DOC_OUT_SPHINX_DIR}"
        VERBATIM
        DEPENDS doxygen
      )

      add_custom_target(sphinx_publish
        COMMAND "${SPHINX_EXECUTABLE}" -b html
                "${DOC_IN_SPHINX_DIR}" "${DOC_PUBLISH_DIR}"
        VERBATIM
        DEPENDS doxygen
      )
    endif()
  else()
    message(WARNING
            "Failed to find doxygen, documentation will not be generated")
  endif()
endif()


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

if(CMAKE_BUILD_TYPE STREQUAL "Profile")
  find_program(PERF_PATH perf)

  if(NOT PERF_PATH)
    message(FATAL_ERROR "Failed to find perf")
  endif()

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

enable_testing()
add_subdirectory("${TEST_SRC_DIR}")

add_subdirectory("${PROFILE_SRC_DIR}")

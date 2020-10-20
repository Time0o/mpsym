cmake_minimum_required(VERSION 3.6)

include_directories("${MPSYM_INCLUDE_DIR}"
                    "${Boost_INCLUDE_DIRS}")

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)

pybind11_add_module("${PYTHON_MODULE}" "pympsym.cpp")

target_link_libraries("${PYTHON_MODULE}"
                      PRIVATE "${MPSYM_LIB}"
                      PRIVATE nlohmann_json::nlohmann_json)

if(PYTHON_LIBRARY_OUTPUT_DIRECTORY)
  set_target_properties("${PYTHON_MODULE}" PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${PYTHON_LIBRARY_OUTPUT_DIRECTORY}")
endif()

if(NOT PYTHON_NO_SETUP_PY AND CMAKE_BUILD_TYPE STREQUAL "${CMAKE_BUILD_TYPE_RELEASE}")
  configure_file("${PYTHON_SETUP_PY_IN}" "${PYTHON_SETUP_PY_OUT}")

  install(CODE "message(STATUS \"Installing Python module\")
                execute_process(COMMAND ${PIP_EXECUTABLE} install --force .
                                WORKING_DIRECTORY ${PYTHON_BIN_DIR}
                                OUTPUT_QUIET
                                ERROR_QUIET)")
endif()

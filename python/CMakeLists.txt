cmake_minimum_required(VERSION 3.6)

include_directories("${MPSYM_INCLUDE_DIR}")

pybind11_add_module("${PYTHON_MODULE}" "pympsym.cpp")

target_link_libraries("${PYTHON_MODULE}" PRIVATE "${MPSYM_LIB}")

if(NOT PYTHON_PYPI_BUILD AND CMAKE_BUILD_TYPE STREQUAL "${CMAKE_BUILD_TYPE_RELEASE}")
  configure_file("${PYTHON_SETUP_PY_IN}" "${PYTHON_SETUP_PY_OUT}")

  install(CODE "message(STATUS \"Installing Python module\")
                execute_process(COMMAND ${PIP_EXECUTABLE} install --force .
                                WORKING_DIRECTORY ${PYTHON_BIN_DIR}
                                OUTPUT_QUIET
                                ERROR_QUIET)")
endif()

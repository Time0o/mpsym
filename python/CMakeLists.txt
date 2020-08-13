cmake_minimum_required(VERSION 3.6)

include_directories("${MPSYM_INCLUDE_DIR}")

pybind11_add_module("${PYTHON_MODULE}" "pympsym.cpp")

target_link_libraries("${PYTHON_MODULE}" PRIVATE "${MPSYM_LIB}")

if(CMAKE_BUILD_TYPE STREQUAL "${CMAKE_BUILD_TYPE_RELEASE}")
  find_program(PIP_PATH pip)
  if (NOT PIP_PATH)
    message("Failed to find pip")
    message("Python module will not be installed")
  else()
    configure_file("${PYTHON_SETUP_PY_IN}" "${PYTHON_SETUP_PY_OUT}")

    install(CODE "message(STATUS \"Installing ${PYMPSYM} Python module\")
                  execute_process(COMMAND ${PIP_PATH} install --force .\
                                  WORKING_DIRECTORY ${PYTHON_BIN_DIR})
                                  OUTPUT_QUIET
                                  ERROR_QUIET)")
  endif()
endif()

cmake_minimum_required(VERSION 3.6)

include_directories("${MPSYM_INCLUDE_DIR}")

pybind11_add_module("${PYTHON_MODULE}" "pympsym.cpp")

target_link_libraries("${PYTHON_MODULE}" PRIVATE "${MPSYM_LIB}")

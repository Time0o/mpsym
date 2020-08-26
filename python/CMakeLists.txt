cmake_minimum_required(VERSION 3.6)

include_directories("${MPSYM_INCLUDE_DIR}")

pybind11_add_module(pympsym "pympsym.cpp")

target_link_libraries(pympsym PRIVATE "${MPSYM_LIB}")

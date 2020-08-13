cmake_minimum_required(VERSION 3.6)

include_directories("${INCLUDE_DIR}")

pybind11_add_module("${MPSYM_PYTHON}_" "mpsym.cpp")

target_link_libraries("${MPSYM_PYTHON}_" PRIVATE "${MPSYM_LIB}")

set_default_target_properties(TARGET "${MPSYM_PYTHON}_")

# Avoid target name clash with static library
set_target_properties("${MPSYM_PYTHON}_" PROPERTIES
  OUTPUT_NAME "${MPSYM_PYTHON}"
)

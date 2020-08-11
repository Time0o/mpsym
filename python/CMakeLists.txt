cmake_minimum_required(VERSION 3.6)

pybind11_add_module("${MPSYM_PYTHON}_" "mpsym.cpp")

set_default_target_properties(TARGET "${MPSYM_PYTHON}_")

# Avoid target name clash with static library
set_target_properties("${MPSYM_PYTHON}_" PROPERTIES
  OUTPUT_NAME "${MPSYM_PYTHON}"
)

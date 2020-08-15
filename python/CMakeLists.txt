cmake_minimum_required(VERSION 3.6)

pybind11_add_module("${PYTHON_MODULE}" "pympsym.cpp")

configure_file("${PYTHON_SETUP_PY_IN}" "${PYTHON_SETUP_PY_OUT}")

cmake_minimum_required(VERSION 3.6)

include_directories("${MPSYM_INCLUDE_DIR}"
                    "${Boost_INCLUDE_DIRS}")

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)

pybind11_add_module("${PYTHON_MODULE}" "_mpsym.cpp")

target_link_libraries("${PYTHON_MODULE}"
                      PRIVATE "${MPSYM_LIB}"
                      PRIVATE nlohmann_json::nlohmann_json)

if(PYTHON_LIBRARY_OUTPUT_DIRECTORY)
  set_target_properties("${PYTHON_MODULE}" PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${PYTHON_LIBRARY_OUTPUT_DIRECTORY}")
endif()

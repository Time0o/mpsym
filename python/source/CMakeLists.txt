cmake_minimum_required(VERSION 3.6)

include_directories("${MPSYM_INCLUDE_DIR}"
                    "${Boost_INCLUDE_DIRS}")

if("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
  link_directories("${BOOST_LIBRARY_DIRS}"
                   "/usr/local/opt/icu4c/lib")
else()
  link_directories("${BOOST_LIBRARY_DIRS}")
endif()

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)

pybind11_add_module("${PYTHON_MODULE_INTERNAL}" "_mpsym.cpp")

target_link_libraries("${PYTHON_MODULE_INTERNAL}"
                      PRIVATE "${MPSYM_LIB}"
                      PRIVATE nlohmann_json::nlohmann_json)

if(NOT PYTHON_MODULE_INTERNAL_OUT_DIR)
  set(PYTHON_MODULE_INTERNAL_OUT_DIR "${PYTHON_MODULE_DIR}")
endif()

set_target_properties("${PYTHON_MODULE_INTERNAL}" PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY "${PYTHON_MODULE_INTERNAL_OUT_DIR}")

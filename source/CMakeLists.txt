cmake_minimum_required(VERSION 3.0)

include_directories("${PROJECT_SOURCE_DIR}/include")

add_library(cgtl "perm.cc")
set_target_properties(cgtl PROPERTIES LINKER_LANGUAGE CXX)

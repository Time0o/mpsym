cmake_minimum_required(VERSION 3.6)

include_directories("${INCLUDE_DIR}")
include_directories("${NAUTY_WORKDIR}")

add_library("${MAIN_LIB}"
            "arch_graph.cc"
            "bsgs.cc"
            "dbg.cc"
            "perm.cc"
            "perm_group.cc"
            "perm_word.cc"
            "pr_randomizer.cc"
            "schreier_sims.cc")

target_link_libraries("${MAIN_LIB}" "${NAUTY_WORKDIR}/nauty.a")

set(CGTL_COMPILER_FLAGS "-Wall -pedantic -Wextra")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CGTL_COMPILER_FLAGS "${CGTL_COMPILER_FLAGS} -g -O0 --coverage")
endif()

set_target_properties("${MAIN_LIB}" PROPERTIES
  LINKER_LANGUAGE CXX
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS "${CGTL_COMPILER_FLAGS}"
)

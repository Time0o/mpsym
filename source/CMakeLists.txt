cmake_minimum_required(VERSION 3.6)

include_directories("${INCLUDE_DIR}"
                    "${Boost_INCLUDE_DIRS}"
                    "${LUA_INCLUDE_DIR}"
                    "${NAUTY_WORKDIR}")

add_library("${MAIN_LIB}"
            "arch_graph.cc"
            "block_system.cc"
            "bsgs.cc"
            "dbg.cc"
            "eemp.cc"
            "partial_perm.cc"
            "perm.cc"
            "perm_group.cc"
            "perm_word.cc"
            "pr_randomizer.cc"
            "schreier_sims.cc")

target_link_libraries("${MAIN_LIB}"
                      "${Boost_LIBRARIES}"
                      "${LUA_LIBRARIES}"
                      "${NAUTY_WORKDIR}/nauty.a")

set(MAIN_COMPILER_FLAGS "-Wall -pedantic -Wextra")
set(MAIN_LINK_FLAGS "")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_link_libraries("${MAIN_LIB}" gcov)

  set(MAIN_COMPILER_FLAGS "${MAIN_COMPILER_FLAGS} -g -O0 --coverage")
  set(MAIN_LINK_FLAGS "${MAIN_LINK_FLAGS} --coverage")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Profile")
  set(MAIN_COMPILER_FLAGS "${MAIN_COMPILER_FLAGS} -g -O2")
endif()

set_target_properties("${MAIN_LIB}" PROPERTIES
  LINKER_LANGUAGE CXX
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS "${MAIN_COMPILER_FLAGS}"
  LINK_FLAGS "${MAIN_LINK_FLAGS}"
)

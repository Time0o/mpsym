cmake_minimum_required(VERSION 3.6)

include_directories("${INCLUDE_DIR}"
                    "${Boost_INCLUDE_DIR}"
                    "${LUA_INCLUDE_DIR}"
                    "${NAUTY_WORKDIR}")

add_library("${MAIN_LIB}"
            "arch_graph.cc"
            "block_system.cc"
            "bsgs.cc"
            "bsgs_schreier_sims.cc"
            "bsgs_solve.cc"
            "dbg.cc"
            "eemp.cc"
            "orbits.cc"
            "partial_perm.cc"
            "partial_perm_inverse_semigroup.cc"
            "perm.cc"
            "perm_group.cc"
            "perm_group_disjoint_decomp.cc"
            "perm_group_wreath_decomp.cc"
            "perm_word.cc"
            "pr_randomizer.cc"
            "schreier_structure.cc"
            "timer.cc")

target_link_libraries("${MAIN_LIB}"
                      "${Boost_LIBRARIES}"
                      "${LUA_LIBRARIES}"
                      "${NAUTY_WORKDIR}/nauty.a")

if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND DO_COVERAGE)
  target_link_libraries("${MAIN_LIB}" gcov)
endif()

set_default_target_properties(TARGET "${MAIN_LIB}")

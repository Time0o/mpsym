cmake_minimum_required(VERSION 3.6)

include_directories("${INCLUDE_DIR}"
                    "${Boost_INCLUDE_DIR}"
                    "${LUA_INCLUDE_DIR}"
                    "${NAUTY_WORKDIR}")

add_library("${MPSYM_LIB}"
            "arch_graph.cpp"
            "arch_graph_nauty.cpp"
            "arch_graph_cluster.cpp"
            "arch_graph_system.cpp"
            "arch_graph_system_lua.cpp"
            "arch_uniform_super_graph.cpp"
            "block_system.cpp"
            "bsgs.cpp"
            "bsgs_base_change.cpp"
            "bsgs_reduce_gens.cpp"
            "bsgs_schreier_sims.cpp"
            "bsgs_solve.cpp"
            "dbg.cpp"
            "eemp.cpp"
            "explicit_transversals.cpp"
            "orbits.cpp"
            "partial_perm.cpp"
            "partial_perm_inverse_semigroup.cpp"
            "perm.cpp"
            "perm_group.cpp"
            "perm_group_disjoint_decomp.cpp"
            "perm_group_wreath_decomp.cpp"
            "perm_set.cpp"
            "pr_randomizer.cpp"
            "schreier_tree.cpp"
            "timer.cpp")

target_link_libraries("${MPSYM_LIB}"
                      "${Boost_LIBRARIES}"
                      "${LUA_LIBRARIES}"
                      "${NAUTY_WORKDIR}/nauty.a")

if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND DO_COVERAGE)
  target_link_libraries("${MPSYM_LIB}" gcov)
endif()

set_default_target_properties(TARGET "${MPSYM_LIB}")

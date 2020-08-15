cmake_minimum_required(VERSION 3.6)

include_directories("${MPSYM_INCLUDE_DIR}"
                    "${Boost_INCLUDE_DIR}"
                    "${LUA_INCLUDE_DIR}"
                    "${NAUTY_WORK_DIR}")

add_library("${MPSYM_LIB}" SHARED
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
                      "${NAUTY_LIB}")

# Coverage
if(DO_COVERAGE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")

  target_link_libraries("${MPSYM_LIB}" "${GCOV_LIB}")
endif()

# Installation
if(CMAKE_BUILD_TYPE STREQUAL "${CMAKE_BUILD_TYPE_RELEASE}")
  file(GLOB PUBLIC_HEADERS "${MPSYM_INCLUDE_DIR}/*.hpp")

  set_target_properties("${MPSYM_LIB}" PROPERTIES
    PUBLIC_HEADER "${PUBLIC_HEADERS}"
  )

  install(TARGETS "${MPSYM_LIB}"
    LIBRARY DESTINATION "${MPSYM_INSTALL_LIB_DIR}"
    PUBLIC_HEADER DESTINATION "${MPSYM_INSTALL_INCLUDE_DIR}"
  )
endif()

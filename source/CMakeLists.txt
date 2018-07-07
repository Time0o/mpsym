cmake_minimum_required(VERSION 3.6)

include_directories("${PROJECT_SOURCE_DIR}/include")

add_library(cgtl "bsgs.cc" "dbg.cc" "perm.cc" "perm_word.cc"
                 "perm_group.cc" "pr_randomizer.cc" "schreier_sims.cc")

set_target_properties(cgtl PROPERTIES
  LINKER_LANGUAGE CXX
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS "-Wall -pedantic -Wextra"
)

cmake_minimum_required(VERSION 3.6)

include_directories("${PROJECT_SOURCE_DIR}/include")

add_library(cgtl "bsgs.cc" "dbg.cc" "perm.cc" "perm_group.cc" "schreier_sims.cc")

set_target_properties(cgtl PROPERTIES
  LINKER_LANGUAGE CXX
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS "-Wall -pedantic -Wextra"
)

#option(ENABLE_CLANG_TIDY "ENABLE_CLANG_TIDY" OFF)
#
#if(ENABLE_CLANG_TIDY)
#  find_program(CLANG_TIDY NAMES "clang-tidy")
#
#  if(CLANG_TIDY)
#    string(CONCAT CLANG_CHECKS
#      "-checks=bugprone-*,"
#      ",cert-*,"
#      "cppcoreguidelines-*,"
#      "misc-*,"
#      "modernize-*,"
#      "performance-*,"
#      "portability*"
#    )
#
#    set(EXEC_CLANG_TIDY "${CLANG_TIDY}" "${CLANG_CHECKS}")
#
#    set_target_properties(cgtl PROPERTIES CXX_CLANG_TIDY "${EXEC_CLANG_TIDY}")
#  else()
#    message(WARNING "failed to find clang-tidy installation")
#  endif()
#endif()

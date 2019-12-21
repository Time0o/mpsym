cmake_minimum_required(VERSION 3.6)

include_directories("${INCLUDE_DIR}" "${PROFILE_INCLUDE_DIR}")
include_directories(SYSTEM "${PERMLIB_INCLUDE_DIR}") # to ignore all warnings

file(GLOB PROFILE_SOURCES "*_profile.cc")

file(GLOB PROFILE_UTIL "../common/profile_*.cc")

foreach(PROFILE_SOURCE ${PROFILE_SOURCES})
  get_filename_component(PROFILE_PROG "${PROFILE_SOURCE}" NAME)
  string(REPLACE "_profile.cc" "" PROFILE_PROG "${PROFILE_PROG}")

  add_executable("${PROFILE_PROG}" "${PROFILE_SOURCE}" "${PROFILE_UTIL}")

  target_link_libraries("${PROFILE_PROG}" "${MAIN_LIB}")

  set_default_target_properties(TARGET "${PROFILE_PROG}")

  set_target_properties("${PROFILE_PROG}" PROPERTIES
    CXX_STANDARD 17
  )
endforeach()

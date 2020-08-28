cmake_minimum_required(VERSION 3.6)

include_directories("${MPSYM_INCLUDE_DIR}"
                    "${MPSYM_PROFILE_INCLUDE_DIR}"
                    "${NAUTY_WORK_DIR}")

include_directories(SYSTEM "${PERMLIB_INCLUDE_DIR}") # to ignore all warnings

file(GLOB PROFILE_SOURCES "*_profile.cpp")

file(GLOB PROFILE_UTIL "${PROFILE_COMMON_DIR}/profile_*.cpp")

foreach(PROFILE_SOURCE ${PROFILE_SOURCES})
  get_filename_component(PROFILE_PROG "${PROFILE_SOURCE}" NAME)
  string(REPLACE "_profile.cpp" "" PROFILE_PROG "${PROFILE_PROG}")

  add_executable("${PROFILE_PROG}" "${PROFILE_SOURCE}" "${PROFILE_UTIL}")

  target_link_libraries("${PROFILE_PROG}" "${MPSYM_LIB}")

  set_target_properties("${PROFILE_PROG}" PROPERTIES
    CXX_STANDARD 17
  )
endforeach()

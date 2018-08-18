cmake_minimum_required(VERSION 3.6)

include_directories("${INCLUDE_DIR}"
                    "${PROFILE_INCLUDE_DIR}")

file(GLOB PROFILE_SOURCES "*_profile.cc")
set(PROFILE_UTILITY "profile_utility.cc")

file(GLOB PROFILE_RESOURCES "${PROFILE_RESOURCE_DIR}/*")
file(COPY "${PROFILE_RESOURCES}" DESTINATION "${RESOURCE_DIR}/profile")

foreach(PROFILE_SOURCE "${PROFILE_SOURCES}")
  get_filename_component(PROFILE_PROG "${PROFILE_SOURCE}" NAME)
  string(REPLACE "_profile.cc" "" PROFILE_PROG "${PROFILE_PROG}")

  add_executable("${PROFILE_PROG}" "${PROFILE_SOURCE}" "${PROFILE_UTILITY}")

  target_link_libraries("${PROFILE_PROG}"
                        "${MAIN_LIB}")

  set_target_properties("${PROFILE_PROG}" PROPERTIES
    LINKER_LANGUAGE CXX
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    COMPILE_FLAGS "-Wall -pedantic -Wextra -g -O2"
  )
endforeach()
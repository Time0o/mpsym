cmake_minimum_required(VERSION 3.6)

include_directories("${INCLUDE_DIR}")

file(GLOB EXAMPLE_SOURCES "*.cc")

foreach(EXAMPLE_SOURCE ${EXAMPLE_SOURCES})
  get_filename_component(EXAMPLE_PROG "${EXAMPLE_SOURCE}" NAME)
  string(REPLACE ".cc" "" EXAMPLE_PROG "${EXAMPLE_PROG}")

  add_executable("${EXAMPLE_PROG}" "${EXAMPLE_SOURCE}")

  target_link_libraries("${EXAMPLE_PROG}" "${MAIN_LIB}")

  set_default_target_properties(TARGET "${EXAMPLE_PROG}")
endforeach()

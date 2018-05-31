cmake_minimum_required(VERSION 3.0)

include_directories("${PROJECT_SOURCE_DIR}/include")

file(GLOB testsources "*_test.cc")
foreach(testsource ${testsources})
  get_filename_component(test "${testsource}" NAME)
  string(REPLACE ".cc" "" test "${test}")

  add_executable("${test}" "${testsource}")
  target_link_libraries("${test}" gmock_main)
  target_link_libraries("${test}" cgtl)
  add_test(NAME "${test}" COMMAND "${test}")
endforeach()

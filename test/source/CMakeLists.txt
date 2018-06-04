cmake_minimum_required(VERSION 3.6)

include_directories("${PROJECT_SOURCE_DIR}/include"
                    "${PROJECT_SOURCE_DIR}/test/include")

set(TEST_UTILITY "test_utility.cc")

file(GLOB testsources "*_test.cc")
foreach(testsource ${testsources})
  get_filename_component(test "${testsource}" NAME)
  string(REPLACE ".cc" "" test "${test}")

  add_executable("${test}" "${testsource}" "${TEST_UTILITY}")
  target_link_libraries("${test}" gmock_main)
  target_link_libraries("${test}" cgtl)
  add_test(NAME "${test}" COMMAND "${test}")
endforeach()

cmake_minimum_required(VERSION 3.6)

include_directories("${INCLUDE_DIR}"
                    "${TEST_INCLUDE_DIR}"
                    "${RESOURCE_INCLUDE_DIR}")

file(GLOB TEST_SOURCES "*_test.cc")
set(TEST_UTILITY "test_utility.cc")

foreach(TEST_SOURCE ${TEST_SOURCES})
  get_filename_component(TEST_PROG "${TEST_SOURCE}" NAME)
  string(REPLACE "_test.cc" "" TEST_PROG "${TEST_PROG}")

  add_test(NAME "${TEST_PROG}"
           COMMAND "${TEST_PROG}"
           WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

  add_executable("${TEST_PROG}" "${TEST_SOURCE}" "${TEST_UTILITY}")

  target_link_libraries("${TEST_PROG}" "${MAIN_LIB}" gmock_main)

  set_default_target_properties(TARGET "${TEST_PROG}")
endforeach()

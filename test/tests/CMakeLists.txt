cmake_minimum_required(VERSION 3.6)

include_directories("${INCLUDE_DIR}"
                    "${TEST_INCLUDE_DIR}"
                    "${TEST_SRC_DIR}"
                    "${NAUTY_WORKDIR}")

file(GLOB TEST_SOURCES "*_test.cpp")

set(TEST_UTILITY "../source/test_utility.cpp")

foreach(TEST_SOURCE ${TEST_SOURCES})
  get_filename_component(TEST_PROG "${TEST_SOURCE}" NAME)
  string(REPLACE "_test.cpp" "" TEST_PROG "${TEST_PROG}")

  add_test(NAME "${TEST_PROG}"
           COMMAND "${TEST_PROG}"
           WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

  add_executable("${TEST_PROG}" "${TEST_SOURCE}" "${TEST_UTILITY}")

  target_link_libraries("${TEST_PROG}" "${MPSYM_LIB}" gmock_main)

  set_default_target_properties(TARGET "${TEST_PROG}")
endforeach()

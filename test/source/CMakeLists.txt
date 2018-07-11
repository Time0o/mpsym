cmake_minimum_required(VERSION 3.6)

include_directories("${INCLUDE_DIR}" "${TEST_INCLUDE_DIR}")

set(TEST_UTILITY "test_utility.cc")

file(GLOB testsources "*_test.cc")
file(COPY "${TEST_RESOURCE_DIR}" DESTINATION "${CMAKE_BINARY_DIR}")

foreach(testsource ${testsources})
  get_filename_component(test "${testsource}" NAME)
  string(REPLACE ".cc" "" test "${test}")

  add_test(NAME "${test}"
           COMMAND "${test}"
           WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

  add_executable("${test}" "${testsource}" "${TEST_UTILITY}")

  target_link_libraries("${test}" "${MAIN_LIB}" gmock_main)

  set_target_properties("${test}" PROPERTIES
    LINKER_LANGUAGE CXX
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    COMPILE_FLAGS "-Wall -pedantic -Wextra"
  )
endforeach()

cmake_minimum_required(VERSION 3.6)

# Valgrind support
set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --error-exitcode=1")

include(CTest)

add_custom_target(test_memcheck
  COMMAND ${CMAKE_CTEST_COMMAND} --force-new-ctest-process --test-action memcheck
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
)

# C++
include_directories("${MPSYM_INCLUDE_DIR}"
                    "${MPSYM_TEST_INCLUDE_DIR}"
                    "${MPSYM_TEST_SRC_DIR}"
                    "${Boost_INCLUDE_DIRS}"
                    "${LUA_INCLUDE_DIR}"
                    "${NAUTY_WORK_DIR}")

file(GLOB TEST_SOURCES "*_test.cpp")

set(TEST_UTILITY "${MPSYM_TEST_SRC_DIR}/test_utility.cpp")

foreach(TEST_SOURCE ${TEST_SOURCES})
  get_filename_component(TEST_PROG "${TEST_SOURCE}" NAME)
  string(REPLACE "_test.cpp" "" TEST_PROG "${TEST_PROG}")

  add_test(NAME "${TEST_PROG}"
           COMMAND "${TEST_PROG}"
           WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

  add_executable("${TEST_PROG}" "${TEST_SOURCE}" "${TEST_UTILITY}")

  target_link_libraries("${TEST_PROG}" "${MPSYM_LIB}" "${GTEST_GMOCK_MAIN}")
endforeach()

# Python
if(PYTHON_BINDINGS)
  add_test(
    NAME python
    COMMAND ${CMAKE_COMMAND} -E env "LUA_PATH=${LUA_SRC_DIR}/?.lua;$ENV{LUA_PATH}" "${PYTHON_EXECUTABLE}" -m pytest python_test.py
    WORKING_DIRECTORY "${MPSYM_TEST_TESTS_DIR}"
  )

  set_tests_properties(
    python
    PROPERTIES ENVIRONMENT "PYTHONPATH=${PYTHON_BIN_DIR}:$ENV{PYTHONPATH}"
  )
endif()

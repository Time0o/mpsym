file(GLOB TOOL_SOURCES "**.cc")

foreach(TOOL_SOURCE ${TOOL_SOURCES})
  get_filename_component(TOOL_PROG "${TOOL_SOURCE}" NAME)
  string(REPLACE ".cc" "" TOOL_PROG "${TOOL_PROG}")

  set(EXECUTABLE_OUTPUT_PATH "${TOOL_BIN_DIR}")
  add_executable("${TOOL_PROG}" "${TOOL_SOURCE}")
endforeach()

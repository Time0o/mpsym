function(message)
    if (NOT MESSAGE_QUIET)
        _message(${ARGN})
    endif()
endfunction()

function(add_external_project)
  set(options NONE)
  set(args NAME CONFIG_FILE WORK_DIR SRC_DIR BIN_DIR)
  cmake_parse_arguments(Arg "${options}" "${args}" "" "${ARGN}")

  message(STATUS "Processing ${Arg_CONFIG_FILE}")

  set(MESSAGE_QUIET ON)

  configure_file("${Arg_CONFIG_FILE}" "${Arg_WORK_DIR}/CMakeLists.txt")

  execute_process(
    COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY "${Arg_WORK_DIR}"
    OUTPUT_QUIET
    ERROR_QUIET
  )

  if(result)
    unset(MESSAGE_QUIET)
    message(FATAL_ERROR "CMake step for ${Arg_NAME} failed")
  endif()

  execute_process(
    COMMAND "${CMAKE_COMMAND}" --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY "${Arg_WORK_DIR}"
    OUTPUT_QUIET
    ERROR_QUIET
  )

  if(result)
    unset(MESSAGE_QUIET)
    message(FATAL_ERROR "Build step for ${Arg_NAME} failed")
  endif()

  if (DEFINED Arg_SRC_DIR AND DEFINED Arg_BIN_DIR)
    add_subdirectory("${Arg_SRC_DIR}" "${Arg_BIN_DIR}" EXCLUDE_FROM_ALL)
  endif()

  unset(MESSAGE_QUIET)
endfunction()

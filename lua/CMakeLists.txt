cmake_minimum_required(VERSION 3.6)

if(NOT LUA_NO_ROCK AND NOT LUA_EMBED AND CMAKE_BUILD_TYPE STREQUAL "${CMAKE_BUILD_TYPE_RELEASE}")
  configure_file("${LUA_ROCKSPEC_IN}" "${LUA_ROCKSPEC_OUT}")

  install(CODE "message(STATUS \"Installing Lua rock\")
                execute_process(COMMAND ${LUAROCKS_EXECUTABLE} make ${LUA_ROCKSPEC_OUT}
                                OUTPUT_QUIET
                                ERROR_QUIET)")
endif()

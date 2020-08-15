cmake_minimum_required(VERSION 3.6)

add_library("${NAUTY_LIB}" SHARED
            "${NAUTY_WORK_DIR}/gtnauty.c"
            "${NAUTY_WORK_DIR}/gtools.c"
            "${NAUTY_WORK_DIR}/gutil1.c"
            "${NAUTY_WORK_DIR}/gutil2.c"
            "${NAUTY_WORK_DIR}/naugraph.c"
            "${NAUTY_WORK_DIR}/naugroup.c"
            "${NAUTY_WORK_DIR}/naurng.c"
            "${NAUTY_WORK_DIR}/nausparse.c"
            "${NAUTY_WORK_DIR}/nautil.c"
            "${NAUTY_WORK_DIR}/nautinv.c"
            "${NAUTY_WORK_DIR}/naututil.c"
            "${NAUTY_WORK_DIR}/nauty.c"
            "${NAUTY_WORK_DIR}/schreier.c"
            "${NAUTY_WORK_DIR}/traces.c")

if (CMAKE_BUILD_TYPE STREQUAL "${CMAKE_BUILD_TYPE_RELEASE}")
  install(TARGETS "${NAUTY_LIB}"
    LIBRARY DESTINATION "${NAUTY_INSTALL_LIB_DIR}"
  )
endif()

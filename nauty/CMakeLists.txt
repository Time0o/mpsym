cmake_minimum_required(VERSION 3.6)

add_library("${NAUTY_LIB}" SHARED
            "${NAUTY_WORKDIR}/gtnauty.c"
            "${NAUTY_WORKDIR}/gtools.c"
            "${NAUTY_WORKDIR}/gutil1.c"
            "${NAUTY_WORKDIR}/gutil2.c"
            "${NAUTY_WORKDIR}/naugraph.c"
            "${NAUTY_WORKDIR}/naugroup.c"
            "${NAUTY_WORKDIR}/naurng.c"
            "${NAUTY_WORKDIR}/nausparse.c"
            "${NAUTY_WORKDIR}/nautil.c"
            "${NAUTY_WORKDIR}/nautinv.c"
            "${NAUTY_WORKDIR}/naututil.c"
            "${NAUTY_WORKDIR}/nauty.c"
            "${NAUTY_WORKDIR}/schreier.c"
            "${NAUTY_WORKDIR}/traces.c")

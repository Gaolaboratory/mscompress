set(ZLIB_SOURCE_DIR ${VENDOR_DIR}/zlib)
# Add a custom target to run nmake and build zlib
# Set ZLIB_EXTRA_CFLAGS before including this file to pass additional flags.
if(MINGW)
# Cross build with GCC (e.g. x86_64-w64-mingw32-gcc from Linux/WSL). Makefile.gcc
# builds libz.a in-tree, same as the Unix branch. HAS_SSE2 for the same reason
# as the nmake branch below (the fork's deflate gates hash-table sliding on it).
add_custom_target(
zlib_build
COMMAND make -C ${ZLIB_SOURCE_DIR} -f win32/Makefile.gcc
CC=${CMAKE_C_COMPILER} AR=${CMAKE_AR} RANLIB=${CMAKE_RANLIB}
"CFLAGS=-O3 -DHAS_SSE2 -include emmintrin.h ${ZLIB_EXTRA_CFLAGS}" libz.a
WORKING_DIRECTORY ${ZLIB_SOURCE_DIR}
BYPRODUCTS ${ZLIB_SOURCE_DIR}/libz.a
VERBATIM
)
elseif(WIN32)
    # The cloudflare fork's Makefile.msc omits HAS_SSE2, without which deflate
    # never slides its hash tables and corrupts any input larger than ~64 KB
    # (configure sets it on Unix). deflate.c only includes intrinsics headers
    # under HAS_SSE42, so force-include emmintrin.h for __m128i.
    # Build out-of-source (TOP points nmake at the sources) so projects that
    # pass different ZLIB_EXTRA_CFLAGS (e.g. node's Z_PREFIX vs. the CLI)
    # don't reuse each other's stale zlib.lib from the shared vendor tree.
    set(ZLIB_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/zlib-build)
    file(MAKE_DIRECTORY ${ZLIB_BUILD_DIR})
    add_custom_target(
        zlib_build
        COMMAND nmake -f ${ZLIB_SOURCE_DIR}/win32/Makefile.msc zlib.lib "TOP=${ZLIB_SOURCE_DIR}" "LOC=-DHAS_SSE2 -FIemmintrin.h ${ZLIB_EXTRA_CFLAGS}"
        WORKING_DIRECTORY ${ZLIB_BUILD_DIR}
        BYPRODUCTS ${ZLIB_BUILD_DIR}/zlib.lib
        VERBATIM
    )
elseif(APPLE)
    # On macOS, define fdopen to prevent zlib's macro redefinition conflict
    add_custom_target(
        zlib_build
        COMMAND env "CFLAGS=-O3 -Dfdopen=fdopen ${ZLIB_EXTRA_CFLAGS}" ${ZLIB_SOURCE_DIR}/configure
        COMMAND make -C ${ZLIB_SOURCE_DIR}
        WORKING_DIRECTORY ${ZLIB_SOURCE_DIR}
    )
else()
    add_custom_target(
        zlib_build
        COMMAND env "CFLAGS=-O3 -fPIC ${ZLIB_EXTRA_CFLAGS}" ${ZLIB_SOURCE_DIR}/configure
        COMMAND make -C ${ZLIB_SOURCE_DIR}
        WORKING_DIRECTORY ${ZLIB_SOURCE_DIR}
    )
endif()
# Set the location of the built zlib library
if(MSVC)
    set(ZLIB_LIBRARY "${ZLIB_BUILD_DIR}/zlib.lib")
elseif(MINGW)
set(ZLIB_LIBRARY "${ZLIB_SOURCE_DIR}/libz.a")
else()
    set(ZLIB_LIBRARY "${ZLIB_SOURCE_DIR}/libz.a")
endif()

include_directories(${ZLIB_SOURCE_DIR})
add_dependencies(mscompress zlib_build)

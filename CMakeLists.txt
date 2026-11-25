cmake_minimum_required(VERSION 3.15)

project(EVG C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(Threads REQUIRED)
find_library(DEFLATE_LIBRARY NAMES deflate)
find_path(DEFLATE_INCLUDE_DIR NAMES libdeflate.h)
find_program(GNU_MAKE NAMES gmake make)
if(NOT GNU_MAKE)
    message(FATAL_ERROR "GNU make is required to build external/htslib")
endif()

# ------ htslib (external/htslib) -> static libhts.a; no libcurl, libdeflate when found ------
set(HTS_SRC "${CMAKE_CURRENT_SOURCE_DIR}/external/htslib")
if(NOT EXISTS "${HTS_SRC}/bgzf.c" OR NOT EXISTS "${HTS_SRC}/htscodecs/configure.ac")
    message(FATAL_ERROR "external/htslib missing. Run: git submodule update --init --recursive")
endif()

# Built in a copy of the submodule so the source tree stays clean; the directory name follows the compression library.
if(DEFLATE_LIBRARY AND DEFLATE_INCLUDE_DIR)
    set(HTS_MODE "libdeflate")
else()
    set(HTS_MODE "zlib")
endif()
message(STATUS "htslib compression: ${HTS_MODE}")
set(HTS_BUILD "${CMAKE_BINARY_DIR}/external_build/htslib-${HTS_MODE}")
set(HTS_LIB "${HTS_BUILD}/libhts.a")

set(HTS_INCLUDE_DIRS ${ZLIB_INCLUDE_DIRS} ${BZIP2_INCLUDE_DIR} ${LIBLZMA_INCLUDE_DIRS})
set(HTS_LINK_LIBRARIES ZLIB::ZLIB BZip2::BZip2 ${LIBLZMA_LIBRARIES} Threads::Threads)
set(HTS_DEFLATE_COMMAND)
if(HTS_MODE STREQUAL "libdeflate")
    list(APPEND HTS_INCLUDE_DIRS ${DEFLATE_INCLUDE_DIR})
    list(INSERT HTS_LINK_LIBRARIES 0 ${DEFLATE_LIBRARY})
    set(HTS_DEFLATE_COMMAND COMMAND /bin/sh -c "echo '#define HAVE_LIBDEFLATE 1' >> '${HTS_BUILD}/config.h'")
endif()

set(HTS_CPPFLAGS "")
foreach(dir ${HTS_INCLUDE_DIRS})
    string(APPEND HTS_CPPFLAGS " -I${dir}")
endforeach()
string(STRIP "${HTS_CPPFLAGS}" HTS_CPPFLAGS)

# htscodecs/htscodecs/version.h normally comes from git; write it from configure.ac (rewritten only when the version changes).
file(STRINGS "${HTS_SRC}/htscodecs/configure.ac" HTSCODECS_AC_INIT LIMIT_COUNT 1 REGEX "^AC_INIT\\(htscodecs,")
string(REGEX MATCH "[0-9]+(\\.[0-9]+)+" HTSCODECS_VERSION "${HTSCODECS_AC_INIT}")
if(NOT HTSCODECS_VERSION)
    message(FATAL_ERROR "Failed to read the htscodecs version from ${HTS_SRC}/htscodecs/configure.ac")
endif()
set(HTSCODECS_VERSION_HEADER "${CMAKE_BINARY_DIR}/generated/htscodecs_version.h")
set(HTSCODECS_VERSION_TEXT "#define HTSCODECS_VERSION_TEXT \"${HTSCODECS_VERSION}\"\n")
set(HTSCODECS_VERSION_OLD "")
if(EXISTS "${HTSCODECS_VERSION_HEADER}")
    file(READ "${HTSCODECS_VERSION_HEADER}" HTSCODECS_VERSION_OLD)
endif()
if(NOT HTSCODECS_VERSION_OLD STREQUAL HTSCODECS_VERSION_TEXT)
    file(WRITE "${HTSCODECS_VERSION_HEADER}" "${HTSCODECS_VERSION_TEXT}")
endif()
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${HTS_SRC}/htscodecs/configure.ac")

# Compile inputs of libhts.a, so a submodule update that changes only sources rebuilds it (re-globbed on every build).
file(GLOB HTS_SOURCES CONFIGURE_DEPENDS
    "${HTS_SRC}/*.c" "${HTS_SRC}/*.h" "${HTS_SRC}/*.mk" "${HTS_SRC}/*.sh" "${HTS_SRC}/Makefile"
    "${HTS_SRC}/htslib/*.h" "${HTS_SRC}/cram/*.c" "${HTS_SRC}/cram/*.h" "${HTS_SRC}/os/*.c" "${HTS_SRC}/os/*.h"
    "${HTS_SRC}/htscodecs/htscodecs/*.c" "${HTS_SRC}/htscodecs/htscodecs/*.h")

set(HTS_MAKE ${CMAKE_COMMAND} -E env MAKEFLAGS=s ${GNU_MAKE} -C "${HTS_BUILD}" "CC=${CMAKE_C_COMPILER}" "CPPFLAGS=${HTS_CPPFLAGS}")

add_custom_command(
    OUTPUT "${HTS_LIB}"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${HTS_BUILD}"
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${HTS_SRC}" "${HTS_BUILD}"
    COMMAND ${CMAKE_COMMAND} -E remove -f "${HTS_BUILD}/.git" "${HTS_BUILD}/htscodecs/.git"
    COMMAND ${HTS_MAKE} distclean
    COMMAND ${CMAKE_COMMAND} -E copy "${HTSCODECS_VERSION_HEADER}" "${HTS_BUILD}/htscodecs/htscodecs/version.h"
    COMMAND ${HTS_MAKE} config.h
    COMMAND /bin/sh -c "sed -i '/^#define HAVE_LIBCURL 1$/d' '${HTS_BUILD}/config.h'"
    ${HTS_DEFLATE_COMMAND}
    COMMAND ${HTS_MAKE} -j 4 "NONCONFIGURE_OBJS=" "htslib_default_libs=-lz -lm -lbz2 -llzma" libhts.a
    DEPENDS "${HTSCODECS_VERSION_HEADER}" ${HTS_SOURCES}
    COMMENT "Building htslib (libhts.a, ${HTS_MODE})"
    VERBATIM
)
add_custom_target(htslib_build DEPENDS "${HTS_LIB}")

add_library(hts_static STATIC IMPORTED GLOBAL)
set_target_properties(hts_static PROPERTIES
    IMPORTED_LOCATION "${HTS_LIB}"
    INTERFACE_INCLUDE_DIRECTORIES "${HTS_SRC}"
)

# ------ Sources ------
set(SRC_FILES_common
    src/hts_file.cpp
    src/save.cpp
    src/strip_split_join.cpp
)

set(SRC_FILES_graphvcf
    src/graphvcf_options.cpp
    src/vcf_open.cpp
    src/vcf_convert.cpp
    src/vcf_merge.cpp
    src/vcf_recall.cpp
    src/vcf_split.cpp
    src/vcf_breakpoint.cpp
    src/vcf_count.cpp
    src/vcf_filter.cpp
    src/delly2vcf.cpp
    src/vcf_win.cpp
)

set(SRC_FILES_fastAQ
    src/fastaq_options.cpp
    src/fastaq_count.cpp
    src/fastaq_sample.cpp
    src/fastaq_convert.cpp
)

# ------ Targets ------
add_executable(graphvcf graphvcf.cpp ${SRC_FILES_graphvcf} ${SRC_FILES_common})
add_executable(fastAQ fastAQ.cpp ${SRC_FILES_fastAQ} ${SRC_FILES_common})

foreach(target graphvcf fastAQ)
    target_include_directories(${target} PRIVATE include)
    target_compile_options(${target} PRIVATE -O3)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(${target} PRIVATE -Wall -Wextra)
    endif()
    target_link_libraries(${target} hts_static ${HTS_LINK_LIBRARIES})
    add_dependencies(${target} htslib_build)
endforeach()

# ------ Install ------
# EVG.py finds graphvcf, fastAQ, src/evg_pipeline/ and include/ProgramMetadata.hpp next to itself.
install(TARGETS graphvcf fastAQ DESTINATION bin)
install(PROGRAMS EVG.py DESTINATION bin)
install(DIRECTORY src/evg_pipeline/ DESTINATION bin/src/evg_pipeline FILES_MATCHING PATTERN "*.py" PATTERN "__pycache__" EXCLUDE)
install(FILES include/ProgramMetadata.hpp DESTINATION bin/include)
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink EVG.py ${CMAKE_CURRENT_BINARY_DIR}/EVG)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/EVG DESTINATION bin)

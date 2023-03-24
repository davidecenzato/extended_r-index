CMAKE_MINIMUM_REQUIRED(VERSION 3.15)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")

# Set a default build type if none was specified
# ------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# About this project
# ------------------------------------------------------------------------------
project(ebwt_r-index)
SET(VERSION_MAJOR "0")
SET(VERSION_MINOR "1")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

# Set environment
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install)

find_package(Git)
if(GIT_FOUND)
    message("git found: ${GIT_EXECUTABLE}")
else()
    message(WARNING "git not found. Cloning of submodules will not work.")
endif()

# Configure thirdparty
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_INCLUDEDIR "include") # This is an hack because include(GUIInstallDirs) doesn't work

add_subdirectory(thirdparty)


# Configure the compiler with the appropriate flags
# ------------------------------------------------------------------------------
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  # using Clang
  include(ConfigureCompilerClang)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  include(ConfigureCompilerGcc)
else ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	message(FATAL_ERROR "Only the compiler gcc and clang are supported")
endif()


# Targets
include_directories(${PROJECT_SOURCE_DIR})

add_executable(er-index main.cpp)
target_link_libraries(er-index malloc_count dl sdsl divsufsort divsufsort64)

add_executable(er-index64 main.cpp)
target_link_libraries(er-index64 malloc_count dl sdsl divsufsort divsufsort64)
target_compile_options(er-index64 PUBLIC "-DM64")

add_executable(genpattern genpattern.cpp)
target_link_libraries(genpattern malloc_count dl sdsl divsufsort divsufsort64)

add_executable(circpfpNT.x pfpebwt/circpfp.cpp pfpebwt/utils.c)
target_link_libraries(circpfpNT.x malloc_count z)
target_compile_options(circpfpNT.x PUBLIC "-DNOTHREADS")

add_executable(parsebwtNT.x pfpebwt/parse.cpp pfpebwt/utils.c pfpebwt/csais.cpp)
target_link_libraries(parsebwtNT.x malloc_count dl sdsl divsufsort divsufsort64)

add_executable(parsebwtNT64.x pfpebwt/parse.cpp pfpebwt/utils.c pfpebwt/csais.cpp)
target_link_libraries(parsebwtNT64.x malloc_count dl sdsl divsufsort divsufsort64)
target_compile_options(parsebwtNT64.x PUBLIC "-DP64")

add_executable(bebwtNT.x pfpebwt/ebwt.cpp pfpebwt/utils.c pfpebwt/gsa/gsacak.c)
target_link_libraries(bebwtNT.x malloc_count dl sdsl divsufsort divsufsort64)

add_executable(bebwtNTp64.x pfpebwt/ebwt.cpp pfpebwt/utils.c pfpebwt/gsa/gsacak.c)
target_link_libraries(bebwtNTp64.x malloc_count dl sdsl divsufsort divsufsort64)
target_compile_options(bebwtNTp64.x PUBLIC "-DP64")

add_executable(bebwtNTd64.x pfpebwt/ebwt.cpp pfpebwt/utils.c pfpebwt/gsa/gsacak.c)
target_link_libraries(bebwtNTd64.x malloc_count dl z sdsl divsufsort divsufsort64)
target_compile_options(bebwtNTd64.x PUBLIC "-DM64")

add_executable(bebwtNT64.x pfpebwt/ebwt.cpp pfpebwt/utils.c pfpebwt/gsa/gsacak.c)
target_link_libraries(bebwtNT64.x malloc_count dl sdsl divsufsort divsufsort64)
target_compile_options(bebwtNT64.x PUBLIC "-DM64")
target_compile_options(bebwtNT64.x PUBLIC "-DP64")


# configure_file(${PROJECT_SOURCE_DIR}/ext_r-index.py ${PROJECT_BINARY_DIR}/ext_r-index.py)
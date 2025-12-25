cmake_minimum_required(VERSION 3.20)

project(slatec-modern
    VERSION 4.2.0
    DESCRIPTION "Modernised SLATEC mathematical library - Fortran 2018"
    HOMEPAGE_URL "https://github.com/Zaneham/SLATEC"
    LANGUAGES Fortran
)

# Options
# Note: Shared library disabled by default due to duplicate symbol issues
# from include file structure. Static library works correctly.
option(SLATEC_BUILD_SHARED "Build shared library" OFF)
option(SLATEC_BUILD_STATIC "Build static library" ON)
option(SLATEC_BUILD_TESTS "Build test suite" ON)
option(SLATEC_USE_EXTERNAL_BLAS "Use external BLAS (otherwise use internal)" OFF)
option(SLATEC_USE_EXTERNAL_LAPACK "Use external LAPACK (otherwise use internal)" OFF)

# Build type
set(CMAKE_CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo" CACHE STRING "" FORCE)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# Compiler flags
# Note: Using -std=gnu (not f2018) for GNU to allow legacy mixed-kind literals
# The codebase uses patterns like 0.5E-18_DP which are GNU extensions
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fimplicit-none")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fcheck=all -fbacktrace -Wall -Wextra")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -march=native")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -stand f18")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check all -traceback -warn all")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -xHost")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Flang" OR CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2018")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
endif()

# External dependencies (optional)
if(SLATEC_USE_EXTERNAL_BLAS)
    find_package(BLAS REQUIRED)
endif()
if(SLATEC_USE_EXTERNAL_LAPACK)
    find_package(LAPACK REQUIRED)
endif()

# Source files
add_subdirectory(src/modern)

# Tests
if(SLATEC_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

# Summary
message(STATUS "")
message(STATUS "SLATEC-Modern ${PROJECT_VERSION} configuration:")
message(STATUS "  Compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Shared library: ${SLATEC_BUILD_SHARED}")
message(STATUS "  Static library: ${SLATEC_BUILD_STATIC}")
message(STATUS "  Tests: ${SLATEC_BUILD_TESTS}")
message(STATUS "")

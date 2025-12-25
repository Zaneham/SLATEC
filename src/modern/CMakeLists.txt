# SLATEC-Modern Library
#
# Builds the complete library from top-level module files.
# The .f90 files include .inc files which contain the actual implementations.
#
# NOTE: We do NOT use GLOB_RECURSE because subdirectory .f90 files are either:
#   - Included via .inc files (would cause duplicate symbols)
#   - Module definition files included via .f90 includes (would cause duplicates)

# Top-level module files that form the library
# Each module file includes its dependencies via include statements
# Order matters: blas -> lapack -> linpack -> linear
set(SLATEC_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/slatec.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/service/service.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/approximation/approximation.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/data_handling/data_handling.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/diff_integ/diff_integ.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/diff_integ_eq/diff_integ_eq.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/integ_trans/integ_trans.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/interpolation/interpolation.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/linear/blas.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/linear/lapack.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/linear/linpack.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/linear/lapack_external.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/linear/linear.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/nonlin_eq/nonlin_eq.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/optimization/optimization.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/special_functions/special_functions.f90"
)

# Static library
if(SLATEC_BUILD_STATIC)
    add_library(slatec_static STATIC ${SLATEC_SOURCES})
    set_target_properties(slatec_static PROPERTIES
        OUTPUT_NAME slatec
        VERSION ${PROJECT_VERSION}
        POSITION_INDEPENDENT_CODE ON
        Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"
    )
    target_include_directories(slatec_static PUBLIC
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
        $<INSTALL_INTERFACE:include/slatec>
    )
    # Alias for use in same build tree
    add_library(slatec::static ALIAS slatec_static)
endif()

# Shared library - must depend on static if both built to avoid race condition
if(SLATEC_BUILD_SHARED)
    add_library(slatec_shared SHARED ${SLATEC_SOURCES})
    set_target_properties(slatec_shared PROPERTIES
        OUTPUT_NAME slatec
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
        Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"
    )
    target_include_directories(slatec_shared PUBLIC
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
        $<INSTALL_INTERFACE:include/slatec>
    )
    add_library(slatec::shared ALIAS slatec_shared)
    # Serialize builds to prevent module file race condition
    if(SLATEC_BUILD_STATIC)
        add_dependencies(slatec_shared slatec_static)
    endif()
endif()

# Default target alias (prefer static)
if(SLATEC_BUILD_STATIC)
    add_library(slatec::slatec ALIAS slatec_static)
elseif(SLATEC_BUILD_SHARED)
    add_library(slatec::slatec ALIAS slatec_shared)
endif()

# Installation
include(GNUInstallDirs)

# Install libraries
set(INSTALL_TARGETS)
if(SLATEC_BUILD_STATIC)
    list(APPEND INSTALL_TARGETS slatec_static)
endif()
if(SLATEC_BUILD_SHARED)
    list(APPEND INSTALL_TARGETS slatec_shared)
endif()

install(TARGETS ${INSTALL_TARGETS}
    EXPORT slatec-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install module files
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/slatec
    FILES_MATCHING PATTERN "*.mod"
)

# CMake package configuration
install(EXPORT slatec-targets
    FILE slatec-targets.cmake
    NAMESPACE slatec::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/slatec
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/slatec-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Generate config file
file(WRITE "${CMAKE_BINARY_DIR}/slatec-config.cmake"
"# SLATEC-Modern CMake Configuration
include(CMakeFindDependencyMacro)
include(\"\${CMAKE_CURRENT_LIST_DIR}/slatec-targets.cmake\")
")

install(FILES
    "${CMAKE_BINARY_DIR}/slatec-config.cmake"
    "${CMAKE_BINARY_DIR}/slatec-config-version.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/slatec
)
# SLATEC-Modern Test Suite
#
# Uses test-drive framework for unit testing.
# https://github.com/fortran-lang/test-drive

include(FetchContent)

# Fetch test-drive if not already available
FetchContent_Declare(
    test-drive
    GIT_REPOSITORY https://github.com/fortran-lang/test-drive
    GIT_TAG v0.4.0
)
FetchContent_MakeAvailable(test-drive)

# test-drive test suite sources (modules only, not standalone programs)
set(TESTDRIVE_SOURCES
    main.f90
    test_utils.f90
    test_bessel_regression.f90
    test_bessel_reference.f90
    test_bessel_portability.f90
    test_specfun_bessel.f90
    test_diff_integ_quadpack.f90
    test_diff_integ_eq.f90
)

# Create test-drive test executable
add_executable(slatec_tests ${TESTDRIVE_SOURCES})

# Link against slatec library and test-drive
if(SLATEC_BUILD_STATIC)
    target_link_libraries(slatec_tests PRIVATE slatec_static test-drive)
elseif(SLATEC_BUILD_SHARED)
    target_link_libraries(slatec_tests PRIVATE slatec_shared test-drive)
endif()

# Add test-drive suite to CTest
add_test(NAME slatec_tests COMMAND slatec_tests)

# Golden tests - standalone programs that don't use test-drive
# These test against A&S reference values and can detect hardware drift
add_executable(test_bessel_golden test_bessel_golden.f90)
if(SLATEC_BUILD_STATIC)
    target_link_libraries(test_bessel_golden PRIVATE slatec_static)
elseif(SLATEC_BUILD_SHARED)
    target_link_libraries(test_bessel_golden PRIVATE slatec_shared)
endif()
add_test(NAME bessel_golden COMMAND test_bessel_golden)

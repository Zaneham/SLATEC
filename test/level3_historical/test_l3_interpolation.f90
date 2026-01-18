!> Level 3: Historical Baseline Tests for Interpolation
!>
!> Purpose: Does our output match what IBM System/360 users saw?
!>
!> These tests compare against "golden" outputs captured from actual
!> IBM System/360 execution via Hercules emulation running MVT/FORTRAN G.
!>
!> If Level 3 fails but Level 2 passes:
!>   - Mathematics is correct
!>   - Platform differs from IBM 360
!>   - Document in DEVIATIONS.md with root cause
!>
!> Platform tested: IBM System/360 (Hercules/TK4-), FORTRAN G (IEYFORT)
!> Date verified: 18 January 2026
!>
!> Note: IBM Hex floating-point differs from IEEE 754:
!>   - Base 16 (not base 2)
!>   - 56-bit mantissa (not 52-bit)
!>   - No gradual underflow, no NaN/Inf
!>   - "Wobbling precision" (0-3 bits)
!>
!> STATUS:
!>   - Polynomial: VERIFIED with actual IBM 360 output
!>   - PCHIP: N/A (algorithm from 1980, post-dates IBM 360 era)
!>   - B-spline: N/A (de Boor algorithms post-date IBM 360 era, plus L2 issues)

module test_interpolation_level3
  use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
  implicit none
  private

  public :: run_all_tests

  ! Tolerances account for IBM Hex vs IEEE 754 differences
  real(dp), parameter :: tol_dp = 1.0e-12_dp
  real(dp), parameter :: tol_hex = 1.0e-10_dp  ! Looser for hex comparisons

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 3: HISTORICAL BASELINE (IBM SYSTEM/360)'
    print '(A)', '================================================================'
    print '(A)', 'Platform: Hercules/TK4-, FORTRAN G (IEYFORT)'
    print '(A)', 'Verified: 18 January 2026'
    print '(A)', ''

    ! Polynomial interpolation - VERIFIED with IBM 360 golden values
    call test_polynomial_ibm360(p, f)
    passed = passed + p
    failed = failed + f

    ! PCHIP - N/A (algorithm post-dates IBM 360)
    call test_pchip_not_applicable(p, f)
    passed = passed + p
    failed = failed + f

    ! B-spline - N/A (algorithm post-dates IBM 360, plus L2 issues)
    call test_bspline_not_applicable(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 3 SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! Polynomial Golden Values from IBM System/360
  !
  ! CAPTURED via: python -m src.fortran360.cli run tests/slatec/interpolation/test_poly.f
  ! Platform: Hercules 3.07, TK4-/MVT, FORTRAN G (IEYFORT)
  ! Date: 18 January 2026
  !
  ! Test case: f(x) = x^3 interpolated through x = 0, 1, 2, 3, 4
  ! Newton divided differences computed explicitly
  ! Evaluation at: 0.5, 1.5, 2.5, 3.5
  !
  ! IBM 360 OUTPUT (verbatim):
  !   NEWTON COEFFICIENTS:
  !     C(1) =  0.0
  !     C(2) =  0.100000000000000D 01
  !     C(3) =  0.300000000000000D 01
  !     C(4) =  0.100000000000000D 01
  !     C(5) =  0.0
  !
  !   TEST 1: p(0.5)
  !     IBM 360: p(0.5) =  0.125000000000000D 00
  !     PASS
  !
  !   TEST 2: p(1.5)
  !     IBM 360: p(1.5) =  0.337500000000000D 01
  !     PASS
  !
  !   TEST 3: p(2.5)
  !     IBM 360: p(2.5) =  0.156250000000000D 02
  !     PASS
  !
  !   TEST 4: p(3.5)
  !     IBM 360: p(3.5) =  0.428750000000000D 02
  !     PASS
  !---------------------------------------------------------------------------
  subroutine test_polynomial_ibm360(passed, failed)
    use interpolation, only: DPLINT, DPOLVL
    integer, intent(out) :: passed, failed
    integer, parameter :: n = 5
    integer, parameter :: ne = 4
    integer :: i, ierr
    real(dp) :: x(n), y(n), c(n)
    real(dp) :: xe(ne), pe(ne), yp(1), work(2*n)

    ! IBM 360 golden values (VERIFIED from Hercules/TK4-)
    ! These are the exact values output by IBM System/360 FORTRAN G
    real(dp), parameter :: pe_ibm360(ne) = [ &
      0.125000000000000e+00_dp, &  ! p(0.5) from IBM 360
      0.337500000000000e+01_dp, &  ! p(1.5) from IBM 360
      0.156250000000000e+02_dp, &  ! p(2.5) from IBM 360
      0.428750000000000e+02_dp  &  ! p(3.5) from IBM 360
    ]

    passed = 0
    failed = 0

    print '(A)', 'Polynomial - IBM 360 Golden Values (VERIFIED)'
    print '(A)', '---------------------------------------------'
    print '(A)', 'Golden values captured from Hercules/TK4-, FORTRAN G'
    print '(A)', 'Test: y = x^3 through x = 0,1,2,3,4'
    print '(A)', ''

    ! Setup: y = x^3 at x = 0, 1, 2, 3, 4
    do i = 1, n
      x(i) = real(i-1, dp)
      y(i) = x(i)**3
    end do

    call DPLINT(n, x, y, c)

    ! Evaluation points
    xe = [0.5_dp, 1.5_dp, 2.5_dp, 3.5_dp]

    do i = 1, ne
      call DPOLVL(0, xe(i), pe(i), yp, n, x, c, work, ierr)

      if (abs(pe(i) - pe_ibm360(i)) < tol_dp) then
        print '(A,F4.1,A,ES22.15)', '  [PASS] p(', xe(i), ') = ', pe(i)
        passed = passed + 1
      else
        print '(A,F4.1,A)', '  [FAIL] p(', xe(i), ') differs from IBM 360'
        print '(A,ES22.15)', '         Modern:  ', pe(i)
        print '(A,ES22.15)', '         IBM 360: ', pe_ibm360(i)
        failed = failed + 1
      end if
    end do

    print '(A)', ''
    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'

  end subroutine test_polynomial_ibm360

  !---------------------------------------------------------------------------
  ! PCHIP - Not Applicable for Level 3
  !
  ! PCHIP (Piecewise Cubic Hermite Interpolating Polynomial) was developed
  ! by Fritsch and Carlson in 1980 (SIAM J. Numer. Anal. 17(2), 238-246).
  !
  ! IBM System/360 FORTRAN G/H compilers are from 1966/1969. PCHIP did not
  ! exist in that era, so there are no historical IBM 360 golden values.
  !
  ! For PCHIP validation, we rely on:
  !   - Level 2 (mathematical verification against Fritsch & Carlson 1980)
  !   - Level 4 (hostile environment/portability testing)
  !---------------------------------------------------------------------------
  subroutine test_pchip_not_applicable(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'PCHIP - Not Applicable (Post-1980 Algorithm)'
    print '(A)', '---------------------------------------------'
    print '(A)', 'PCHIP was developed by Fritsch & Carlson (1980).'
    print '(A)', 'No IBM 360 golden values exist for this algorithm.'
    print '(A)', 'Validation via L2 (mathematical) and L4 (hostile) tests.'
    print '(A)', ''
    print '(A)', '  [N/A] PCHIP tests skipped for Level 3'
    print '(A)', ''

  end subroutine test_pchip_not_applicable

  !---------------------------------------------------------------------------
  ! B-spline - Not Applicable for Level 3
  !
  ! The SLATEC B-spline routines (DBINT4, DBVALU, etc.) are based on
  ! Carl de Boor's work from the 1970s, culminating in "A Practical Guide
  ! to Splines" (1978).
  !
  ! Additionally, Level 2 tests reveal that DBINT4 has mathematical issues:
  !   - Natural spline does not interpolate exactly (max error: 3.56e-03)
  !   - Natural spline boundary condition S''(1) = 60.1 instead of 0
  !
  ! Until DBINT4 is fixed, Level 3 testing is not meaningful.
  !---------------------------------------------------------------------------
  subroutine test_bspline_not_applicable(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'B-spline - Not Applicable (Post-1978 Algorithm + L2 Issues)'
    print '(A)', '------------------------------------------------------------'
    print '(A)', 'B-spline routines based on de Boor (1978).'
    print '(A)', 'No IBM 360 golden values exist for these algorithms.'
    print '(A)', ''
    print '(A)', 'Additionally, L2 tests show DBINT4 has mathematical issues:'
    print '(A)', '  - Natural spline does not interpolate exactly'
    print '(A)', '  - S''''(1) = 60.1 instead of expected 0'
    print '(A)', ''
    print '(A)', '  [N/A] B-spline tests skipped for Level 3'
    print '(A)', ''

  end subroutine test_bspline_not_applicable

end module test_interpolation_level3

!> Main program for Level 3 Interpolation tests
program run_level3_interpolation
  use test_interpolation_level3
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level3_interpolation

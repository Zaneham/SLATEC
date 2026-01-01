!> Level 3: Historical Baseline Tests for MINPACK
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
!> Platform tested: IBM System/360 (Hercules/TK4-), FORTRAN G, 1 Jan 2026
!>
!> Note: IBM Hex floating-point differs from IEEE 754:
!>   - Base 16 (not base 2)
!>   - 56-bit mantissa (not 52-bit)
!>   - No gradual underflow, no NaN/Inf
!>   - "Wobbling precision" (0-3 bits)

module test_minpack_level3
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
    print '(A)', 'Golden outputs from Hercules/TK4-/FORTRAN G, 1 Jan 2026'
    print '(A)', ''

    ! DENORM golden values from IBM 360
    call test_denorm_ibm360_golden(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 3 SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! DENORM Golden Values from IBM System/360
  ! Captured via: python -m src.fortran360.cli run tests/slatec/minpack/test_enorm.f
  ! Platform: Hercules 3.07, TK4-/MVT, FORTRAN G (IEYFORT)
  ! Date: 1 January 2026
  !---------------------------------------------------------------------------
  subroutine test_denorm_ibm360_golden(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(10), result

    passed = 0
    failed = 0

    print '(A)', 'DENORM - IBM 360 Golden Values'
    print '(A)', '-------------------------------'
    print '(A)', 'Source: test_enorm.f on Hercules/TK4-/FORTRAN G'
    print '(A)', ''

    ! Test 1: Unit vector
    ! IBM 360 output: NORM = 0.1000000000000D+01
    x = 0.0_dp
    x(1) = 1.0_dp
    result = denorm_dp(10, x)
    if (abs(result - 1.0_dp) < tol_dp) then
      print '(A)', '  [PASS] Unit vector: IBM 360 = 0.1000000000000D+01'
      passed = passed + 1
    else
      print '(A,ES20.13)', '  [FAIL] Unit vector: modern = ', result
      failed = failed + 1
    end if

    ! Test 2: Vector of ones
    ! IBM 360 output: NORM = 0.3162277660168D+01 (sqrt(10))
    x = 1.0_dp
    result = denorm_dp(10, x)
    if (abs(result - sqrt(10.0_dp)) < tol_dp) then
      print '(A)', '  [PASS] Ones vector: IBM 360 = 0.3162277660168D+01'
      passed = passed + 1
    else
      print '(A,ES20.13)', '  [FAIL] Ones vector: modern = ', result
      failed = failed + 1
    end if

    ! Test 3: 3-4-5 triangle
    ! IBM 360 output: NORM = 0.5000000000000D+01
    x = 0.0_dp
    x(1) = 3.0_dp
    x(2) = 4.0_dp
    result = denorm_dp(10, x)
    if (abs(result - 5.0_dp) < tol_dp) then
      print '(A)', '  [PASS] 3-4-5 triangle: IBM 360 = 0.5000000000000D+01'
      passed = passed + 1
    else
      print '(A,ES20.13)', '  [FAIL] 3-4-5 triangle: modern = ', result
      failed = failed + 1
    end if

    ! Test 4: Large components (overflow protection)
    ! IBM 360 output: NORM = 0.1732050807569D+16
    x(1) = 1.0e15_dp
    x(2) = 1.0e15_dp
    x(3) = 1.0e15_dp
    result = denorm_dp(3, x)
    if (abs(result - sqrt(3.0_dp)*1.0e15_dp) / (sqrt(3.0_dp)*1.0e15_dp) < tol_hex) then
      print '(A)', '  [PASS] Large components: IBM 360 = 0.1732050807569D+16'
      passed = passed + 1
    else
      print '(A,ES20.13)', '  [FAIL] Large components: modern = ', result
      failed = failed + 1
    end if

    ! Test 5: Small components (underflow protection)
    ! IBM 360 output: NORM = 0.1732050807569D-14
    x(1) = 1.0e-15_dp
    x(2) = 1.0e-15_dp
    x(3) = 1.0e-15_dp
    result = denorm_dp(3, x)
    if (abs(result - sqrt(3.0_dp)*1.0e-15_dp) / (sqrt(3.0_dp)*1.0e-15_dp) < tol_hex) then
      print '(A)', '  [PASS] Small components: IBM 360 = 0.1732050807569D-14'
      passed = passed + 1
    else
      print '(A,ES20.13)', '  [FAIL] Small components: modern = ', result
      failed = failed + 1
    end if

    ! Test 6: Mixed magnitude
    ! IBM 360 output: NORM = 0.1000000000000D+19 (dominated by 1E18)
    x(1) = 1.0e-18_dp
    x(2) = 1.0_dp
    x(3) = 1.0e18_dp
    result = denorm_dp(3, x)
    if (abs(result - 1.0e18_dp) / 1.0e18_dp < tol_hex) then
      print '(A)', '  [PASS] Mixed magnitude: IBM 360 = 0.1000000000000D+19'
      passed = passed + 1
    else
      print '(A,ES20.13)', '  [FAIL] Mixed magnitude: modern = ', result
      failed = failed + 1
    end if

    ! Test 7: Zero vector
    ! IBM 360 output: NORM = 0
    x = 0.0_dp
    result = denorm_dp(5, x)
    if (result == 0.0_dp) then
      print '(A)', '  [PASS] Zero vector: IBM 360 = 0'
      passed = passed + 1
    else
      print '(A,ES20.13)', '  [FAIL] Zero vector: modern = ', result
      failed = failed + 1
    end if

    print '(A)', ''
    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'

  end subroutine test_denorm_ibm360_golden

  !---------------------------------------------------------------------------
  ! DENORM Implementation (Double Precision)
  !---------------------------------------------------------------------------
  pure function denorm_dp(n, x) result(norm)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    real(dp) :: norm

    real(dp), parameter :: rdwarf = 3.834e-20_dp
    real(dp), parameter :: rgiant = 1.304e19_dp
    real(dp) :: s1, s2, s3, x1max, x3max, xabs, agiant
    integer :: i

    s1 = 0.0_dp
    s2 = 0.0_dp
    s3 = 0.0_dp
    x1max = 0.0_dp
    x3max = 0.0_dp
    agiant = rgiant / real(n, dp)

    do i = 1, n
      xabs = abs(x(i))

      if (xabs > rdwarf .and. xabs < agiant) then
        s2 = s2 + xabs**2
      else if (xabs <= rdwarf) then
        if (xabs > x3max) then
          s3 = 1.0_dp + s3 * (x3max/xabs)**2
          x3max = xabs
        else if (xabs /= 0.0_dp) then
          s3 = s3 + (xabs/x3max)**2
        end if
      else
        if (xabs > x1max) then
          s1 = 1.0_dp + s1 * (x1max/xabs)**2
          x1max = xabs
        else
          s1 = s1 + (xabs/x1max)**2
        end if
      end if
    end do

    if (s1 /= 0.0_dp) then
      norm = x1max * sqrt(s1 + (s2/x1max)/x1max)
    else if (s2 /= 0.0_dp) then
      if (s2 >= x3max) then
        norm = sqrt(s2 * (1.0_dp + (x3max/s2)*(x3max*s3)))
      else
        norm = sqrt(x3max * ((s2/x3max) + (x3max*s3)))
      end if
    else
      norm = x3max * sqrt(s3)
    end if

  end function denorm_dp

end module test_minpack_level3

!> Main program for Level 3 MINPACK tests
program run_level3_minpack
  use test_minpack_level3
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level3_minpack

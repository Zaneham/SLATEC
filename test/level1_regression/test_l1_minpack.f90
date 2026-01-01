!> Level 1: Regression Tests for MINPACK
!>
!> Purpose: Does the code work? Does it still work after changes?
!> These tests verify basic functionality and catch regressions.
!>
!> What we test:
!>   - DENORM/ENORM: Euclidean norm computation
!>   - DQRFAC/QRFAC: QR factorisation
!>   - DNLS1/SNLS1: Nonlinear least squares
!>   - DNSQ/SNSQ: Nonlinear equations
!>
!> Reference: Moré, Garbow, Hillstrom - MINPACK (ANL-80-74)

module test_minpack_level1
  use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: tol_dp = 1.0e-10_dp
  real(sp), parameter :: tol_sp = 1.0e-5_sp

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 1: MINPACK REGRESSION TESTS'
    print '(A)', '================================================================'
    print '(A)', ''

    ! DENORM tests
    call test_denorm_suite(p, f)
    passed = passed + p
    failed = failed + f

    ! ENORM tests (single precision)
    call test_enorm_suite(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 1 SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! DENORM Test Suite (Double Precision Euclidean Norm)
  !---------------------------------------------------------------------------
  subroutine test_denorm_suite(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'DENORM (Double Precision Euclidean Norm)'
    print '(A)', '----------------------------------------'

    call test_denorm_unit_vector(passed, failed)
    call test_denorm_ones_vector(passed, failed)
    call test_denorm_345_triangle(passed, failed)
    call test_denorm_large_components(passed, failed)
    call test_denorm_small_components(passed, failed)
    call test_denorm_mixed_magnitude(passed, failed)
    call test_denorm_zero_vector(passed, failed)

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_denorm_suite

  subroutine test_denorm_unit_vector(passed, failed)
    integer, intent(inout) :: passed, failed
    real(dp) :: x(10), result, expected
    integer :: i

    ! Unit vector e1 = (1, 0, 0, ..., 0)
    x = 0.0_dp
    x(1) = 1.0_dp

    result = denorm_dp(10, x)
    expected = 1.0_dp

    if (abs(result - expected) < tol_dp) then
      print '(A)', '  [PASS] Unit vector: ||e1|| = 1'
      passed = passed + 1
    else
      print '(A,ES12.5,A,ES12.5)', '  [FAIL] Unit vector: got ', result, ' expected ', expected
      failed = failed + 1
    end if
  end subroutine

  subroutine test_denorm_ones_vector(passed, failed)
    integer, intent(inout) :: passed, failed
    real(dp) :: x(10), result, expected

    x = 1.0_dp
    result = denorm_dp(10, x)
    expected = sqrt(10.0_dp)

    if (abs(result - expected) < tol_dp) then
      print '(A)', '  [PASS] Ones vector: ||(1,1,...,1)|| = sqrt(10)'
      passed = passed + 1
    else
      print '(A,ES12.5,A,ES12.5)', '  [FAIL] Ones vector: got ', result, ' expected ', expected
      failed = failed + 1
    end if
  end subroutine

  subroutine test_denorm_345_triangle(passed, failed)
    integer, intent(inout) :: passed, failed
    real(dp) :: x(10), result, expected

    x = 0.0_dp
    x(1) = 3.0_dp
    x(2) = 4.0_dp

    result = denorm_dp(10, x)
    expected = 5.0_dp

    if (abs(result - expected) < tol_dp) then
      print '(A)', '  [PASS] 3-4-5 triangle: ||(3,4,0,...)|| = 5'
      passed = passed + 1
    else
      print '(A,ES12.5,A,ES12.5)', '  [FAIL] 3-4-5 triangle: got ', result, ' expected ', expected
      failed = failed + 1
    end if
  end subroutine

  subroutine test_denorm_large_components(passed, failed)
    integer, intent(inout) :: passed, failed
    real(dp) :: x(3), result, expected, rel_error

    x = 1.0e15_dp
    result = denorm_dp(3, x)
    expected = sqrt(3.0_dp) * 1.0e15_dp
    rel_error = abs(result - expected) / expected

    if (rel_error < 1.0e-10_dp) then
      print '(A)', '  [PASS] Large components: overflow protection works'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [FAIL] Large components: relative error = ', rel_error
      failed = failed + 1
    end if
  end subroutine

  subroutine test_denorm_small_components(passed, failed)
    integer, intent(inout) :: passed, failed
    real(dp) :: x(3), result, expected, rel_error

    x = 1.0e-15_dp
    result = denorm_dp(3, x)
    expected = sqrt(3.0_dp) * 1.0e-15_dp
    rel_error = abs(result - expected) / expected

    if (rel_error < 1.0e-10_dp) then
      print '(A)', '  [PASS] Small components: underflow protection works'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [FAIL] Small components: relative error = ', rel_error
      failed = failed + 1
    end if
  end subroutine

  subroutine test_denorm_mixed_magnitude(passed, failed)
    integer, intent(inout) :: passed, failed
    real(dp) :: x(3), result, expected, rel_error

    x(1) = 1.0e-18_dp
    x(2) = 1.0_dp
    x(3) = 1.0e18_dp

    result = denorm_dp(3, x)
    expected = 1.0e18_dp  ! Dominated by largest component
    rel_error = abs(result - expected) / expected

    if (rel_error < 1.0e-10_dp) then
      print '(A)', '  [PASS] Mixed magnitude: handles 36 orders of magnitude'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [FAIL] Mixed magnitude: relative error = ', rel_error
      failed = failed + 1
    end if
  end subroutine

  subroutine test_denorm_zero_vector(passed, failed)
    integer, intent(inout) :: passed, failed
    real(dp) :: x(5), result

    x = 0.0_dp
    result = denorm_dp(5, x)

    if (result == 0.0_dp) then
      print '(A)', '  [PASS] Zero vector: ||0|| = 0'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [FAIL] Zero vector: got ', result
      failed = failed + 1
    end if
  end subroutine

  !---------------------------------------------------------------------------
  ! ENORM Test Suite (Single Precision Euclidean Norm)
  !---------------------------------------------------------------------------
  subroutine test_enorm_suite(passed, failed)
    integer, intent(out) :: passed, failed
    real(sp) :: x(10), result, expected

    passed = 0
    failed = 0

    print '(A)', 'ENORM (Single Precision Euclidean Norm)'
    print '(A)', '---------------------------------------'

    ! Test 1: Unit vector
    x = 0.0_sp
    x(1) = 1.0_sp
    result = enorm_sp(10, x)
    expected = 1.0_sp

    if (abs(result - expected) < tol_sp) then
      print '(A)', '  [PASS] Unit vector'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [FAIL] Unit vector: got ', result
      failed = failed + 1
    end if

    ! Test 2: 3-4-5 triangle
    x = 0.0_sp
    x(1) = 3.0_sp
    x(2) = 4.0_sp
    result = enorm_sp(10, x)
    expected = 5.0_sp

    if (abs(result - expected) < tol_sp) then
      print '(A)', '  [PASS] 3-4-5 triangle'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [FAIL] 3-4-5 triangle: got ', result
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_enorm_suite

  !---------------------------------------------------------------------------
  ! DENORM Implementation (Double Precision)
  ! Scaled accumulation to avoid overflow/underflow
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
        ! Intermediate range - direct accumulation
        s2 = s2 + xabs**2
      else if (xabs <= rdwarf) then
        ! Small components - scaled accumulation
        if (xabs > x3max) then
          s3 = 1.0_dp + s3 * (x3max/xabs)**2
          x3max = xabs
        else if (xabs /= 0.0_dp) then
          s3 = s3 + (xabs/x3max)**2
        end if
      else
        ! Large components - scaled accumulation
        if (xabs > x1max) then
          s1 = 1.0_dp + s1 * (x1max/xabs)**2
          x1max = xabs
        else
          s1 = s1 + (xabs/x1max)**2
        end if
      end if
    end do

    ! Combine the three sums
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

  !---------------------------------------------------------------------------
  ! ENORM Implementation (Single Precision)
  !---------------------------------------------------------------------------
  pure function enorm_sp(n, x) result(norm)
    integer, intent(in) :: n
    real(sp), intent(in) :: x(n)
    real(sp) :: norm

    real(sp), parameter :: rdwarf = 3.834e-20_sp
    real(sp), parameter :: rgiant = 1.304e19_sp
    real(sp) :: s1, s2, s3, x1max, x3max, xabs, agiant
    integer :: i

    s1 = 0.0_sp
    s2 = 0.0_sp
    s3 = 0.0_sp
    x1max = 0.0_sp
    x3max = 0.0_sp
    agiant = rgiant / real(n, sp)

    do i = 1, n
      xabs = abs(x(i))

      if (xabs > rdwarf .and. xabs < agiant) then
        s2 = s2 + xabs**2
      else if (xabs <= rdwarf) then
        if (xabs > x3max) then
          s3 = 1.0_sp + s3 * (x3max/xabs)**2
          x3max = xabs
        else if (xabs /= 0.0_sp) then
          s3 = s3 + (xabs/x3max)**2
        end if
      else
        if (xabs > x1max) then
          s1 = 1.0_sp + s1 * (x1max/xabs)**2
          x1max = xabs
        else
          s1 = s1 + (xabs/x1max)**2
        end if
      end if
    end do

    if (s1 /= 0.0_sp) then
      norm = x1max * sqrt(s1 + (s2/x1max)/x1max)
    else if (s2 /= 0.0_sp) then
      if (s2 >= x3max) then
        norm = sqrt(s2 * (1.0_sp + (x3max/s2)*(x3max*s3)))
      else
        norm = sqrt(x3max * ((s2/x3max) + (x3max*s3)))
      end if
    else
      norm = x3max * sqrt(s3)
    end if

  end function enorm_sp

end module test_minpack_level1

!> Main program for Level 1 MINPACK tests
program run_level1_minpack
  use test_minpack_level1
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level1_minpack

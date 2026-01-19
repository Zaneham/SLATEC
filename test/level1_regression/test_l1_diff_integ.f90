!> Level 1: Regression Tests for Differentiation/Integration (diff_integ)
!>
!> Purpose: Did the last code change break anything?
!>
!> These tests lock in current behaviour. If a test fails after code changes,
!> either the change introduced a bug, or the test needs updating (with justification).
!>
!> Routines tested:
!>   - DGAUS8: Adaptive 8-point Legendre-Gauss quadrature
!>   - DQAGS: QUADPACK general-purpose adaptive integrator
!>   - DQAGI: QUADPACK infinite interval integrator
!>   - DQNG: QUADPACK non-adaptive Gauss-Kronrod

module test_diff_integ_level1
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use diff_integ, only: DGAUS8, DQAGS, DQAGI, DQNG
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: tol = 1.0e-10_dp
  real(dp), parameter :: pi = 3.14159265358979323846_dp

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 1: REGRESSION TESTS (diff_integ)'
    print '(A)', '================================================================'
    print '(A)', ''

    call test_dgaus8_basic(p, f)
    passed = passed + p
    failed = failed + f

    call test_dqags_basic(p, f)
    passed = passed + p
    failed = failed + f

    call test_dqagi_basic(p, f)
    passed = passed + p
    failed = failed + f

    call test_dqng_basic(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 1 SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! DGAUS8 Basic Tests
  !---------------------------------------------------------------------------
  subroutine test_dgaus8_basic(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: result, err, exact
    integer :: ierr

    passed = 0
    failed = 0

    print '(A)', 'DGAUS8 - Adaptive Gauss-Legendre Quadrature'
    print '(A)', '-------------------------------------------'

    ! Test 1: Integral of x^2 from 0 to 1 = 1/3
    err = 1.0e-12_dp
    call DGAUS8(f_x_squared, 0.0_dp, 1.0_dp, err, result, ierr)
    exact = 1.0_dp / 3.0_dp
    if (abs(result - exact) < tol .and. ierr == 1) then
      print '(A,ES15.8)', '  [PASS] integral(x^2, 0, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8)', '  [FAIL] integral(x^2, 0, 1) = ', result, ' expected ', exact
      failed = failed + 1
    end if

    ! Test 2: Integral of sin(x) from 0 to pi = 2
    err = 1.0e-12_dp
    call DGAUS8(f_sin, 0.0_dp, pi, err, result, ierr)
    exact = 2.0_dp
    if (abs(result - exact) < tol .and. ierr == 1) then
      print '(A,ES15.8)', '  [PASS] integral(sin(x), 0, pi) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8)', '  [FAIL] integral(sin(x), 0, pi) = ', result, ' expected ', exact
      failed = failed + 1
    end if

    ! Test 3: Integral of exp(-x) from 0 to 1 = 1 - 1/e
    err = 1.0e-12_dp
    call DGAUS8(f_exp_neg, 0.0_dp, 1.0_dp, err, result, ierr)
    exact = 1.0_dp - exp(-1.0_dp)
    if (abs(result - exact) < tol .and. ierr == 1) then
      print '(A,ES15.8)', '  [PASS] integral(exp(-x), 0, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8)', '  [FAIL] integral(exp(-x), 0, 1) = ', result, ' expected ', exact
      failed = failed + 1
    end if

    print '(A)', ''

  end subroutine test_dgaus8_basic

  !---------------------------------------------------------------------------
  ! DQAGS Basic Tests (QUADPACK general-purpose)
  !---------------------------------------------------------------------------
  subroutine test_dqags_basic(passed, failed)
    integer, intent(out) :: passed, failed
    integer, parameter :: limit = 500
    integer, parameter :: lenw = 4 * limit
    real(dp) :: result, abserr, exact
    integer :: neval, ier, last
    integer :: iwork(limit)
    real(dp) :: work(lenw)

    passed = 0
    failed = 0

    print '(A)', 'DQAGS - QUADPACK Adaptive Integrator'
    print '(A)', '-------------------------------------'

    ! Test 1: Integral of 1/sqrt(x) from 0 to 1 = 2 (integrable singularity)
    call DQAGS(f_inv_sqrt, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-10_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = 2.0_dp
    if (abs(result - exact) < 1.0e-8_dp .and. ier == 0) then
      print '(A,ES15.8)', '  [PASS] integral(1/sqrt(x), 0, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8,A,I2)', '  [FAIL] integral(1/sqrt(x), 0, 1) = ', result, &
            ' expected ', exact, ' ier=', ier
      failed = failed + 1
    end if

    ! Test 2: Integral of log(x) from 0 to 1 = -1 (integrable singularity)
    call DQAGS(f_log, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-10_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = -1.0_dp
    if (abs(result - exact) < 1.0e-8_dp .and. ier == 0) then
      print '(A,ES15.8)', '  [PASS] integral(log(x), 0, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8,A,I2)', '  [FAIL] integral(log(x), 0, 1) = ', result, &
            ' expected ', exact, ' ier=', ier
      failed = failed + 1
    end if

    ! Test 3: Integral of x*exp(-x^2) from 0 to 1 = (1 - 1/e)/2
    call DQAGS(f_x_exp_neg_x2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-10_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = (1.0_dp - exp(-1.0_dp)) / 2.0_dp
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES15.8)', '  [PASS] integral(x*exp(-x^2), 0, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8,A,I2)', '  [FAIL] integral(x*exp(-x^2), 0, 1) = ', result, &
            ' expected ', exact, ' ier=', ier
      failed = failed + 1
    end if

    print '(A)', ''

  end subroutine test_dqags_basic

  !---------------------------------------------------------------------------
  ! DQAGI Basic Tests (QUADPACK infinite intervals)
  !---------------------------------------------------------------------------
  subroutine test_dqagi_basic(passed, failed)
    integer, intent(out) :: passed, failed
    integer, parameter :: limit = 500
    integer, parameter :: lenw = 4 * limit
    real(dp) :: result, abserr, exact
    integer :: neval, ier, last
    integer :: iwork(limit)
    real(dp) :: work(lenw)

    passed = 0
    failed = 0

    print '(A)', 'DQAGI - QUADPACK Infinite Interval Integrator'
    print '(A)', '----------------------------------------------'

    ! Test 1: Integral of exp(-x) from 0 to infinity = 1
    ! inf = 1 means (bound, +infinity)
    call DQAGI(f_exp_neg, 0.0_dp, 1, 0.0_dp, 1.0e-10_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = 1.0_dp
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES15.8)', '  [PASS] integral(exp(-x), 0, inf) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8,A,I2)', '  [FAIL] integral(exp(-x), 0, inf) = ', result, &
            ' expected ', exact, ' ier=', ier
      failed = failed + 1
    end if

    ! Test 2: Integral of exp(-x^2) from -infinity to infinity = sqrt(pi)
    ! inf = 2 means (-infinity, +infinity)
    call DQAGI(f_exp_neg_x2, 0.0_dp, 2, 0.0_dp, 1.0e-10_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = sqrt(pi)
    if (abs(result - exact) < 1.0e-8_dp .and. ier == 0) then
      print '(A,ES15.8)', '  [PASS] integral(exp(-x^2), -inf, inf) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8,A,I2)', '  [FAIL] integral(exp(-x^2), -inf, inf) = ', result, &
            ' expected ', exact, ' ier=', ier
      failed = failed + 1
    end if

    ! Test 3: Integral of 1/(1+x^2) from 0 to infinity = pi/2
    call DQAGI(f_lorentzian, 0.0_dp, 1, 0.0_dp, 1.0e-10_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = pi / 2.0_dp
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES15.8)', '  [PASS] integral(1/(1+x^2), 0, inf) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8,A,I2)', '  [FAIL] integral(1/(1+x^2), 0, inf) = ', result, &
            ' expected ', exact, ' ier=', ier
      failed = failed + 1
    end if

    print '(A)', ''

  end subroutine test_dqagi_basic

  !---------------------------------------------------------------------------
  ! DQNG Basic Tests (QUADPACK non-adaptive)
  !---------------------------------------------------------------------------
  subroutine test_dqng_basic(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: result, abserr, exact
    integer :: neval, ier

    passed = 0
    failed = 0

    print '(A)', 'DQNG - QUADPACK Non-adaptive Gauss-Kronrod'
    print '(A)', '-------------------------------------------'

    ! Test 1: Integral of x^3 from 0 to 1 = 1/4
    call DQNG(f_x_cubed, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-10_dp, &
              result, abserr, neval, ier)
    exact = 0.25_dp
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES15.8)', '  [PASS] integral(x^3, 0, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8,A,I2)', '  [FAIL] integral(x^3, 0, 1) = ', result, &
            ' expected ', exact, ' ier=', ier
      failed = failed + 1
    end if

    ! Test 2: Integral of cos(x) from 0 to pi/2 = 1
    call DQNG(f_cos, 0.0_dp, pi/2.0_dp, 0.0_dp, 1.0e-10_dp, &
              result, abserr, neval, ier)
    exact = 1.0_dp
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES15.8)', '  [PASS] integral(cos(x), 0, pi/2) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8,A,I2)', '  [FAIL] integral(cos(x), 0, pi/2) = ', result, &
            ' expected ', exact, ' ier=', ier
      failed = failed + 1
    end if

    ! Test 3: Integral of 1/(1+x^2) from 0 to 1 = pi/4
    call DQNG(f_lorentzian, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-10_dp, &
              result, abserr, neval, ier)
    exact = pi / 4.0_dp
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES15.8)', '  [PASS] integral(1/(1+x^2), 0, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8,A,I2)', '  [FAIL] integral(1/(1+x^2), 0, 1) = ', result, &
            ' expected ', exact, ' ier=', ier
      failed = failed + 1
    end if

    print '(A)', ''

  end subroutine test_dqng_basic

  !---------------------------------------------------------------------------
  ! Test functions (integrands)
  !---------------------------------------------------------------------------

  pure real(dp) function f_x_squared(x)
    real(dp), intent(in) :: x
    f_x_squared = x * x
  end function f_x_squared

  pure real(dp) function f_x_cubed(x)
    real(dp), intent(in) :: x
    f_x_cubed = x * x * x
  end function f_x_cubed

  pure real(dp) function f_sin(x)
    real(dp), intent(in) :: x
    f_sin = sin(x)
  end function f_sin

  pure real(dp) function f_cos(x)
    real(dp), intent(in) :: x
    f_cos = cos(x)
  end function f_cos

  pure real(dp) function f_exp_neg(x)
    real(dp), intent(in) :: x
    f_exp_neg = exp(-x)
  end function f_exp_neg

  pure real(dp) function f_exp_neg_x2(x)
    real(dp), intent(in) :: x
    f_exp_neg_x2 = exp(-x*x)
  end function f_exp_neg_x2

  pure real(dp) function f_x_exp_neg_x2(x)
    real(dp), intent(in) :: x
    f_x_exp_neg_x2 = x * exp(-x*x)
  end function f_x_exp_neg_x2

  pure real(dp) function f_inv_sqrt(x)
    real(dp), intent(in) :: x
    if (x > 0.0_dp) then
      f_inv_sqrt = 1.0_dp / sqrt(x)
    else
      f_inv_sqrt = 0.0_dp  ! Avoid division by zero at endpoint
    end if
  end function f_inv_sqrt

  pure real(dp) function f_log(x)
    real(dp), intent(in) :: x
    if (x > 0.0_dp) then
      f_log = log(x)
    else
      f_log = 0.0_dp  ! Avoid log(0) at endpoint
    end if
  end function f_log

  pure real(dp) function f_lorentzian(x)
    real(dp), intent(in) :: x
    f_lorentzian = 1.0_dp / (1.0_dp + x*x)
  end function f_lorentzian

end module test_diff_integ_level1

!> Main program for Level 1 diff_integ tests
program run_level1_diff_integ
  use test_diff_integ_level1
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level1_diff_integ

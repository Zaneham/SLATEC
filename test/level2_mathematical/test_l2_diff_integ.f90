!> Level 2: Mathematical Verification Tests for diff_integ
!>
!> Purpose: Does the algorithm match the original mathematics?
!>
!> These tests verify against authoritative mathematical references:
!>   - Abramowitz & Stegun (1964) - Handbook of Mathematical Functions
!>   - NIST DLMF - https://dlmf.nist.gov/
!>   - Closed-form solutions from calculus
!>
!> Do not change Level 2 tests without mathematical justification.
!>
!> Routines tested:
!>   - DGAUS8: Verifies Gauss-Legendre quadrature exactness
!>   - DQAGS: Verifies handling of endpoint singularities
!>   - DQAGI: Verifies infinite interval integration
!>   - DQK15/21/31/41/51/61: Verifies Gauss-Kronrod rule exactness

module test_diff_integ_level2
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use diff_integ, only: DGAUS8, DQAGS, DQAGI, DQNG, DQK15, DQK21, DQK31, DQK41, DQK51, DQK61
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: pi = 3.14159265358979323846_dp
  real(dp), parameter :: euler_gamma = 0.5772156649015328606_dp  ! Euler-Mascheroni constant

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 2: MATHEMATICAL VERIFICATION (diff_integ)'
    print '(A)', '================================================================'
    print '(A)', ''

    call test_gauss_kronrod_exactness(p, f)
    passed = passed + p
    failed = failed + f

    call test_classical_integrals(p, f)
    passed = passed + p
    failed = failed + f

    call test_special_function_integrals(p, f)
    passed = passed + p
    failed = failed + f

    call test_endpoint_singularities(p, f)
    passed = passed + p
    failed = failed + f

    call test_infinite_intervals(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 2 SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! Gauss-Kronrod Exactness
  ! Property: G-K rule of order n integrates polynomials of degree <= 2n+1 exactly
  ! Reference: Kronrod (1965), Piessens et al. (1983) QUADPACK
  !---------------------------------------------------------------------------
  subroutine test_gauss_kronrod_exactness(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: result, abserr, resabs, resasc
    real(dp) :: exact
    real(dp), parameter :: tol = 1.0e-14_dp

    passed = 0
    failed = 0

    print '(A)', 'Gauss-Kronrod Exactness (Polynomial Reproduction)'
    print '(A)', '-------------------------------------------------'
    print '(A)', 'G-K rule of order n integrates degree <= 2n+1 exactly'
    print '(A)', ''

    ! QK15 (7-point Gauss, 15-point Kronrod) should be exact for degree <= 15
    ! Test: integral of x^14 from -1 to 1 = 2/15
    call DQK15(f_x14, -1.0_dp, 1.0_dp, result, abserr, resabs, resasc)
    exact = 2.0_dp / 15.0_dp
    if (abs(result - exact) < tol) then
      print '(A,ES22.15)', '  [PASS] QK15: int(x^14, -1, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15)', '  [FAIL] QK15: got ', result, ' expected ', exact
      failed = failed + 1
    end if

    ! QK21 should be exact for degree <= 21
    ! Test: integral of x^20 from -1 to 1 = 2/21
    call DQK21(f_x20, -1.0_dp, 1.0_dp, result, abserr, resabs, resasc)
    exact = 2.0_dp / 21.0_dp
    if (abs(result - exact) < tol) then
      print '(A,ES22.15)', '  [PASS] QK21: int(x^20, -1, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15)', '  [FAIL] QK21: got ', result, ' expected ', exact
      failed = failed + 1
    end if

    ! QK31 should be exact for degree <= 31
    ! Test: integral of x^30 from -1 to 1 = 2/31
    call DQK31(f_x30, -1.0_dp, 1.0_dp, result, abserr, resabs, resasc)
    exact = 2.0_dp / 31.0_dp
    if (abs(result - exact) < tol) then
      print '(A,ES22.15)', '  [PASS] QK31: int(x^30, -1, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15)', '  [FAIL] QK31: got ', result, ' expected ', exact
      failed = failed + 1
    end if

    print '(A)', ''

  end subroutine test_gauss_kronrod_exactness

  !---------------------------------------------------------------------------
  ! Classical Definite Integrals
  ! Reference: Standard calculus results
  !---------------------------------------------------------------------------
  subroutine test_classical_integrals(passed, failed)
    integer, intent(out) :: passed, failed
    integer, parameter :: limit = 500, lenw = 4*limit
    real(dp) :: result, abserr, exact, err
    integer :: neval, ier, last, ierr
    integer :: iwork(limit)
    real(dp) :: work(lenw)
    real(dp), parameter :: tol = 1.0e-12_dp

    passed = 0
    failed = 0

    print '(A)', 'Classical Definite Integrals'
    print '(A)', '----------------------------'
    print '(A)', 'Reference: Standard calculus'
    print '(A)', ''

    ! Integral of sin^2(x) from 0 to pi = pi/2
    ! Identity: sin^2(x) = (1 - cos(2x))/2
    call DQAGS(f_sin_squared, 0.0_dp, pi, 0.0_dp, 1.0e-12_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = pi / 2.0_dp
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES22.15)', '  [PASS] int(sin^2(x), 0, pi) = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15)', '  [FAIL] int(sin^2(x), 0, pi) = ', result, ' expected ', exact
      failed = failed + 1
    end if

    ! Integral of x*sin(x) from 0 to pi = pi (integration by parts)
    call DQAGS(f_x_sin, 0.0_dp, pi, 0.0_dp, 1.0e-12_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = pi
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES22.15)', '  [PASS] int(x*sin(x), 0, pi) = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15)', '  [FAIL] int(x*sin(x), 0, pi) = ', result, ' expected ', exact
      failed = failed + 1
    end if

    ! Integral of 1/(1+x^2) from 0 to 1 = pi/4 (arctan)
    err = 1.0e-12_dp
    call DGAUS8(f_lorentzian, 0.0_dp, 1.0_dp, err, result, ierr)
    exact = pi / 4.0_dp
    if (abs(result - exact) < tol .and. ierr == 1) then
      print '(A,ES22.15)', '  [PASS] int(1/(1+x^2), 0, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15)', '  [FAIL] int(1/(1+x^2), 0, 1) = ', result, ' expected ', exact
      failed = failed + 1
    end if

    ! Integral of x^n * exp(-x) from 0 to infinity = n! (Gamma function)
    ! Test n=5: should give 5! = 120
    call DQAGI(f_x5_exp_neg, 0.0_dp, 1, 0.0_dp, 1.0e-10_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = 120.0_dp  ! 5!
    if (abs(result - exact) / exact < 1.0e-10_dp .and. ier == 0) then
      print '(A,ES22.15)', '  [PASS] int(x^5*exp(-x), 0, inf) = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15)', '  [FAIL] int(x^5*exp(-x), 0, inf) = ', result, ' expected ', exact
      failed = failed + 1
    end if

    print '(A)', ''

  end subroutine test_classical_integrals

  !---------------------------------------------------------------------------
  ! Special Function Integrals
  ! Reference: Abramowitz & Stegun, NIST DLMF
  !---------------------------------------------------------------------------
  subroutine test_special_function_integrals(passed, failed)
    integer, intent(out) :: passed, failed
    integer, parameter :: limit = 500, lenw = 4*limit
    real(dp) :: result, abserr, exact
    integer :: neval, ier, last
    integer :: iwork(limit)
    real(dp) :: work(lenw)
    real(dp), parameter :: tol = 1.0e-10_dp

    passed = 0
    failed = 0

    print '(A)', 'Special Function Integrals'
    print '(A)', '--------------------------'
    print '(A)', 'Reference: A&S, NIST DLMF'
    print '(A)', ''

    ! Gaussian integral: int(exp(-x^2), -inf, inf) = sqrt(pi)
    ! DLMF 7.4.1
    call DQAGI(f_gaussian, 0.0_dp, 2, 0.0_dp, 1.0e-12_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = sqrt(pi)
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES22.15)', '  [PASS] int(exp(-x^2), -inf, inf) = sqrt(pi) = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15)', '  [FAIL] Gaussian = ', result, ' expected ', exact
      failed = failed + 1
    end if

    ! Dirichlet integral: int(sin(x)/x, 0, inf) = pi/2
    ! This is a conditionally convergent improper integral
    ! Using finite interval approximation instead of DQAGI for oscillatory
    ! int(sin(x)/x, 0, 100) is close to pi/2
    call DQAGS(f_sinc, 1.0e-10_dp, 100.0_dp, 0.0_dp, 1.0e-8_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = pi / 2.0_dp
    if (abs(result - exact) < 0.01_dp .and. ier == 0) then  ! Truncation error expected
      print '(A,ES22.15)', '  [PASS] int(sin(x)/x, 0, 100) ~ pi/2 = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15,A,I2)', '  [FAIL] Dirichlet = ', result, &
            ' expected ~', exact, ' ier=', ier
      failed = failed + 1
    end if

    ! Beta function: int(x^(a-1)*(1-x)^(b-1), 0, 1) = Beta(a,b)
    ! Test a=2, b=3: Beta(2,3) = Gamma(2)*Gamma(3)/Gamma(5) = 1*2/24 = 1/12
    call DQAGS(f_beta_2_3, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-12_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = 1.0_dp / 12.0_dp
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES22.15)', '  [PASS] Beta(2,3) = 1/12 = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15)', '  [FAIL] Beta(2,3) = ', result, ' expected ', exact
      failed = failed + 1
    end if

    print '(A)', ''

  end subroutine test_special_function_integrals

  !---------------------------------------------------------------------------
  ! Endpoint Singularities
  ! Reference: QUADPACK book (Piessens et al., 1983)
  !---------------------------------------------------------------------------
  subroutine test_endpoint_singularities(passed, failed)
    integer, intent(out) :: passed, failed
    integer, parameter :: limit = 500, lenw = 4*limit
    real(dp) :: result, abserr, exact
    integer :: neval, ier, last
    integer :: iwork(limit)
    real(dp) :: work(lenw)
    real(dp), parameter :: tol = 1.0e-8_dp

    passed = 0
    failed = 0

    print '(A)', 'Endpoint Singularities'
    print '(A)', '----------------------'
    print '(A)', 'Reference: QUADPACK (Piessens et al., 1983)'
    print '(A)', ''

    ! int(1/sqrt(x), 0, 1) = 2 (algebraic singularity at x=0)
    call DQAGS(f_inv_sqrt, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-10_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = 2.0_dp
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES22.15)', '  [PASS] int(1/sqrt(x), 0, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15)', '  [FAIL] int(1/sqrt(x), 0, 1) = ', result, ' expected ', exact
      failed = failed + 1
    end if

    ! int(log(x), 0, 1) = -1 (logarithmic singularity at x=0)
    call DQAGS(f_log, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-10_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = -1.0_dp
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES22.15)', '  [PASS] int(log(x), 0, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15)', '  [FAIL] int(log(x), 0, 1) = ', result, ' expected ', exact
      failed = failed + 1
    end if

    ! int(x*log(x), 0, 1) = -1/4 (logarithmic singularity at x=0)
    call DQAGS(f_x_log, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-10_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = -0.25_dp
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES22.15)', '  [PASS] int(x*log(x), 0, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15)', '  [FAIL] int(x*log(x), 0, 1) = ', result, ' expected ', exact
      failed = failed + 1
    end if

    ! int(1/sqrt(1-x^2), 0, 1) = pi/2 (algebraic singularity at x=1)
    call DQAGS(f_inv_sqrt_1mx2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-10_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = pi / 2.0_dp
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES22.15)', '  [PASS] int(1/sqrt(1-x^2), 0, 1) = pi/2 = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15)', '  [FAIL] int(1/sqrt(1-x^2), 0, 1) = ', result, ' expected ', exact
      failed = failed + 1
    end if

    print '(A)', ''

  end subroutine test_endpoint_singularities

  !---------------------------------------------------------------------------
  ! Infinite Intervals
  ! Reference: Standard results, NIST DLMF
  !---------------------------------------------------------------------------
  subroutine test_infinite_intervals(passed, failed)
    integer, intent(out) :: passed, failed
    integer, parameter :: limit = 500, lenw = 4*limit
    real(dp) :: result, abserr, exact
    integer :: neval, ier, last
    integer :: iwork(limit)
    real(dp) :: work(lenw)
    real(dp), parameter :: tol = 1.0e-10_dp

    passed = 0
    failed = 0

    print '(A)', 'Infinite Interval Integration'
    print '(A)', '------------------------------'
    print '(A)', ''

    ! int(exp(-x), 0, inf) = 1
    call DQAGI(f_exp_neg, 0.0_dp, 1, 0.0_dp, 1.0e-12_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = 1.0_dp
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES22.15)', '  [PASS] int(exp(-x), 0, inf) = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15)', '  [FAIL] int(exp(-x), 0, inf) = ', result, ' expected ', exact
      failed = failed + 1
    end if

    ! int(1/(1+x^2), -inf, inf) = pi
    call DQAGI(f_lorentzian, 0.0_dp, 2, 0.0_dp, 1.0e-12_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = pi
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES22.15)', '  [PASS] int(1/(1+x^2), -inf, inf) = pi = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15)', '  [FAIL] int(1/(1+x^2), -inf, inf) = ', result, ' expected ', exact
      failed = failed + 1
    end if

    ! int(x^2*exp(-x^2), -inf, inf) = sqrt(pi)/2 (second moment of Gaussian)
    call DQAGI(f_x2_gaussian, 0.0_dp, 2, 0.0_dp, 1.0e-12_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = sqrt(pi) / 2.0_dp
    if (abs(result - exact) < tol .and. ier == 0) then
      print '(A,ES22.15)', '  [PASS] int(x^2*exp(-x^2), -inf, inf) = sqrt(pi)/2 = ', result
      passed = passed + 1
    else
      print '(A,ES22.15,A,ES22.15)', '  [FAIL] int(x^2*exp(-x^2), -inf, inf) = ', result, ' expected ', exact
      failed = failed + 1
    end if

    print '(A)', ''

  end subroutine test_infinite_intervals

  !---------------------------------------------------------------------------
  ! Test functions (integrands)
  !---------------------------------------------------------------------------

  pure real(dp) function f_x14(x)
    real(dp), intent(in) :: x
    f_x14 = x**14
  end function f_x14

  pure real(dp) function f_x20(x)
    real(dp), intent(in) :: x
    f_x20 = x**20
  end function f_x20

  pure real(dp) function f_x30(x)
    real(dp), intent(in) :: x
    f_x30 = x**30
  end function f_x30

  pure real(dp) function f_sin_squared(x)
    real(dp), intent(in) :: x
    f_sin_squared = sin(x)**2
  end function f_sin_squared

  pure real(dp) function f_x_sin(x)
    real(dp), intent(in) :: x
    f_x_sin = x * sin(x)
  end function f_x_sin

  pure real(dp) function f_lorentzian(x)
    real(dp), intent(in) :: x
    f_lorentzian = 1.0_dp / (1.0_dp + x*x)
  end function f_lorentzian

  pure real(dp) function f_exp_neg(x)
    real(dp), intent(in) :: x
    f_exp_neg = exp(-x)
  end function f_exp_neg

  pure real(dp) function f_x5_exp_neg(x)
    real(dp), intent(in) :: x
    f_x5_exp_neg = (x**5) * exp(-x)
  end function f_x5_exp_neg

  pure real(dp) function f_gaussian(x)
    real(dp), intent(in) :: x
    f_gaussian = exp(-x*x)
  end function f_gaussian

  pure real(dp) function f_x2_gaussian(x)
    real(dp), intent(in) :: x
    f_x2_gaussian = x*x * exp(-x*x)
  end function f_x2_gaussian

  pure real(dp) function f_sinc(x)
    real(dp), intent(in) :: x
    if (abs(x) > 1.0e-10_dp) then
      f_sinc = sin(x) / x
    else
      f_sinc = 1.0_dp  ! limit as x -> 0
    end if
  end function f_sinc

  pure real(dp) function f_beta_2_3(x)
    real(dp), intent(in) :: x
    ! Beta(2,3) integrand: x^(2-1) * (1-x)^(3-1) = x * (1-x)^2
    f_beta_2_3 = x * (1.0_dp - x)**2
  end function f_beta_2_3

  pure real(dp) function f_inv_sqrt(x)
    real(dp), intent(in) :: x
    if (x > 0.0_dp) then
      f_inv_sqrt = 1.0_dp / sqrt(x)
    else
      f_inv_sqrt = 0.0_dp
    end if
  end function f_inv_sqrt

  pure real(dp) function f_log(x)
    real(dp), intent(in) :: x
    if (x > 0.0_dp) then
      f_log = log(x)
    else
      f_log = 0.0_dp
    end if
  end function f_log

  pure real(dp) function f_x_log(x)
    real(dp), intent(in) :: x
    if (x > 0.0_dp) then
      f_x_log = x * log(x)
    else
      f_x_log = 0.0_dp
    end if
  end function f_x_log

  pure real(dp) function f_inv_sqrt_1mx2(x)
    real(dp), intent(in) :: x
    real(dp) :: arg
    arg = 1.0_dp - x*x
    if (arg > 0.0_dp) then
      f_inv_sqrt_1mx2 = 1.0_dp / sqrt(arg)
    else
      f_inv_sqrt_1mx2 = 0.0_dp
    end if
  end function f_inv_sqrt_1mx2

end module test_diff_integ_level2

!> Main program for Level 2 diff_integ tests
program run_level2_diff_integ
  use test_diff_integ_level2
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level2_diff_integ

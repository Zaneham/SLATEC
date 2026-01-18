!> Level 1: Regression Tests for Interpolation
!>
!> Purpose: Does the code work? Does it still work after changes?
!> These tests verify basic functionality and catch regressions.
!>
!> What we test:
!>   - DBVALU/BVALU: B-spline evaluation
!>   - DBINT4/BINT4: B-spline interpolation with cubic splines
!>   - DPCHIM/PCHIM: Piecewise cubic Hermite interpolation (monotone)
!>   - DPCHFE/PCHFE: PCHIP function evaluation
!>   - DPCHIA/PCHIA: PCHIP integration
!>   - DPLINT/POLINT: Polynomial interpolation
!>   - DPOLVL/POLYVL: Polynomial evaluation with derivatives
!>
!> Reference: SLATEC Library documentation

module test_interpolation_level1
  use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: tol_dp = 1.0e-10_dp
  real(sp), parameter :: tol_sp = 1.0e-5_sp
  real(dp), parameter :: pi = 3.14159265358979324_dp

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 1: INTERPOLATION REGRESSION TESTS'
    print '(A)', '================================================================'
    print '(A)', ''

    ! B-spline tests
    call test_bspline_suite(p, f)
    passed = passed + p
    failed = failed + f

    ! PCHIP tests
    call test_pchip_suite(p, f)
    passed = passed + p
    failed = failed + f

    ! Polynomial interpolation tests
    call test_polynomial_suite(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 1 SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! B-spline Test Suite
  !---------------------------------------------------------------------------
  subroutine test_bspline_suite(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'B-SPLINE INTERPOLATION (DBINT4, DBVALU)'
    print '(A)', '---------------------------------------'

    call test_dbint4_runs(passed, failed)
    call test_dbvalu_endpoints(passed, failed)
    call test_dbvalu_derivative(passed, failed)

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_bspline_suite

  subroutine test_dbint4_runs(passed, failed)
    !> Test that DBINT4 runs without error and returns valid output
    use interpolation, only: DBINT4, DBVALU
    integer, intent(inout) :: passed, failed

    integer, parameter :: ndata = 11
    integer, parameter :: ncoef = ndata + 2
    real(dp) :: x(ndata), y(ndata), t(ncoef+4), bc(ncoef), w(5*ncoef)
    real(dp) :: bv
    integer :: i, n, k

    ! Set up data: sin(pi*x) on [0,1]
    do i = 1, ndata
      x(i) = real(i-1, dp) / real(ndata-1, dp)
      y(i) = sin(pi * x(i))
    end do

    ! Call DBINT4 - should complete without error
    call DBINT4(x, y, ndata, 1, 2, pi, 0.0_dp, 1, t, bc, n, k, w)

    ! Check that n and k have expected values
    if (n == ncoef .and. k == 4) then
      ! Check that we can evaluate at first and last points
      bv = DBVALU(t, bc, n, k, 0, x(1))
      if (abs(y(1) - bv) < 1.0e-10_dp) then
        bv = DBVALU(t, bc, n, k, 0, x(ndata))
        if (abs(y(ndata) - bv) < 1.0e-10_dp) then
          print '(A)', '  [PASS] DBINT4: runs and returns valid output'
          passed = passed + 1
          return
        end if
      end if
    end if

    print '(A)', '  [FAIL] DBINT4: unexpected output'
    failed = failed + 1
  end subroutine

  subroutine test_dbvalu_endpoints(passed, failed)
    !> Test DBVALU at endpoints reproduces boundary data values
    use interpolation, only: DBINT4, DBVALU
    integer, intent(inout) :: passed, failed

    integer, parameter :: ndata = 5
    integer, parameter :: ncoef = ndata + 2
    real(dp) :: x(ndata), y(ndata), t(ncoef+4), bc(ncoef), w(5*ncoef)
    real(dp) :: bv_left, bv_right, err_left, err_right
    integer :: n, k, i

    ! Linear function: y = 2x + 1
    do i = 1, ndata
      x(i) = real(i-1, dp)
      y(i) = 2.0_dp * x(i) + 1.0_dp
    end do

    call DBINT4(x, y, ndata, 1, 1, 2.0_dp, 2.0_dp, 1, t, bc, n, k, w)

    bv_left = DBVALU(t, bc, n, k, 0, x(1))
    bv_right = DBVALU(t, bc, n, k, 0, x(ndata))

    err_left = abs(y(1) - bv_left)
    err_right = abs(y(ndata) - bv_right)

    ! Endpoints should be exact
    if (err_left < 1.0e-10_dp .and. err_right < 1.0e-10_dp) then
      print '(A)', '  [PASS] DBVALU: endpoint values match'
      passed = passed + 1
    else
      print '(A,2ES12.5)', '  [FAIL] DBVALU: endpoint errors = ', err_left, err_right
      failed = failed + 1
    end if
  end subroutine

  subroutine test_dbvalu_derivative(passed, failed)
    !> Test DBVALU can compute derivatives (returns finite value)
    use interpolation, only: DBINT4, DBVALU
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer, intent(inout) :: passed, failed

    integer, parameter :: ndata = 5
    integer, parameter :: ncoef = ndata + 2
    real(dp) :: x(ndata), y(ndata), t(ncoef+4), bc(ncoef), w(5*ncoef)
    real(dp) :: deriv
    integer :: n, k, i

    ! Quadratic function: y = x^2
    do i = 1, ndata
      x(i) = real(i-1, dp)
      y(i) = x(i)**2
    end do

    ! Use second derivative BC (natural spline)
    call DBINT4(x, y, ndata, 2, 2, 0.0_dp, 0.0_dp, 1, t, bc, n, k, w)

    ! Get first derivative - just check it returns a finite value
    deriv = DBVALU(t, bc, n, k, 1, 2.0_dp)

    if (ieee_is_finite(deriv)) then
      print '(A)', '  [PASS] DBVALU: derivative returns finite value'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DBVALU: derivative returned non-finite'
      failed = failed + 1
    end if
  end subroutine

  !---------------------------------------------------------------------------
  ! PCHIP Test Suite
  !---------------------------------------------------------------------------
  subroutine test_pchip_suite(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'PCHIP INTERPOLATION (DPCHIM, DPCHFE, DPCHIA)'
    print '(A)', '--------------------------------------------'

    call test_dpchim_monotone(passed, failed)
    call test_dpchfe_evaluation(passed, failed)
    call test_dpchia_integration(passed, failed)

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_pchip_suite

  subroutine test_dpchim_monotone(passed, failed)
    !> Test DPCHIM preserves monotonicity
    use interpolation, only: DPCHIM
    integer, intent(inout) :: passed, failed

    integer, parameter :: n = 5
    integer, parameter :: incfd = 1
    real(dp) :: x(n), f(incfd,n), d(incfd,n)
    integer :: ierr, i
    logical :: all_nonneg

    ! Strictly increasing function: x^2
    x = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    do i = 1, n
      f(1,i) = x(i)**2
    end do

    call DPCHIM(n, x, f, d, incfd, ierr)

    ! For monotone increasing data, derivatives should be non-negative
    all_nonneg = all(d(1,:) >= -1.0e-10_dp)

    if (ierr == 0 .and. all_nonneg) then
      print '(A)', '  [PASS] DPCHIM: monotonicity preserved'
      passed = passed + 1
    else
      print '(A,I3)', '  [FAIL] DPCHIM: ierr = ', ierr
      failed = failed + 1
    end if
  end subroutine

  subroutine test_dpchfe_evaluation(passed, failed)
    !> Test DPCHFE evaluates correctly at data points
    use interpolation, only: DPCHIM, DPCHFE
    integer, intent(inout) :: passed, failed

    integer, parameter :: n = 5, ne = 5
    integer, parameter :: incfd = 1
    real(dp) :: x(n), f(incfd,n), d(incfd,n), xe(ne), fe(ne)
    integer :: ierr, i
    logical :: skip
    real(dp) :: max_err, err

    ! Set up cubic: f(x) = x^3
    x = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    do i = 1, n
      f(1,i) = x(i)**3
    end do

    call DPCHIM(n, x, f, d, incfd, ierr)

    ! Evaluate at data points
    xe = x
    skip = .false.
    call DPCHFE(n, x, f, d, incfd, skip, ne, xe, fe, ierr)

    max_err = 0.0_dp
    do i = 1, ne
      err = abs(f(1,i) - fe(i))
      max_err = max(max_err, err)
    end do

    if (ierr == 0 .and. max_err < 1.0e-10_dp) then
      print '(A)', '  [PASS] DPCHFE: evaluation at data points'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [FAIL] DPCHFE: max error = ', max_err
      failed = failed + 1
    end if
  end subroutine

  subroutine test_dpchia_integration(passed, failed)
    !> Test DPCHIA integration of linear function
    use interpolation, only: DPCHIM, DPCHIA
    integer, intent(inout) :: passed, failed

    integer, parameter :: n = 3
    integer, parameter :: incfd = 1
    real(dp) :: x(n), f(incfd,n), d(incfd,n)
    real(dp) :: integral, exact, err
    integer :: ierr, i

    ! Linear function: f(x) = 2x on [0, 2]
    ! Integral from 0 to 2 = x^2 |_0^2 = 4
    x = [0.0_dp, 1.0_dp, 2.0_dp]
    do i = 1, n
      f(1,i) = 2.0_dp * x(i)
    end do

    call DPCHIM(n, x, f, d, incfd, ierr)

    integral = DPCHIA(n, x, f, d, incfd, 0.0_dp, 2.0_dp)
    exact = 4.0_dp
    err = abs(exact - integral)

    if (err < 1.0e-10_dp) then
      print '(A)', '  [PASS] DPCHIA: linear function integration'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [FAIL] DPCHIA: error = ', err
      failed = failed + 1
    end if
  end subroutine

  !---------------------------------------------------------------------------
  ! Polynomial Interpolation Test Suite
  !---------------------------------------------------------------------------
  subroutine test_polynomial_suite(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'POLYNOMIAL INTERPOLATION (DPLINT, DPOLCF, DPOLVL)'
    print '(A)', '-------------------------------------------------'

    call test_dplint_basic(passed, failed)
    call test_dpolcf_taylor(passed, failed)
    call test_dpolvl_derivatives(passed, failed)

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_polynomial_suite

  subroutine test_dplint_basic(passed, failed)
    !> Test DPLINT basic interpolation
    use interpolation, only: DPLINT, DPOLVL
    integer, intent(inout) :: passed, failed

    integer, parameter :: n = 3
    real(dp) :: x(n), y(n), c(n), d(n), w(2*n)
    real(dp) :: yf, err
    integer :: ierr

    ! Quadratic through (0,0), (1,1), (2,4) => y = x^2
    x = [0.0_dp, 1.0_dp, 2.0_dp]
    y = [0.0_dp, 1.0_dp, 4.0_dp]

    call DPLINT(n, x, y, c)

    ! Evaluate at x = 1.5, expect 2.25
    call DPOLVL(0, 1.5_dp, yf, d, n, x, c, w, ierr)
    err = abs(2.25_dp - yf)

    if (err < tol_dp) then
      print '(A)', '  [PASS] DPLINT: quadratic interpolation'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [FAIL] DPLINT: error at x=1.5 = ', err
      failed = failed + 1
    end if
  end subroutine

  subroutine test_dpolcf_taylor(passed, failed)
    !> Test DPOLCF Taylor coefficients
    use interpolation, only: DPLINT, DPOLCF
    integer, intent(inout) :: passed, failed

    integer, parameter :: n = 5
    real(dp) :: x(n), y(n), c(n), d(n), w(2*n)
    real(dp) :: err
    integer :: i

    ! Quartic: f(x) = x^4 - 2x^2 + 1 = (x^2-1)^2
    ! Taylor at x=0: 1 + 0*x - 2*x^2 + 0*x^3 + 1*x^4
    real(dp), parameter :: expected(n) = [1.0_dp, 0.0_dp, -2.0_dp, 0.0_dp, 1.0_dp]

    x = [real(dp) :: -2, -1, 0, 1, 2]
    do i = 1, n
      y(i) = (x(i)**2 - 1.0_dp)**2
    end do

    call DPLINT(n, x, y, c)
    call DPOLCF(0.0_dp, n, x, c, d, w)

    err = maxval(abs(expected - d))

    if (err < 1.0e-8_dp) then
      print '(A)', '  [PASS] DPOLCF: Taylor coefficients'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [FAIL] DPOLCF: max coefficient error = ', err
      failed = failed + 1
    end if
  end subroutine

  subroutine test_dpolvl_derivatives(passed, failed)
    !> Test DPOLVL derivative computation
    use interpolation, only: DPLINT, DPOLVL
    integer, intent(inout) :: passed, failed

    integer, parameter :: n = 4
    real(dp) :: x(n), y(n), c(n), d(n), w(2*n)
    real(dp) :: yf, err
    integer :: ierr

    ! Cubic: f(x) = x^3, f'(x) = 3x^2, f''(x) = 6x, f'''(x) = 6
    x = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    y = [0.0_dp, 1.0_dp, 8.0_dp, 27.0_dp]

    call DPLINT(n, x, y, c)

    ! Get all derivatives at x = 1
    call DPOLVL(3, 1.0_dp, yf, d, n, x, c, w, ierr)

    ! d(1) = f'(1) = 3, d(2) = f''(1) = 6, d(3) = f'''(1) = 6
    err = abs(3.0_dp - d(1)) + abs(6.0_dp - d(2)) + abs(6.0_dp - d(3))

    if (err < tol_dp) then
      print '(A)', '  [PASS] DPOLVL: derivatives of x^3 at x=1'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [FAIL] DPOLVL: derivative error = ', err
      failed = failed + 1
    end if
  end subroutine

end module test_interpolation_level1

!---------------------------------------------------------------------------
! Main program
!---------------------------------------------------------------------------
program test_l1_interpolation
  use test_interpolation_level1, only: run_all_tests
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed == 0) then
    print '(A)', ''
    print '(A)', 'SUCCESS: All Level 1 interpolation tests passed!'
    stop 0
  else
    print '(A)', ''
    print '(A,I3,A)', 'FAILURE: ', failed, ' test(s) failed'
    stop 1
  end if

end program test_l1_interpolation

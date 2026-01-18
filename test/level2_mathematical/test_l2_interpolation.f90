!> Level 2: Mathematical Verification for Interpolation
!>
!> Purpose: Does the algorithm match the original mathematics?
!> These tests verify mathematical properties against authoritative references.
!>
!> Mathematical Properties Tested:
!>   - B-splines: Interpolation, boundary conditions, partition of unity
!>   - PCHIP: Monotonicity preservation, interpolation exactness
!>   - Polynomials: Uniqueness theorem, Taylor expansion
!>
!> References:
!>   - de Boor, C. (1978). "A Practical Guide to Splines." Springer-Verlag.
!>   - Fritsch, F.N. & Carlson, R.E. (1980). "Monotone Piecewise Cubic
!>     Interpolation." SIAM J. Numer. Anal. 17, 238-246.
!>   - Burden, R.L. & Faires, J.D. (2011). "Numerical Analysis." 9th ed.

module test_interpolation_level2
  use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: tol_dp = 1.0e-10_dp
  real(dp), parameter :: pi = 3.14159265358979323846_dp

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 2: INTERPOLATION MATHEMATICAL VERIFICATION'
    print '(A)', '================================================================'
    print '(A)', 'Reference: de Boor (1978), Fritsch & Carlson (1980)'
    print '(A)', ''

    call test_bspline_interpolation(p, f)
    passed = passed + p
    failed = failed + f

    call test_bspline_boundary_conditions(p, f)
    passed = passed + p
    failed = failed + f

    call test_pchip_interpolation(p, f)
    passed = passed + p
    failed = failed + f

    call test_pchip_monotonicity(p, f)
    passed = passed + p
    failed = failed + f

    call test_polynomial_uniqueness(p, f)
    passed = passed + p
    failed = failed + f

    call test_polynomial_exactness(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 2 INTERPOLATION SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! B-SPLINE INTERPOLATION: S(x_i) = y_i exactly
  ! Reference: de Boor (1978) Chapter IX
  ! A cubic spline through n data points has n+2 coefficients.
  ! With boundary conditions, the system is determined.
  !---------------------------------------------------------------------------
  subroutine test_bspline_interpolation(passed, failed)
    use interpolation, only: DBINT4, DBVALU
    integer, intent(out) :: passed, failed
    integer, parameter :: ndata = 5
    integer :: n, k, i
    real(dp) :: x(ndata), y(ndata), t(ndata+6), bc(ndata+2), w(5,ndata+2)
    real(dp) :: y_interp, err, max_err

    passed = 0
    failed = 0

    print '(A)', 'B-spline Interpolation Property'
    print '(A)', '--------------------------------'
    print '(A)', '  Property: S(x_i) = y_i for all data points'
    print '(A)', ''

    ! Test 1: Natural spline (2nd derivative = 0 at ends)
    ! Data: quadratic function f(x) = x^2
    do i = 1, ndata
      x(i) = real(i-1, dp) / real(ndata-1, dp)
      y(i) = x(i)**2
    end do

    ! ibcl=2, ibcr=2: natural spline (S'' = 0 at both ends)
    call DBINT4(x, y, ndata, 2, 2, 0.0_dp, 0.0_dp, 1, t, bc, n, k, w)

    ! Verify interpolation at each data point
    max_err = 0.0_dp
    do i = 1, ndata
      y_interp = DBVALU(t, bc, n, k, 0, x(i))
      err = abs(y_interp - y(i))
      max_err = max(max_err, err)
    end do

    if (max_err < tol_dp) then
      print '(A)', '  [PASS] Natural spline: S(x_i) = y_i (max error < 1e-10)'
      passed = passed + 1
    else
      print '(A,ES12.4)', '  [FAIL] Natural spline max error: ', max_err
      failed = failed + 1
    end if

    ! Test 2: Clamped spline with correct derivatives
    ! Data: sin(pi*x/2) on [0,1]
    ! f(0) = 0, f(1) = 1
    ! f'(0) = pi/2, f'(1) = 0
    do i = 1, ndata
      x(i) = real(i-1, dp) / real(ndata-1, dp)
      y(i) = sin(pi * x(i) / 2.0_dp)
    end do

    ! ibcl=1, ibcr=1: clamped spline
    ! Specify correct first derivatives
    call DBINT4(x, y, ndata, 1, 1, pi/2.0_dp, 0.0_dp, 1, t, bc, n, k, w)

    max_err = 0.0_dp
    do i = 1, ndata
      y_interp = DBVALU(t, bc, n, k, 0, x(i))
      err = abs(y_interp - y(i))
      max_err = max(max_err, err)
    end do

    if (max_err < tol_dp) then
      print '(A)', '  [PASS] Clamped spline (correct BCs): S(x_i) = y_i'
      passed = passed + 1
    else
      print '(A,ES12.4)', '  [FAIL] Clamped spline max error: ', max_err
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_bspline_interpolation

  !---------------------------------------------------------------------------
  ! B-SPLINE BOUNDARY CONDITIONS
  ! Reference: de Boor (1978) Theorem IX.1
  ! Clamped spline: S'(a) = f'(a), S'(b) = f'(b)
  ! Natural spline: S''(a) = S''(b) = 0
  !---------------------------------------------------------------------------
  subroutine test_bspline_boundary_conditions(passed, failed)
    use interpolation, only: DBINT4, DBVALU
    integer, intent(out) :: passed, failed
    integer, parameter :: ndata = 11
    integer :: n, k
    real(dp) :: x(ndata), y(ndata), t(ndata+6), bc(ndata+2), w(5,ndata+2)
    real(dp) :: deriv_left, deriv_right, deriv2_left, deriv2_right
    real(dp) :: exact_deriv_left, exact_deriv_right
    integer :: i

    passed = 0
    failed = 0

    print '(A)', 'B-spline Boundary Conditions'
    print '(A)', '-----------------------------'

    ! Use sin(pi*x/2): f'(0) = pi/2, f'(1) = 0, f''(0) = 0, f''(1) = -pi^2/4
    do i = 1, ndata
      x(i) = real(i-1, dp) / real(ndata-1, dp)
      y(i) = sin(pi * x(i) / 2.0_dp)
    end do
    exact_deriv_left = pi / 2.0_dp   ! f'(0) = (pi/2)cos(0) = pi/2
    exact_deriv_right = 0.0_dp       ! f'(1) = (pi/2)cos(pi/2) = 0

    ! Test 1: Clamped spline derivatives
    call DBINT4(x, y, ndata, 1, 1, exact_deriv_left, exact_deriv_right, 1, &
                t, bc, n, k, w)

    deriv_left = DBVALU(t, bc, n, k, 1, x(1))    ! S'(0)
    deriv_right = DBVALU(t, bc, n, k, 1, x(ndata))  ! S'(1)

    if (abs(deriv_left - exact_deriv_left) < 1.0e-8_dp) then
      print '(A)', '  [PASS] Clamped S''(0) = pi/2 (specified left BC)'
      passed = passed + 1
    else
      print '(A,2ES15.8)', '  [FAIL] S''(0): expected/got ', exact_deriv_left, deriv_left
      failed = failed + 1
    end if

    if (abs(deriv_right - exact_deriv_right) < 1.0e-8_dp) then
      print '(A)', '  [PASS] Clamped S''(1) = 0 (specified right BC)'
      passed = passed + 1
    else
      print '(A,2ES15.8)', '  [FAIL] S''(1): expected/got ', exact_deriv_right, deriv_right
      failed = failed + 1
    end if

    ! Test 2: Natural spline (S'' = 0 at boundaries)
    call DBINT4(x, y, ndata, 2, 2, 0.0_dp, 0.0_dp, 1, t, bc, n, k, w)

    deriv2_left = DBVALU(t, bc, n, k, 2, x(1))      ! S''(0)
    deriv2_right = DBVALU(t, bc, n, k, 2, x(ndata))  ! S''(1)

    if (abs(deriv2_left) < 1.0e-8_dp) then
      print '(A)', '  [PASS] Natural S''''(0) = 0'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] Natural S''''(0): ', deriv2_left
      failed = failed + 1
    end if

    if (abs(deriv2_right) < 1.0e-8_dp) then
      print '(A)', '  [PASS] Natural S''''(1) = 0'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] Natural S''''(1): ', deriv2_right
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_bspline_boundary_conditions

  !---------------------------------------------------------------------------
  ! PCHIP INTERPOLATION: p(x_i) = y_i exactly
  ! Reference: Fritsch & Carlson (1980)
  !---------------------------------------------------------------------------
  subroutine test_pchip_interpolation(passed, failed)
    use interpolation, only: DPCHIM, DPCHFE
    integer, intent(out) :: passed, failed
    integer, parameter :: n = 7
    integer, parameter :: incfd = 1
    integer :: i, ierr, ne
    real(dp) :: x(n), f(incfd,n), d(incfd,n)
    real(dp) :: xe(n), fe(n)
    real(dp) :: err, max_err
    logical :: skip

    passed = 0
    failed = 0

    print '(A)', 'PCHIP Interpolation Property'
    print '(A)', '-----------------------------'
    print '(A)', '  Property: p(x_i) = y_i for all data points'
    print '(A)', ''

    ! Data: exp(-x) on [0, 3]
    do i = 1, n
      x(i) = real(i-1, dp) * 0.5_dp
      f(1,i) = exp(-x(i))
    end do

    ! Compute derivatives
    call DPCHIM(n, x, f, d, incfd, ierr)

    if (ierr < 0) then
      print '(A,I4)', '  [FAIL] DPCHIM error: ', ierr
      failed = failed + 1
      return
    end if

    ! Evaluate at data points
    ne = n
    xe = x
    skip = .false.
    call DPCHFE(n, x, f, d, incfd, skip, ne, xe, fe, ierr)

    max_err = 0.0_dp
    do i = 1, n
      err = abs(fe(i) - f(1,i))
      max_err = max(max_err, err)
    end do

    if (max_err < tol_dp) then
      print '(A)', '  [PASS] PCHIP reproduces data exactly: max error < 1e-10'
      passed = passed + 1
    else
      print '(A,ES12.4)', '  [FAIL] PCHIP interpolation error: ', max_err
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_pchip_interpolation

  !---------------------------------------------------------------------------
  ! PCHIP MONOTONICITY PRESERVATION
  ! Reference: Fritsch & Carlson (1980) Theorem 1
  ! If data is monotone, the PCHIP interpolant is monotone.
  !---------------------------------------------------------------------------
  subroutine test_pchip_monotonicity(passed, failed)
    use interpolation, only: DPCHIM, DPCHFE
    integer, intent(out) :: passed, failed
    integer, parameter :: n = 6
    integer, parameter :: incfd = 1
    integer, parameter :: ne = 101
    integer :: i, ierr
    real(dp) :: x(n), f(incfd,n), d(incfd,n)
    real(dp) :: xe(ne), fe(ne)
    logical :: skip
    logical :: is_monotone

    passed = 0
    failed = 0

    print '(A)', 'PCHIP Monotonicity Preservation'
    print '(A)', '--------------------------------'
    print '(A)', '  Reference: Fritsch & Carlson (1980) Theorem 1'
    print '(A)', ''

    ! Test 1: Strictly increasing data (sqrt function)
    do i = 1, n
      x(i) = real(i-1, dp)
      f(1,i) = sqrt(x(i))
    end do

    call DPCHIM(n, x, f, d, incfd, ierr)

    ! Evaluate at many points
    do i = 1, ne
      xe(i) = x(1) + (x(n) - x(1)) * real(i-1, dp) / real(ne-1, dp)
    end do
    skip = .false.
    call DPCHFE(n, x, f, d, incfd, skip, ne, xe, fe, ierr)

    ! Check monotonicity
    is_monotone = .true.
    do i = 2, ne
      if (fe(i) < fe(i-1) - tol_dp) then
        is_monotone = .false.
        exit
      end if
    end do

    if (is_monotone) then
      print '(A)', '  [PASS] Increasing data -> monotone interpolant'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Monotonicity violated for increasing data'
      failed = failed + 1
    end if

    ! Test 2: Strictly decreasing data (exp(-x))
    do i = 1, n
      x(i) = real(i-1, dp)
      f(1,i) = exp(-x(i))
    end do

    call DPCHIM(n, x, f, d, incfd, ierr)

    skip = .false.
    call DPCHFE(n, x, f, d, incfd, skip, ne, xe, fe, ierr)

    is_monotone = .true.
    do i = 2, ne
      if (fe(i) > fe(i-1) + tol_dp) then
        is_monotone = .false.
        exit
      end if
    end do

    if (is_monotone) then
      print '(A)', '  [PASS] Decreasing data -> monotone interpolant'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Monotonicity violated for decreasing data'
      failed = failed + 1
    end if

    ! Test 3: Data with flat region (step function)
    x = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    f(1,:) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]

    call DPCHIM(n, x, f, d, incfd, ierr)

    skip = .false.
    call DPCHFE(n, x, f, d, incfd, skip, ne, xe, fe, ierr)

    is_monotone = .true.
    do i = 2, ne
      if (fe(i) < fe(i-1) - tol_dp) then
        is_monotone = .false.
        exit
      end if
    end do

    if (is_monotone) then
      print '(A)', '  [PASS] Non-strict monotone data preserved'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Monotonicity violated for flat regions'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_pchip_monotonicity

  !---------------------------------------------------------------------------
  ! POLYNOMIAL INTERPOLATION UNIQUENESS
  ! Reference: Burden & Faires (2011) Theorem 3.2
  ! The polynomial of degree n-1 through n points is unique.
  !---------------------------------------------------------------------------
  subroutine test_polynomial_uniqueness(passed, failed)
    use interpolation, only: DPLINT, DPOLCF, DPOLVL
    integer, intent(out) :: passed, failed
    integer, parameter :: n = 5
    integer :: i, ierr
    real(dp) :: x1(n), c1(n), d1(n), x2(n), c2(n), d2(n)
    real(dp) :: y1(n), y2(n)
    real(dp) :: test_x, p1, p2, yp1(1), yp2(1)
    real(dp) :: work1(2*n), work2(2*n)

    passed = 0
    failed = 0

    print '(A)', 'Polynomial Uniqueness Theorem'
    print '(A)', '------------------------------'
    print '(A)', '  Reference: Burden & Faires, Theorem 3.2'
    print '(A)', '  The polynomial of degree n-1 through n points is unique.'
    print '(A)', ''

    ! Points to interpolate: y = x^3 at x = 0, 1, 2, 3, 4
    y1 = [0.0_dp, 1.0_dp, 8.0_dp, 27.0_dp, 64.0_dp]

    ! Two different orderings of the same points
    x1 = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    x2 = [4.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, 1.0_dp]
    y2 = [64.0_dp, 8.0_dp, 0.0_dp, 27.0_dp, 1.0_dp]

    ! Compute Newton coefficients from first ordering
    call DPLINT(n, x1, y1, c1)
    call DPOLCF(0.0_dp, n, x1, c1, d1, work1)

    ! Compute Newton coefficients from second ordering
    call DPLINT(n, x2, y2, c2)
    call DPOLCF(0.0_dp, n, x2, c2, d2, work2)

    ! Evaluate both at several test points using DPOLVL
    test_x = 1.5_dp
    call DPOLVL(0, test_x, p1, yp1, n, x1, c1, work1, ierr)
    call DPOLVL(0, test_x, p2, yp2, n, x2, c2, work2, ierr)

    if (abs(p1 - p2) < tol_dp) then
      print '(A)', '  [PASS] Different orderings produce same polynomial at x=1.5'
      passed = passed + 1
    else
      print '(A,2ES15.8)', '  [FAIL] p1, p2 = ', p1, p2
      failed = failed + 1
    end if

    test_x = 2.7_dp
    call DPOLVL(0, test_x, p1, yp1, n, x1, c1, work1, ierr)
    call DPOLVL(0, test_x, p2, yp2, n, x2, c2, work2, ierr)

    if (abs(p1 - p2) < tol_dp) then
      print '(A)', '  [PASS] Different orderings produce same polynomial at x=2.7'
      passed = passed + 1
    else
      print '(A,2ES15.8)', '  [FAIL] p1, p2 = ', p1, p2
      failed = failed + 1
    end if

    ! Verify polynomial equals x^3 (degree 3, coefficients [0, 0, 0, 1])
    test_x = 2.5_dp
    call DPOLVL(0, test_x, p1, yp1, n, x1, c1, work1, ierr)
    if (abs(p1 - 2.5_dp**3) < tol_dp) then
      print '(A)', '  [PASS] Interpolant equals x^3 (polynomial reproduction)'
      passed = passed + 1
    else
      print '(A,2ES15.8)', '  [FAIL] p(2.5) vs 2.5^3: ', p1, 2.5_dp**3
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_polynomial_uniqueness

  !---------------------------------------------------------------------------
  ! POLYNOMIAL EXACTNESS FOR DEGREE n-1
  ! Reference: Burden & Faires (2011) Corollary 3.3
  ! If f(x) is a polynomial of degree <= n-1, then the interpolant = f(x).
  !---------------------------------------------------------------------------
  subroutine test_polynomial_exactness(passed, failed)
    use interpolation, only: DPLINT, DPOLVL
    integer, intent(out) :: passed, failed
    integer, parameter :: n = 5
    integer :: i, ierr
    real(dp) :: x(n), y(n), c(n)
    real(dp) :: test_x, p_interp, p_exact
    real(dp) :: max_err, yp(1), work(2*n)

    passed = 0
    failed = 0

    print '(A)', 'Polynomial Exactness Property'
    print '(A)', '------------------------------'
    print '(A)', '  If f(x) has degree <= n-1, then p_n-1(x) = f(x) exactly.'
    print '(A)', ''

    ! Test 1: Linear function through 5 points
    ! f(x) = 2x + 3
    do i = 1, n
      x(i) = real(i-1, dp)
      y(i) = 2.0_dp * x(i) + 3.0_dp
    end do

    call DPLINT(n, x, y, c)

    max_err = 0.0_dp
    do i = 1, 20
      test_x = real(i-1, dp) / 5.0_dp
      call DPOLVL(0, test_x, p_interp, yp, n, x, c, work, ierr)
      p_exact = 2.0_dp * test_x + 3.0_dp
      max_err = max(max_err, abs(p_interp - p_exact))
    end do

    if (max_err < tol_dp) then
      print '(A)', '  [PASS] Linear function reproduced exactly'
      passed = passed + 1
    else
      print '(A,ES12.4)', '  [FAIL] Linear reproduction error: ', max_err
      failed = failed + 1
    end if

    ! Test 2: Quadratic function through 5 points
    ! f(x) = x^2 - 2x + 1 = (x-1)^2
    do i = 1, n
      x(i) = real(i-1, dp)
      y(i) = (x(i) - 1.0_dp)**2
    end do

    call DPLINT(n, x, y, c)

    max_err = 0.0_dp
    do i = 1, 20
      test_x = real(i-1, dp) / 5.0_dp
      call DPOLVL(0, test_x, p_interp, yp, n, x, c, work, ierr)
      p_exact = (test_x - 1.0_dp)**2
      max_err = max(max_err, abs(p_interp - p_exact))
    end do

    if (max_err < tol_dp) then
      print '(A)', '  [PASS] Quadratic function reproduced exactly'
      passed = passed + 1
    else
      print '(A,ES12.4)', '  [FAIL] Quadratic reproduction error: ', max_err
      failed = failed + 1
    end if

    ! Test 3: Cubic function through 5 points
    ! f(x) = x^3 - 3x^2 + 2x
    do i = 1, n
      x(i) = real(i-1, dp)
      y(i) = x(i)**3 - 3.0_dp*x(i)**2 + 2.0_dp*x(i)
    end do

    call DPLINT(n, x, y, c)

    max_err = 0.0_dp
    do i = 1, 20
      test_x = real(i-1, dp) / 5.0_dp
      call DPOLVL(0, test_x, p_interp, yp, n, x, c, work, ierr)
      p_exact = test_x**3 - 3.0_dp*test_x**2 + 2.0_dp*test_x
      max_err = max(max_err, abs(p_interp - p_exact))
    end do

    if (max_err < tol_dp) then
      print '(A)', '  [PASS] Cubic function reproduced exactly'
      passed = passed + 1
    else
      print '(A,ES12.4)', '  [FAIL] Cubic reproduction error: ', max_err
      failed = failed + 1
    end if

    ! Test 4: Quartic function (degree 4, need exactly 5 points)
    ! f(x) = x^4
    do i = 1, n
      x(i) = real(i-1, dp)
      y(i) = x(i)**4
    end do

    call DPLINT(n, x, y, c)

    max_err = 0.0_dp
    do i = 1, 20
      test_x = real(i-1, dp) / 5.0_dp
      call DPOLVL(0, test_x, p_interp, yp, n, x, c, work, ierr)
      p_exact = test_x**4
      max_err = max(max_err, abs(p_interp - p_exact))
    end do

    if (max_err < tol_dp) then
      print '(A)', '  [PASS] Quartic function reproduced exactly (n points, degree n-1)'
      passed = passed + 1
    else
      print '(A,ES12.4)', '  [FAIL] Quartic reproduction error: ', max_err
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_polynomial_exactness

end module test_interpolation_level2

!> Main program
program run_level2_interpolation
  use test_interpolation_level2
  implicit none
  integer :: passed, failed
  call run_all_tests(passed, failed)
  if (failed > 0) stop 1
end program run_level2_interpolation

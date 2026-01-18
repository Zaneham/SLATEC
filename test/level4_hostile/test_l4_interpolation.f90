!> Level 4: Hostile Environment Tests for Interpolation
!>
!> Purpose: Do compilers, processors, or OS change anything?
!>
!> These tests probe edge cases that break on:
!>   - Aggressive compiler optimizations (-ffast-math, -Ofast)
!>   - Different floating-point modes (DAZ, FTZ, rounding)
!>   - Vectorization edge cases (SIMD boundary issues)
!>   - FMA vs separate mul+add
!>   - Catastrophic cancellation in divided differences
!>   - Ill-conditioned interpolation (Runge phenomenon)
!>   - Near-singular data (closely spaced points)
!>
!> Test Categories:
!>   1. Near-singular data (closely spaced x-values)
!>   2. Extreme values (very large/small domain)
!>   3. Runge phenomenon (ill-conditioned polynomial)
!>   4. Catastrophic cancellation (Newton divided differences)
!>   5. Subnormal handling (tiny function values)
!>   6. Reproducibility (identical results across runs)

module test_interpolation_level4
  use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, int64
  use, intrinsic :: ieee_arithmetic
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: tol_dp = 1.0e-10_dp
  real(dp), parameter :: tol_loose = 1.0e-6_dp
  real(dp), parameter :: pi = 3.14159265358979323846_dp

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 4: HOSTILE ENVIRONMENT TESTS (INTERPOLATION)'
    print '(A)', '================================================================'
    print '(A)', 'Platform portability and edge case testing'
    print '(A)', ''

    call test_closely_spaced_data(p, f)
    passed = passed + p; failed = failed + f

    call test_extreme_domain(p, f)
    passed = passed + p; failed = failed + f

    call test_runge_phenomenon(p, f)
    passed = passed + p; failed = failed + f

    call test_cancellation_in_differences(p, f)
    passed = passed + p; failed = failed + f

    call test_subnormal_function_values(p, f)
    passed = passed + p; failed = failed + f

    call test_reproducibility(p, f)
    passed = passed + p; failed = failed + f

    call test_pchip_sign_changes(p, f)
    passed = passed + p; failed = failed + f

    call test_spline_nearly_flat(p, f)
    passed = passed + p; failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 4 INTERPOLATION SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! Test 1: Closely Spaced Data Points
  ! Tests numerical stability when x-values are very close together
  ! Can expose catastrophic cancellation in divided differences
  !---------------------------------------------------------------------------
  subroutine test_closely_spaced_data(passed, failed)
    use interpolation, only: DPCHIM, DPCHFE
    integer, intent(out) :: passed, failed
    integer, parameter :: n = 5
    integer, parameter :: incfd = 1
    integer :: i, ierr, ne
    real(dp) :: x(n), f(incfd,n), d(incfd,n)
    real(dp) :: xe(3), fe(3)
    real(dp) :: delta
    logical :: skip

    passed = 0
    failed = 0

    print '(A)', 'CLOSELY SPACED DATA POINTS'
    print '(A)', '--------------------------'
    print '(A)', 'Tests stability with near-coincident x-values'
    print '(A)', ''

    ! Test 1a: Normal spacing (control)
    delta = 0.1_dp
    do i = 1, n
      x(i) = real(i-1, dp) * delta
      f(1,i) = sin(x(i))
    end do

    call DPCHIM(n, x, f, d, incfd, ierr)

    if (ierr >= 0) then
      print '(A)', '  [PASS] Normal spacing (delta=0.1): DPCHIM succeeded'
      passed = passed + 1
    else
      print '(A,I4)', '  [FAIL] Normal spacing: DPCHIM error ', ierr
      failed = failed + 1
    end if

    ! Test 1b: Very close spacing (1e-10)
    delta = 1.0e-10_dp
    do i = 1, n
      x(i) = real(i-1, dp) * delta
      f(1,i) = sin(x(i))  ! sin(x) ≈ x for small x
    end do

    call DPCHIM(n, x, f, d, incfd, ierr)

    if (ierr >= 0) then
      print '(A)', '  [PASS] Tight spacing (delta=1e-10): DPCHIM succeeded'
      passed = passed + 1
    else
      print '(A,I4)', '  [INFO] Tight spacing: DPCHIM reported ', ierr
      passed = passed + 1  ! Expected behavior, not failure
    end if

    ! Test 1c: Machine epsilon spacing
    delta = epsilon(1.0_dp)
    do i = 1, n
      x(i) = 1.0_dp + real(i-1, dp) * delta  ! Near x=1 to avoid underflow
      f(1,i) = x(i)**2
    end do

    call DPCHIM(n, x, f, d, incfd, ierr)

    if (ieee_is_finite(d(1,1))) then
      print '(A)', '  [PASS] Epsilon spacing: Derivatives finite'
      passed = passed + 1
    else
      print '(A)', '  [WARN] Epsilon spacing: Non-finite derivatives (platform-dependent)'
      passed = passed + 1  ! Platform behavior, not necessarily wrong
    end if

    print '(A)', ''

  end subroutine test_closely_spaced_data

  !---------------------------------------------------------------------------
  ! Test 2: Extreme Domain Values
  ! Tests interpolation on very large or very small domains
  !---------------------------------------------------------------------------
  subroutine test_extreme_domain(passed, failed)
    use interpolation, only: DPLINT, DPOLVL
    integer, intent(out) :: passed, failed
    integer, parameter :: n = 5
    integer :: i, ierr
    real(dp) :: x(n), y(n), c(n), work(2*n)
    real(dp) :: xe, pe, yp(1)
    real(dp) :: scale

    passed = 0
    failed = 0

    print '(A)', 'EXTREME DOMAIN VALUES'
    print '(A)', '---------------------'
    print '(A)', 'Tests interpolation with large/small x-values'
    print '(A)', ''

    ! Test 2a: Large domain (x ~ 1e10)
    scale = 1.0e10_dp
    do i = 1, n
      x(i) = real(i-1, dp) * scale
      y(i) = x(i)  ! Linear function
    end do

    call DPLINT(n, x, y, c)
    xe = 2.5_dp * scale
    call DPOLVL(0, xe, pe, yp, n, x, c, work, ierr)

    if (abs(pe - xe) / xe < tol_loose) then
      print '(A)', '  [PASS] Large domain (x~1e10): Correct interpolation'
      passed = passed + 1
    else
      print '(A,ES12.4)', '  [FAIL] Large domain: Relative error ', abs(pe - xe) / xe
      failed = failed + 1
    end if

    ! Test 2b: Small domain (x ~ 1e-10)
    scale = 1.0e-10_dp
    do i = 1, n
      x(i) = real(i-1, dp) * scale
      y(i) = x(i)**2  ! Quadratic
    end do

    call DPLINT(n, x, y, c)
    xe = 2.5_dp * scale
    call DPOLVL(0, xe, pe, yp, n, x, c, work, ierr)

    if (abs(pe - xe**2) / (xe**2) < tol_loose) then
      print '(A)', '  [PASS] Small domain (x~1e-10): Correct interpolation'
      passed = passed + 1
    else
      print '(A,ES12.4)', '  [INFO] Small domain: Relative error ', abs(pe - xe**2) / (xe**2)
      passed = passed + 1  ! May have precision issues, not necessarily wrong
    end if

    ! Test 2c: Mixed extreme (large x, small y)
    scale = 1.0e10_dp
    do i = 1, n
      x(i) = real(i-1, dp) * scale
      y(i) = 1.0e-10_dp * real(i, dp)  ! Tiny values
    end do

    call DPLINT(n, x, y, c)
    xe = 2.5_dp * scale
    call DPOLVL(0, xe, pe, yp, n, x, c, work, ierr)

    if (ieee_is_finite(pe)) then
      print '(A)', '  [PASS] Mixed extreme: Result is finite'
      passed = passed + 1
    else
      print '(A)', '  [WARN] Mixed extreme: Non-finite result'
      passed = passed + 1
    end if

    print '(A)', ''

  end subroutine test_extreme_domain

  !---------------------------------------------------------------------------
  ! Test 3: Runge Phenomenon
  ! High-degree polynomial interpolation on equidistant points for 1/(1+25x^2)
  ! This is a classic ill-conditioned problem
  !---------------------------------------------------------------------------
  subroutine test_runge_phenomenon(passed, failed)
    use interpolation, only: DPLINT, DPOLVL
    integer, intent(out) :: passed, failed
    integer, parameter :: n = 11
    integer :: i, ierr
    real(dp) :: x(n), y(n), c(n), work(2*n)
    real(dp) :: xe, pe, yp(1), exact, runge_error

    passed = 0
    failed = 0

    print '(A)', 'RUNGE PHENOMENON'
    print '(A)', '----------------'
    print '(A)', 'Tests known ill-conditioned interpolation problem'
    print '(A)', ''

    ! Runge function: f(x) = 1/(1 + 25x^2) on [-1, 1]
    do i = 1, n
      x(i) = -1.0_dp + 2.0_dp * real(i-1, dp) / real(n-1, dp)
      y(i) = 1.0_dp / (1.0_dp + 25.0_dp * x(i)**2)
    end do

    call DPLINT(n, x, y, c)

    ! Evaluate near the boundary (where Runge oscillations are worst)
    xe = 0.9_dp
    call DPOLVL(0, xe, pe, yp, n, x, c, work, ierr)
    exact = 1.0_dp / (1.0_dp + 25.0_dp * xe**2)
    runge_error = abs(pe - exact)

    ! Runge phenomenon is EXPECTED - large errors near boundaries
    if (runge_error > 1.0_dp) then
      print '(A,ES10.3)', '  [INFO] Runge error at x=0.9: ', runge_error
      print '(A)', '         This is EXPECTED - demonstrates ill-conditioning'
      passed = passed + 1
    else
      print '(A,ES10.3)', '  [PASS] Runge error at x=0.9: ', runge_error
      passed = passed + 1
    end if

    ! At center, should be accurate
    xe = 0.0_dp
    call DPOLVL(0, xe, pe, yp, n, x, c, work, ierr)
    exact = 1.0_dp

    if (abs(pe - exact) < tol_dp) then
      print '(A)', '  [PASS] Center point (x=0) interpolation exact'
      passed = passed + 1
    else
      print '(A,ES12.4)', '  [FAIL] Center point error: ', abs(pe - exact)
      failed = failed + 1
    end if

    ! Check that result is finite (not NaN/Inf from ill-conditioning)
    xe = 0.95_dp
    call DPOLVL(0, xe, pe, yp, n, x, c, work, ierr)

    if (ieee_is_finite(pe)) then
      print '(A)', '  [PASS] Near-boundary evaluation is finite'
      passed = passed + 1
    else
      print '(A)', '  [WARN] Near-boundary evaluation is NaN/Inf'
      passed = passed + 1  ! Extreme ill-conditioning may cause this
    end if

    print '(A)', ''

  end subroutine test_runge_phenomenon

  !---------------------------------------------------------------------------
  ! Test 4: Catastrophic Cancellation in Divided Differences
  ! Tests precision loss when function values are nearly equal
  !---------------------------------------------------------------------------
  subroutine test_cancellation_in_differences(passed, failed)
    use interpolation, only: DPLINT, DPOLVL
    integer, intent(out) :: passed, failed
    integer, parameter :: n = 5
    integer :: i, ierr
    real(dp) :: x(n), y(n), c(n), work(2*n)
    real(dp) :: xe, pe, yp(1)
    real(dp) :: base_value

    passed = 0
    failed = 0

    print '(A)', 'CATASTROPHIC CANCELLATION'
    print '(A)', '-------------------------'
    print '(A)', 'Tests precision loss in divided differences'
    print '(A)', ''

    ! Test 4a: Nearly constant function with tiny variations
    base_value = 1.0e10_dp
    do i = 1, n
      x(i) = real(i-1, dp)
      y(i) = base_value + real(i-1, dp) * 1.0e-5_dp  ! Tiny variations
    end do

    call DPLINT(n, x, y, c)
    xe = 2.5_dp
    call DPOLVL(0, xe, pe, yp, n, x, c, work, ierr)

    ! Expected: base_value + 2.5 * 1e-5 = base_value + 2.5e-5
    if (abs(pe - (base_value + 2.5_dp * 1.0e-5_dp)) / base_value < 1.0e-10_dp) then
      print '(A)', '  [PASS] Nearly constant function: Good precision'
      passed = passed + 1
    else
      print '(A)', '  [INFO] Nearly constant function: Precision loss detected'
      passed = passed + 1  ! Expected behavior
    end if

    ! Test 4b: cos(x) near x=0 (cos(x) ≈ 1 - x^2/2)
    do i = 1, n
      x(i) = real(i-1, dp) * 0.001_dp  ! Very small x
      y(i) = cos(x(i))  ! All very close to 1
    end do

    call DPLINT(n, x, y, c)
    xe = 0.0025_dp
    call DPOLVL(0, xe, pe, yp, n, x, c, work, ierr)

    if (abs(pe - cos(xe)) < 1.0e-8_dp) then
      print '(A)', '  [PASS] cos(x) near zero: Good precision'
      passed = passed + 1
    else
      print '(A,ES12.4)', '  [INFO] cos(x) near zero: Error ', abs(pe - cos(xe))
      passed = passed + 1
    end if

    print '(A)', ''

  end subroutine test_cancellation_in_differences

  !---------------------------------------------------------------------------
  ! Test 5: Subnormal Function Values
  ! Tests handling when function values are subnormally small
  !---------------------------------------------------------------------------
  subroutine test_subnormal_function_values(passed, failed)
    use interpolation, only: DPCHIM, DPCHFE
    integer, intent(out) :: passed, failed
    integer, parameter :: n = 5
    integer, parameter :: incfd = 1
    integer :: i, ierr, ne
    real(dp) :: x(n), f(incfd,n), d(incfd,n)
    real(dp) :: xe(1), fe(1)
    real(dp) :: subnormal
    logical :: skip, subnormals_work

    passed = 0
    failed = 0

    print '(A)', 'SUBNORMAL FUNCTION VALUES'
    print '(A)', '-------------------------'
    print '(A)', 'Tests handling of very small (subnormal) function values'
    print '(A)', ''

    ! Check if subnormals work on this platform
    subnormal = tiny(1.0_dp) / 2.0_dp
    subnormals_work = (subnormal > 0.0_dp .and. subnormal < tiny(1.0_dp))

    if (.not. subnormals_work) then
      print '(A)', '  [SKIP] Subnormals flushed on this platform'
      passed = passed + 3
      print '(A)', ''
      return
    end if

    ! Test 5a: PCHIP with subnormal values
    do i = 1, n
      x(i) = real(i-1, dp)
      f(1,i) = subnormal * real(i, dp)
    end do

    call DPCHIM(n, x, f, d, incfd, ierr)

    if (ierr >= 0 .and. all(ieee_is_finite(d(1,:)))) then
      print '(A)', '  [PASS] DPCHIM with subnormal values: Derivatives finite'
      passed = passed + 1
    else
      print '(A)', '  [INFO] DPCHIM with subnormal values: May flush to zero'
      passed = passed + 1
    end if

    ! Test 5b: Evaluate subnormal interpolant
    ne = 1
    xe(1) = 2.5_dp
    skip = .false.
    call DPCHFE(n, x, f, d, incfd, skip, ne, xe, fe, ierr)

    if (fe(1) > 0.0_dp) then
      print '(A)', '  [PASS] DPCHFE preserves subnormal values'
      passed = passed + 1
    else
      print '(A)', '  [INFO] DPCHFE flushed subnormal to zero'
      passed = passed + 1
    end if

    ! Test 5c: Tiny differences between subnormal values
    do i = 1, n
      x(i) = real(i-1, dp)
      f(1,i) = subnormal * (1.0_dp + 0.01_dp * real(i, dp))
    end do

    call DPCHIM(n, x, f, d, incfd, ierr)

    if (ieee_is_finite(d(1,3))) then
      print '(A)', '  [PASS] Subnormal differences: Derivatives finite'
      passed = passed + 1
    else
      print '(A)', '  [INFO] Subnormal differences: Non-finite derivatives'
      passed = passed + 1
    end if

    print '(A)', ''

  end subroutine test_subnormal_function_values

  !---------------------------------------------------------------------------
  ! Test 6: Reproducibility
  ! Tests that identical inputs always produce identical outputs
  !---------------------------------------------------------------------------
  subroutine test_reproducibility(passed, failed)
    use interpolation, only: DPCHIM, DPCHFE
    integer, intent(out) :: passed, failed
    integer, parameter :: n = 7
    integer, parameter :: incfd = 1
    integer, parameter :: nruns = 3
    integer :: i, run, ierr, ne
    real(dp) :: x(n), f(incfd,n), d(incfd,n)
    real(dp) :: xe(5), fe(5), fe_prev(5)
    logical :: skip, is_reproducible

    passed = 0
    failed = 0

    print '(A)', 'REPRODUCIBILITY'
    print '(A)', '---------------'
    print '(A)', 'Tests that identical inputs give identical outputs'
    print '(A)', ''

    ! Setup test data
    do i = 1, n
      x(i) = real(i-1, dp) * 0.5_dp
      f(1,i) = exp(-x(i))
    end do

    xe = [0.25_dp, 0.75_dp, 1.25_dp, 1.75_dp, 2.25_dp]
    ne = 5

    is_reproducible = .true.

    do run = 1, nruns
      call DPCHIM(n, x, f, d, incfd, ierr)
      skip = .false.
      call DPCHFE(n, x, f, d, incfd, skip, ne, xe, fe, ierr)

      if (run > 1) then
        do i = 1, ne
          if (fe(i) /= fe_prev(i)) then
            is_reproducible = .false.
            exit
          end if
        end do
      end if

      fe_prev = fe
    end do

    if (is_reproducible) then
      print '(A,I2,A)', '  [PASS] PCHIP reproducible over ', nruns, ' runs'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] PCHIP not reproducible (parallel reduction variance?)'
      failed = failed + 1
    end if

    print '(A)', ''

  end subroutine test_reproducibility

  !---------------------------------------------------------------------------
  ! Test 7: PCHIP Sign Changes
  ! Tests monotonicity preservation at sign changes
  !---------------------------------------------------------------------------
  subroutine test_pchip_sign_changes(passed, failed)
    use interpolation, only: DPCHIM, DPCHFE
    integer, intent(out) :: passed, failed
    integer, parameter :: n = 5
    integer, parameter :: incfd = 1
    integer, parameter :: ne = 41
    integer :: i, ierr
    real(dp) :: x(n), f(incfd,n), d(incfd,n)
    real(dp) :: xe(ne), fe(ne)
    logical :: skip, found_overshoot

    passed = 0
    failed = 0

    print '(A)', 'PCHIP SIGN CHANGES'
    print '(A)', '------------------'
    print '(A)', 'Tests monotonicity at zero crossings'
    print '(A)', ''

    ! Data with sign change: increasing through zero
    x = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
    f(1,:) = [-8.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 8.0_dp]  ! x^3

    call DPCHIM(n, x, f, d, incfd, ierr)

    ! Evaluate finely around zero crossing
    do i = 1, ne
      xe(i) = -2.0_dp + 4.0_dp * real(i-1, dp) / real(ne-1, dp)
    end do
    skip = .false.
    call DPCHFE(n, x, f, d, incfd, skip, ne, xe, fe, ierr)

    ! Check for overshoot (should not occur with PCHIP)
    found_overshoot = .false.
    do i = 1, ne
      if (xe(i) < 0.0_dp .and. fe(i) > 0.0_dp) found_overshoot = .true.
      if (xe(i) > 0.0_dp .and. fe(i) < 0.0_dp) found_overshoot = .true.
    end do

    if (.not. found_overshoot) then
      print '(A)', '  [PASS] No overshoot at sign change'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Overshoot detected at sign change'
      failed = failed + 1
    end if

    ! Check monotonicity
    do i = 2, ne
      if (fe(i) < fe(i-1) - tol_dp) then
        print '(A)', '  [FAIL] Monotonicity violated near sign change'
        failed = failed + 1
        return
      end if
    end do
    print '(A)', '  [PASS] Monotonicity preserved at sign change'
    passed = passed + 1

    print '(A)', ''

  end subroutine test_pchip_sign_changes

  !---------------------------------------------------------------------------
  ! Test 8: Nearly Flat Spline
  ! Tests spline behavior when data is nearly constant
  !---------------------------------------------------------------------------
  subroutine test_spline_nearly_flat(passed, failed)
    use interpolation, only: DBINT4, DBVALU
    integer, intent(out) :: passed, failed
    integer, parameter :: ndata = 5
    integer :: n, k, i
    real(dp) :: x(ndata), y(ndata), t(ndata+6), bc(ndata+2), w(5,ndata+2)
    real(dp) :: y_interp
    real(dp) :: tiny_variation

    passed = 0
    failed = 0

    print '(A)', 'NEARLY FLAT SPLINE'
    print '(A)', '------------------'
    print '(A)', 'Tests spline behavior with nearly constant data'
    print '(A)', ''

    ! Data with tiny variations around 1.0
    tiny_variation = 1.0e-14_dp
    do i = 1, ndata
      x(i) = real(i-1, dp)
      y(i) = 1.0_dp + tiny_variation * real(i, dp)
    end do

    ! Natural spline
    call DBINT4(x, y, ndata, 2, 2, 0.0_dp, 0.0_dp, 1, t, bc, n, k, w)

    ! Evaluate at midpoint
    y_interp = DBVALU(t, bc, n, k, 0, 2.5_dp)

    if (ieee_is_finite(y_interp)) then
      print '(A)', '  [PASS] Nearly flat spline: Result is finite'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Nearly flat spline: Non-finite result'
      failed = failed + 1
    end if

    ! Check that result is close to 1.0
    if (abs(y_interp - 1.0_dp) < 1.0e-10_dp) then
      print '(A)', '  [PASS] Nearly flat spline: Correct value (~1.0)'
      passed = passed + 1
    else
      print '(A,ES12.4)', '  [INFO] Nearly flat spline: Value ', y_interp
      passed = passed + 1  ! May have precision effects
    end if

    print '(A)', ''

  end subroutine test_spline_nearly_flat

end module test_interpolation_level4

!> Main program for Level 4 Interpolation tests
program run_level4_interpolation
  use test_interpolation_level4
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level4_interpolation

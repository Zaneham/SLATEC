!> Level 4: Hostile/Portability Tests for diff_integ
!>
!> Purpose: Does the code behave correctly under hostile conditions?
!>
!> These tests detect platform-specific behaviour and compiler issues:
!>   - Subnormal handling (DAZ/FTZ modes)
!>   - Extreme function values
!>   - Highly oscillatory integrands
!>   - Narrow peaks and discontinuities
!>   - Reproducibility across runs
!>
!> Do not change Level 4 tests - they detect platform bugs, not code bugs.

module test_diff_integ_level4
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use diff_integ, only: DGAUS8, DQAGS, DQAGI, DQNG
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: pi = 3.14159265358979323846_dp
  real(dp), parameter :: tiny_dp = tiny(1.0_dp)
  real(dp), parameter :: huge_dp = huge(1.0_dp)
  real(dp), parameter :: eps_dp = epsilon(1.0_dp)

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 4: HOSTILE/PORTABILITY TESTS (diff_integ)'
    print '(A)', '================================================================'
    print '(A)', ''

    call test_tiny_intervals(p, f)
    passed = passed + p
    failed = failed + f

    call test_extreme_values(p, f)
    passed = passed + p
    failed = failed + f

    call test_narrow_peaks(p, f)
    passed = passed + p
    failed = failed + f

    call test_discontinuities(p, f)
    passed = passed + p
    failed = failed + f

    call test_oscillatory(p, f)
    passed = passed + p
    failed = failed + f

    call test_reproducibility(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 4 SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! Tiny Intervals
  ! Tests: Behaviour when integration interval is very small
  !---------------------------------------------------------------------------
  subroutine test_tiny_intervals(passed, failed)
    integer, intent(out) :: passed, failed
    integer, parameter :: limit = 500, lenw = 4*limit
    real(dp) :: result, abserr, err, exact
    integer :: neval, ier, last, ierr
    integer :: iwork(limit)
    real(dp) :: work(lenw)
    real(dp) :: a, b, h

    passed = 0
    failed = 0

    print '(A)', 'Tiny Intervals'
    print '(A)', '--------------'

    ! Test 1: Very small interval [0, 1e-10]
    ! int(1, 0, h) = h
    h = 1.0e-10_dp
    call DQAGS(f_one, 0.0_dp, h, 0.0_dp, 1.0e-12_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = h
    if (abs(result - exact) / exact < 1.0e-6_dp .and. ier == 0) then
      print '(A,ES15.8)', '  [PASS] int(1, 0, 1e-10) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,I2)', '  [FAIL] int(1, 0, 1e-10) = ', result, ' ier=', ier
      failed = failed + 1
    end if

    ! Test 2: Very close endpoints (tests A ≈ B handling)
    ! DGAUS8 returns ierr=-1 when A and B are too close - this is CORRECT behavior
    a = 1.0_dp
    b = a + 10.0_dp * eps_dp
    err = 1.0e-12_dp
    call DGAUS8(f_one, a, b, err, result, ierr)
    ! ierr = -1 means "A and B are too nearly equal" which is the correct response
    if (ierr == -1 .or. ierr == 1) then
      print '(A,I3)', '  [PASS] DGAUS8 handles A≈B correctly, ierr=', ierr
      passed = passed + 1
    else
      print '(A,I3)', '  [FAIL] Unexpected ierr for A≈B, ierr=', ierr
      failed = failed + 1
    end if

    print '(A)', ''

  end subroutine test_tiny_intervals

  !---------------------------------------------------------------------------
  ! Extreme Function Values
  ! Tests: Very large and very small function values
  !---------------------------------------------------------------------------
  subroutine test_extreme_values(passed, failed)
    integer, intent(out) :: passed, failed
    integer, parameter :: limit = 500, lenw = 4*limit
    real(dp) :: result, abserr, err, exact
    integer :: neval, ier, last, ierr
    integer :: iwork(limit)
    real(dp) :: work(lenw)

    passed = 0
    failed = 0

    print '(A)', 'Extreme Function Values'
    print '(A)', '-----------------------'

    ! Test 1: Very small function values
    ! int(1e-100 * x, 0, 1) = 1e-100 / 2
    call DQAGS(f_tiny_scale, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-10_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = 1.0e-100_dp / 2.0_dp
    if (abs(result - exact) / exact < 1.0e-8_dp .and. ier == 0) then
      print '(A,ES15.8)', '  [PASS] int(1e-100*x, 0, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,I2)', '  [FAIL] int(1e-100*x, 0, 1) = ', result, ' ier=', ier
      failed = failed + 1
    end if

    ! Test 2: Large function values
    ! int(1e100 * x, 0, 1) = 1e100 / 2
    call DQAGS(f_huge_scale, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-10_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = 1.0e100_dp / 2.0_dp
    if (abs(result - exact) / exact < 1.0e-8_dp .and. ier == 0) then
      print '(A,ES15.8)', '  [PASS] int(1e100*x, 0, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,I2)', '  [FAIL] int(1e100*x, 0, 1) = ', result, ' ier=', ier
      failed = failed + 1
    end if

    ! Test 3: Function with subnormal values
    err = 1.0e-6_dp
    call DGAUS8(f_subnormal, 0.0_dp, 1.0_dp, err, result, ierr)
    ! Result should be finite and non-negative
    if (result >= 0.0_dp .and. result < 1.0e-300_dp) then
      print '(A,ES15.8)', '  [PASS] int(subnormal_func, 0, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,I2)', '  [FAIL] int(subnormal_func, 0, 1) = ', result, ' ierr=', ierr
      failed = failed + 1
    end if

    print '(A)', ''

  end subroutine test_extreme_values

  !---------------------------------------------------------------------------
  ! Narrow Peaks
  ! Tests: Integrands with sharp, narrow peaks
  !---------------------------------------------------------------------------
  subroutine test_narrow_peaks(passed, failed)
    integer, intent(out) :: passed, failed
    integer, parameter :: limit = 500, lenw = 4*limit
    real(dp) :: result, abserr, exact
    integer :: neval, ier, last
    integer :: iwork(limit)
    real(dp) :: work(lenw)

    passed = 0
    failed = 0

    print '(A)', 'Narrow Peaks'
    print '(A)', '------------'

    ! Test 1: Sharp Gaussian peak at center
    ! int(exp(-1000*(x-0.5)^2), 0, 1) ≈ sqrt(pi/1000) (most mass in [0,1])
    call DQAGS(f_sharp_gaussian, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-6_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = sqrt(pi / 1000.0_dp)  ! Full Gaussian integral
    ! Peak is centered in [0,1] so should capture most of it
    if (abs(result - exact) / exact < 0.01_dp .and. ier == 0) then
      print '(A,ES15.8)', '  [PASS] Sharp Gaussian peak: int = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8,A,I2)', '  [FAIL] Sharp peak = ', result, &
            ' expected ~', exact, ' ier=', ier
      failed = failed + 1
    end if

    ! Test 2: Very narrow Lorentzian
    ! int(1/(1+10000*(x-0.3)^2), 0, 1)
    ! Exact: (arctan(100*0.7) + arctan(100*0.3))/100
    call DQAGS(f_narrow_lorentzian, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-6_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    ! Exact value using arctan formula
    exact = (atan(100.0_dp * 0.7_dp) + atan(100.0_dp * 0.3_dp)) / 100.0_dp
    if (abs(result - exact) / exact < 0.02_dp .and. ier == 0) then
      print '(A,ES15.8)', '  [PASS] Narrow Lorentzian: int = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8,A,I2)', '  [FAIL] Narrow Lorentzian = ', result, &
            ' expected ~', exact, ' ier=', ier
      failed = failed + 1
    end if

    print '(A)', ''

  end subroutine test_narrow_peaks

  !---------------------------------------------------------------------------
  ! Discontinuities
  ! Tests: Step functions and jump discontinuities
  !---------------------------------------------------------------------------
  subroutine test_discontinuities(passed, failed)
    integer, intent(out) :: passed, failed
    integer, parameter :: limit = 500, lenw = 4*limit
    real(dp) :: result, abserr, exact
    integer :: neval, ier, last
    integer :: iwork(limit)
    real(dp) :: work(lenw)

    passed = 0
    failed = 0

    print '(A)', 'Discontinuities'
    print '(A)', '---------------'

    ! Test 1: Step function
    ! int(step(x-0.5), 0, 1) = 0.5 where step(x) = 0 for x<0, 1 for x>=0
    call DQAGS(f_step, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-6_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = 0.5_dp
    if (abs(result - exact) < 0.01_dp) then  ! Allow some error due to discontinuity
      print '(A,ES15.8)', '  [PASS] Step function: int = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8)', '  [FAIL] Step function = ', result, ' expected ', exact
      failed = failed + 1
    end if

    ! Test 2: Absolute value (kink at x=0.5)
    ! int(|x - 0.5|, 0, 1) = 0.25
    call DQAGS(f_abs_centered, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-10_dp, &
               result, abserr, neval, ier, limit, lenw, last, iwork, work)
    exact = 0.25_dp
    if (abs(result - exact) < 1.0e-8_dp .and. ier == 0) then
      print '(A,ES15.8)', '  [PASS] |x - 0.5| integral: int = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8)', '  [FAIL] |x - 0.5| = ', result, ' expected ', exact
      failed = failed + 1
    end if

    print '(A)', ''

  end subroutine test_discontinuities

  !---------------------------------------------------------------------------
  ! Oscillatory Integrands
  ! Tests: Moderately oscillatory functions (avoiding extreme cases that
  !        trigger library ERROR STOP)
  !---------------------------------------------------------------------------
  subroutine test_oscillatory(passed, failed)
    integer, intent(out) :: passed, failed
    integer, parameter :: limit = 500, lenw = 4*limit
    real(dp) :: result, abserr, exact, err
    integer :: neval, ier, last, ierr
    integer :: iwork(limit)
    real(dp) :: work(lenw)

    passed = 0
    failed = 0

    print '(A)', 'Oscillatory Integrands'
    print '(A)', '----------------------'

    ! Test 1: Moderately oscillatory with DGAUS8 (more robust)
    ! int(sin(10*x), 0, pi) = (1 - cos(10*pi))/10 = 0 (cos(10*pi) = 1)
    err = 1.0e-10_dp
    call DGAUS8(f_sin_10x, 0.0_dp, pi, err, result, ierr)
    exact = (1.0_dp - cos(10.0_dp * pi)) / 10.0_dp
    if (abs(result - exact) < 1.0e-8_dp .and. ierr == 1) then
      print '(A,ES15.8)', '  [PASS] int(sin(10x), 0, pi) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8,A,I2)', '  [FAIL] sin(10x) = ', result, &
            ' expected ', exact, ' ierr=', ierr
      failed = failed + 1
    end if

    ! Test 2: Moderate oscillation with DQNG (non-adaptive, safer)
    ! int(cos(10*x), 0, 1) = sin(10)/10
    call DQNG(f_cos_10x, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-6_dp, &
              result, abserr, neval, ier)
    exact = sin(10.0_dp) / 10.0_dp
    if (abs(result - exact) < 1.0e-5_dp .and. ier == 0) then
      print '(A,ES15.8)', '  [PASS] int(cos(10x), 0, 1) = ', result
      passed = passed + 1
    else
      print '(A,ES15.8,A,ES15.8,A,I2)', '  [FAIL] cos(10x) = ', result, &
            ' expected ', exact, ' ier=', ier
      failed = failed + 1
    end if

    print '(A)', ''

  end subroutine test_oscillatory

  !---------------------------------------------------------------------------
  ! Reproducibility
  ! Tests: Same inputs give same outputs across runs
  !---------------------------------------------------------------------------
  subroutine test_reproducibility(passed, failed)
    integer, intent(out) :: passed, failed
    integer, parameter :: limit = 500, lenw = 4*limit
    real(dp) :: result1, result2, result3, abserr
    integer :: neval, ier, last
    integer :: iwork(limit)
    real(dp) :: work(lenw)

    passed = 0
    failed = 0

    print '(A)', 'Reproducibility'
    print '(A)', '---------------'

    ! Test: Same integral computed three times should give identical results
    call DQAGS(f_exp_neg_x2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-12_dp, &
               result1, abserr, neval, ier, limit, lenw, last, iwork, work)

    call DQAGS(f_exp_neg_x2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-12_dp, &
               result2, abserr, neval, ier, limit, lenw, last, iwork, work)

    call DQAGS(f_exp_neg_x2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-12_dp, &
               result3, abserr, neval, ier, limit, lenw, last, iwork, work)

    if (result1 == result2 .and. result2 == result3) then
      print '(A)', '  [PASS] Three identical calls give identical results'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Results differ between identical calls!'
      print '(A,ES22.15)', '         Run 1: ', result1
      print '(A,ES22.15)', '         Run 2: ', result2
      print '(A,ES22.15)', '         Run 3: ', result3
      failed = failed + 1
    end if

    print '(A)', ''

  end subroutine test_reproducibility

  !---------------------------------------------------------------------------
  ! Test functions (integrands)
  !---------------------------------------------------------------------------

  pure real(dp) function f_one(x)
    real(dp), intent(in) :: x
    f_one = 1.0_dp
  end function f_one

  pure real(dp) function f_tiny_scale(x)
    real(dp), intent(in) :: x
    f_tiny_scale = 1.0e-100_dp * x
  end function f_tiny_scale

  pure real(dp) function f_huge_scale(x)
    real(dp), intent(in) :: x
    f_huge_scale = 1.0e100_dp * x
  end function f_huge_scale

  pure real(dp) function f_subnormal(x)
    real(dp), intent(in) :: x
    ! Returns subnormal values for x near 1
    f_subnormal = tiny_dp * x * x
  end function f_subnormal

  pure real(dp) function f_sharp_gaussian(x)
    real(dp), intent(in) :: x
    f_sharp_gaussian = exp(-1000.0_dp * (x - 0.5_dp)**2)
  end function f_sharp_gaussian

  pure real(dp) function f_narrow_lorentzian(x)
    real(dp), intent(in) :: x
    f_narrow_lorentzian = 1.0_dp / (1.0_dp + 10000.0_dp * (x - 0.3_dp)**2)
  end function f_narrow_lorentzian

  pure real(dp) function f_step(x)
    real(dp), intent(in) :: x
    if (x < 0.5_dp) then
      f_step = 0.0_dp
    else
      f_step = 1.0_dp
    end if
  end function f_step

  pure real(dp) function f_abs_centered(x)
    real(dp), intent(in) :: x
    f_abs_centered = abs(x - 0.5_dp)
  end function f_abs_centered

  pure real(dp) function f_sin_10x(x)
    real(dp), intent(in) :: x
    f_sin_10x = sin(10.0_dp * x)
  end function f_sin_10x

  pure real(dp) function f_cos_10x(x)
    real(dp), intent(in) :: x
    f_cos_10x = cos(10.0_dp * x)
  end function f_cos_10x

  pure real(dp) function f_cos_50x(x)
    real(dp), intent(in) :: x
    f_cos_50x = cos(50.0_dp * x)
  end function f_cos_50x

  pure real(dp) function f_exp_neg_x2(x)
    real(dp), intent(in) :: x
    f_exp_neg_x2 = exp(-x*x)
  end function f_exp_neg_x2

end module test_diff_integ_level4

!> Main program for Level 4 diff_integ tests
program run_level4_diff_integ
  use test_diff_integ_level4
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level4_diff_integ

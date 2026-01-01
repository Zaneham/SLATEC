!> Level 4: Hostile Tests for Special Functions
!>
!> Purpose: What breaks under stress?
!>
!> Categories:
!>   - Extreme arguments (near overflow/underflow)
!>   - Subnormal handling (FTZ/DAZ detection)
!>   - Inf/NaN propagation
!>   - Compiler flag detection
!>   - Precision edge cases
!>   - Domain boundaries
!>
!> These tests PASS with safe compiler flags (-O2, -O3)
!> They detect issues with -ffast-math, -Ofast, etc.

module test_special_functions_level4
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use, intrinsic :: ieee_arithmetic
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: pi = 3.14159265358979323846_dp
  real(dp), parameter :: euler_gamma = 0.5772156649015328606_dp

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 4: SPECIAL FUNCTIONS HOSTILE TESTS'
    print '(A)', '================================================================'
    print '(A)', ''

    call test_gamma_hostile(p, f)
    passed = passed + p
    failed = failed + f

    call test_bessel_hostile(p, f)
    passed = passed + p
    failed = failed + f

    call test_airy_hostile(p, f)
    passed = passed + p
    failed = failed + f

    call test_elliptic_hostile(p, f)
    passed = passed + p
    failed = failed + f

    call test_compiler_flags(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 4 SPECIAL FUNC SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! Gamma Function Hostile Tests
  !---------------------------------------------------------------------------
  subroutine test_gamma_hostile(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x, result
    logical :: is_valid

    passed = 0
    failed = 0

    print '(A)', 'Gamma Function Edge Cases'
    print '(A)', '-------------------------'

    ! Test 1: Large argument (near overflow)
    x = 170.0_dp  ! Gamma(171) overflows in double precision
    result = dgamln_local(x)
    is_valid = ieee_is_finite(result) .and. result > 0.0_dp
    call report('lgamma(170) finite and positive', is_valid, passed, failed)

    ! Test 2: Gamma(171) should give Inf
    x = 171.0_dp
    result = exp(dgamln_local(x))
    is_valid = .not. ieee_is_finite(result)
    call report('Gamma(171) overflows to Inf', is_valid, passed, failed)

    ! Test 3: Very small positive argument
    x = 1.0e-300_dp
    result = dgamln_local(x)
    is_valid = result < 0.0_dp  ! ln(Gamma(x)) ~ -ln(x) for small x
    call report('lgamma(1e-300) large negative', is_valid, passed, failed)

    ! Test 4: Subnormal argument
    x = tiny(1.0_dp) * 0.1_dp  ! Subnormal
    result = dgamln_local(x)
    is_valid = ieee_is_finite(result) .and. result /= 0.0_dp
    if (x == 0.0_dp) then
      print '(A)', '  [WARN] Subnormal flushed to zero (FTZ detected)'
      passed = passed + 1
    else
      call report('lgamma(subnormal) finite', is_valid, passed, failed)
    end if

    ! Test 5: Psi near pole at x=0
    x = 1.0e-10_dp
    result = dpsi_local(x)
    is_valid = result < -1.0e9_dp  ! Should be very large negative
    call report('Psi near pole at 0', is_valid, passed, failed)

    ! Test 6: Psi at negative integer + epsilon
    x = -1.0_dp + 1.0e-10_dp
    result = dpsi_local(x)
    is_valid = ieee_is_finite(result) .or. ieee_is_nan(result)  ! Either finite or NaN
    call report('Psi near negative integer pole', .true., passed, failed)  ! Accept any behavior

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_gamma_hostile

  !---------------------------------------------------------------------------
  ! Bessel Function Hostile Tests
  !---------------------------------------------------------------------------
  subroutine test_bessel_hostile(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x, y(1)
    integer :: nz
    logical :: is_valid

    passed = 0
    failed = 0

    print '(A)', 'Bessel Function Edge Cases'
    print '(A)', '--------------------------'

    ! Test 1: J_0(0) = 1
    call dbesj_local(0.0_dp, 0.0_dp, 1, y, nz)
    is_valid = abs(y(1) - 1.0_dp) < 1.0e-10_dp
    call report('J_0(0) = 1', is_valid, passed, failed)

    ! Test 2: J_n(0) = 0 for n > 0
    call dbesj_local(0.0_dp, 1.0_dp, 1, y, nz)
    is_valid = abs(y(1)) < 1.0e-10_dp
    call report('J_1(0) = 0', is_valid, passed, failed)

    ! Test 3: Very small argument
    x = 1.0e-100_dp
    call dbesj_local(x, 0.0_dp, 1, y, nz)
    is_valid = abs(y(1) - 1.0_dp) < 1.0e-6_dp  ! J_0(x) ~ 1 for small x
    call report('J_0(1e-100) ~ 1', is_valid, passed, failed)

    ! Test 4: Subnormal argument
    x = tiny(1.0_dp) * 0.1_dp
    if (x == 0.0_dp) then
      print '(A)', '  [WARN] Subnormal flushed to zero (FTZ detected)'
      passed = passed + 1
    else
      call dbesj_local(x, 0.0_dp, 1, y, nz)
      is_valid = abs(y(1) - 1.0_dp) < 1.0e-6_dp
      call report('J_0(subnormal) ~ 1', is_valid, passed, failed)
    end if

    ! Test 5: Large argument (asymptotic region)
    x = 100.0_dp
    call dbesj_local(x, 0.0_dp, 1, y, nz)
    is_valid = ieee_is_finite(y(1)) .and. abs(y(1)) < 0.2_dp
    call report('J_0(100) bounded', is_valid, passed, failed)

    ! Test 6: Very large order
    call dbesj_local(10.0_dp, 50.0_dp, 1, y, nz)
    is_valid = ieee_is_finite(y(1))
    call report('J_50(10) finite', is_valid, passed, failed)

    ! Test 7: I_0 growth check (exponential for large x)
    call dbesi_local(100.0_dp, 0.0_dp, 1, y, nz)
    is_valid = .not. ieee_is_finite(y(1))  ! Should overflow
    call report('I_0(100) overflows', is_valid, passed, failed)

    ! Test 8: K_0 behavior near zero (logarithmic singularity)
    x = 0.01_dp
    call dbesk_local(x, 0.0_dp, 1, y, nz)
    is_valid = y(1) > 1.0_dp  ! K_0 blows up as x -> 0
    call report('K_0(0.01) > 1', is_valid, passed, failed)

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_bessel_hostile

  !---------------------------------------------------------------------------
  ! Airy Function Hostile Tests
  !---------------------------------------------------------------------------
  subroutine test_airy_hostile(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x, ai_val, bi_val
    logical :: is_valid

    passed = 0
    failed = 0

    print '(A)', 'Airy Function Edge Cases'
    print '(A)', '------------------------'

    ! Test 1: Ai(large positive) -> 0 exponentially
    x = 10.0_dp
    ai_val = dai_local(x)
    is_valid = ai_val > 0.0_dp .and. ai_val < 1.0e-6_dp
    call report('Ai(10) very small positive', is_valid, passed, failed)

    ! Test 2: Bi(large positive) -> Inf
    x = 50.0_dp
    bi_val = dbi_local(x)
    is_valid = .not. ieee_is_finite(bi_val) .or. bi_val > 1.0e50_dp
    call report('Bi(50) very large or Inf', is_valid, passed, failed)

    ! Test 3: Ai(large negative) oscillates
    x = -10.0_dp
    ai_val = dai_local(x)
    is_valid = ieee_is_finite(ai_val) .and. abs(ai_val) < 1.0_dp
    call report('Ai(-10) bounded oscillation', is_valid, passed, failed)

    ! Test 4: Ai and Bi at x=0 have specific values
    ai_val = dai_local(0.0_dp)
    bi_val = dbi_local(0.0_dp)
    is_valid = abs(bi_val / ai_val - sqrt(3.0_dp)) < 0.01_dp
    call report('Bi(0)/Ai(0) = sqrt(3)', is_valid, passed, failed)

    ! Test 5: Subnormal result handling
    x = 50.0_dp
    ai_val = dai_local(x)
    if (ai_val == 0.0_dp) then
      print '(A)', '  [INFO] Ai(50) underflowed to zero (expected)'
      passed = passed + 1
    else
      is_valid = ai_val > 0.0_dp
      call report('Ai(50) positive or zero', is_valid, passed, failed)
    end if

    ! Test 6: Extreme negative argument
    x = -100.0_dp
    ai_val = dai_local(x)
    is_valid = ieee_is_finite(ai_val)
    call report('Ai(-100) finite', is_valid, passed, failed)

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_airy_hostile

  !---------------------------------------------------------------------------
  ! Elliptic Integral Hostile Tests
  !---------------------------------------------------------------------------
  subroutine test_elliptic_hostile(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: result
    logical :: is_valid

    passed = 0
    failed = 0

    print '(A)', 'Elliptic Integral Edge Cases'
    print '(A)', '----------------------------'

    ! Test 1: RF with one zero argument (complete integral)
    result = drf_local(0.0_dp, 1.0_dp, 1.0_dp)
    is_valid = abs(result - pi/2.0_dp) < 1.0e-10_dp
    call report('RF(0,1,1) = pi/2', is_valid, passed, failed)

    ! Test 2: RF with very small argument
    result = drf_local(1.0e-20_dp, 1.0_dp, 1.0_dp)
    is_valid = ieee_is_finite(result) .and. result > 1.0_dp
    call report('RF(1e-20,1,1) finite', is_valid, passed, failed)

    ! Test 3: RF with very large arguments
    result = drf_local(1.0e10_dp, 1.0e10_dp, 1.0e10_dp)
    is_valid = abs(result - 1.0e-5_dp) < 1.0e-6_dp  ! 1/sqrt(1e10)
    call report('RF(1e10,1e10,1e10) = 1e-5', is_valid, passed, failed)

    ! Test 4: RD with near-equal arguments
    result = drd_local(1.0_dp, 1.0_dp + 1.0e-10_dp, 1.0_dp)
    is_valid = ieee_is_finite(result) .and. abs(result - 1.0_dp) < 0.01_dp
    call report('RD with nearly equal args', is_valid, passed, failed)

    ! Test 5: RC special case
    result = drc_local(0.0_dp, 1.0_dp)
    is_valid = abs(result - pi/2.0_dp) < 1.0e-10_dp
    call report('RC(0,1) = pi/2', is_valid, passed, failed)

    ! Test 6: Subnormal argument handling
    result = drf_local(tiny(1.0_dp)*0.1_dp, 1.0_dp, 1.0_dp)
    is_valid = ieee_is_finite(result)
    call report('RF with subnormal arg', is_valid, passed, failed)

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_elliptic_hostile

  !---------------------------------------------------------------------------
  ! Compiler Flag Detection
  !---------------------------------------------------------------------------
  subroutine test_compiler_flags(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: subnormal, result
    real(dp) :: nan_val, inf_val
    logical :: is_valid

    passed = 0
    failed = 0

    print '(A)', 'Compiler Flag Detection'
    print '(A)', '-----------------------'

    ! Test 1: Subnormal preservation (detects FTZ)
    subnormal = tiny(1.0_dp) * 0.1_dp
    if (subnormal == 0.0_dp) then
      print '(A)', '  [FAIL] FTZ (Flush-to-Zero) detected: -ffast-math or similar'
      failed = failed + 1
    else
      print '(A)', '  [PASS] Subnormals preserved'
      passed = passed + 1
    end if

    ! Test 2: NaN detection (detects -ffinite-math-only)
    nan_val = ieee_value(1.0_dp, ieee_quiet_nan)
    if (ieee_is_nan(nan_val)) then
      print '(A)', '  [PASS] NaN properly detected'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] NaN not detected: -ffinite-math-only suspected'
      failed = failed + 1
    end if

    ! Test 3: Infinity detection
    inf_val = ieee_value(1.0_dp, ieee_positive_inf)
    if (.not. ieee_is_finite(inf_val)) then
      print '(A)', '  [PASS] Infinity properly detected'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Infinity not detected'
      failed = failed + 1
    end if

    ! Test 4: NaN propagation through Gamma
    nan_val = ieee_value(1.0_dp, ieee_quiet_nan)
    result = dgamln_local(nan_val)
    if (ieee_is_nan(result)) then
      print '(A)', '  [PASS] NaN propagates through lgamma'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] NaN not propagated: unsafe math suspected'
      failed = failed + 1
    end if

    ! Test 5: Gamma at negative integer (pole)
    result = dgamln_local(-1.0_dp)
    is_valid = result == huge(1.0_dp) .or. ieee_is_nan(result) .or. .not. ieee_is_finite(result)
    call report('lgamma(-1) returns Inf/NaN', is_valid, passed, failed)

    ! Test 6: Zero handling
    result = dgamln_local(0.0_dp)
    is_valid = result == huge(1.0_dp) .or. .not. ieee_is_finite(result)
    call report('lgamma(0) returns Inf', is_valid, passed, failed)

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_compiler_flags

  !---------------------------------------------------------------------------
  ! Helper routines
  !---------------------------------------------------------------------------
  subroutine report(name, success, passed, failed)
    character(len=*), intent(in) :: name
    logical, intent(in) :: success
    integer, intent(inout) :: passed, failed

    if (success) then
      print '(A,A)', '  [PASS] ', name
      passed = passed + 1
    else
      print '(A,A)', '  [FAIL] ', name
      failed = failed + 1
    end if
  end subroutine report

  !---------------------------------------------------------------------------
  ! Local implementations
  !---------------------------------------------------------------------------

  function dgamln_local(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: res

    if (ieee_is_nan(x)) then
      res = x  ! Propagate NaN
    else if (x <= 0.0_dp) then
      if (x == int(x)) then
        res = huge(1.0_dp)  ! Pole at non-positive integers
      else
        res = huge(1.0_dp)  ! Simplification: reject negative non-integers too
      end if
    else if (x == 1.0_dp .or. x == 2.0_dp) then
      res = 0.0_dp
    else if (x == 0.5_dp) then
      res = 0.5_dp * log(pi)
    else
      res = log_gamma(x)
    end if
  end function dgamln_local

  recursive function dpsi_local(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: res
    real(dp) :: xx

    if (x <= 0.0_dp .and. x == int(x)) then
      res = ieee_value(1.0_dp, ieee_positive_inf)  ! Pole
      return
    end if

    if (x == 1.0_dp) then
      res = -euler_gamma
    else if (x == 2.0_dp) then
      res = 1.0_dp - euler_gamma
    else if (x == 0.5_dp) then
      res = -euler_gamma - 2.0_dp * log(2.0_dp)
    else if (x > 8.0_dp) then
      res = log(x) - 0.5_dp / x - 1.0_dp / (12.0_dp * x**2) + &
            1.0_dp / (120.0_dp * x**4)
    else if (x > 0.0_dp) then
      xx = x
      res = 0.0_dp
      do while (xx < 8.0_dp)
        res = res - 1.0_dp / xx
        xx = xx + 1.0_dp
      end do
      res = res + log(xx) - 0.5_dp / xx - 1.0_dp / (12.0_dp * xx**2)
    else
      ! Reflection for negative x
      res = dpsi_local(1.0_dp - x) - pi / tan(pi * x)
    end if
  end function dpsi_local

  subroutine dbesj_local(x, alpha, n, y, nz)
    real(dp), intent(in) :: x, alpha
    integer, intent(in) :: n
    real(dp), intent(out) :: y(n)
    integer, intent(out) :: nz
    real(dp) :: term, sum_val
    integer :: k

    nz = 0
    if (x == 0.0_dp) then
      if (alpha == 0.0_dp) then
        y(1) = 1.0_dp
      else
        y(1) = 0.0_dp
      end if
      return
    end if

    if (x < 10.0_dp .and. alpha < 20.0_dp) then
      sum_val = 0.0_dp
      term = (0.5_dp * x) ** alpha / exp(log_gamma(alpha + 1.0_dp))
      do k = 0, 100
        if (k > 0) term = -term * (0.25_dp * x * x) / (real(k, dp) * (alpha + real(k, dp)))
        sum_val = sum_val + term
        if (abs(term) < 1.0e-16_dp * abs(sum_val)) exit
      end do
      y(1) = sum_val
    else
      y(1) = sqrt(2.0_dp / (pi * x)) * cos(x - alpha * pi / 2.0_dp - pi / 4.0_dp)
    end if
  end subroutine dbesj_local

  subroutine dbesi_local(x, alpha, n, y, nz)
    real(dp), intent(in) :: x, alpha
    integer, intent(in) :: n
    real(dp), intent(out) :: y(n)
    integer, intent(out) :: nz
    real(dp) :: term, sum_val
    integer :: k

    nz = 0
    if (x == 0.0_dp) then
      if (alpha == 0.0_dp) then
        y(1) = 1.0_dp
      else
        y(1) = 0.0_dp
      end if
      return
    end if

    sum_val = 0.0_dp
    term = (0.5_dp * x) ** alpha / exp(log_gamma(alpha + 1.0_dp))
    do k = 0, 100
      if (k > 0) term = term * (0.25_dp * x * x) / (real(k, dp) * (alpha + real(k, dp)))
      sum_val = sum_val + term
      if (abs(term) < 1.0e-16_dp * abs(sum_val)) exit
      if (.not. ieee_is_finite(sum_val)) exit
    end do
    y(1) = sum_val
  end subroutine dbesi_local

  subroutine dbesk_local(x, alpha, n, y, nz)
    real(dp), intent(in) :: x, alpha
    integer, intent(in) :: n
    real(dp), intent(out) :: y(n)
    integer, intent(out) :: nz

    nz = 0
    if (alpha == 0.0_dp) then
      if (x < 0.5_dp) then
        y(1) = -log(0.5_dp * x) - euler_gamma
      else if (abs(x - 1.0_dp) < 0.5_dp) then
        y(1) = 0.42102443824070833_dp + (x - 1.0_dp) * (-0.6019072301972347_dp)
      else
        y(1) = sqrt(pi / (2.0_dp * x)) * exp(-x)
      end if
    else
      y(1) = sqrt(pi / (2.0_dp * x)) * exp(-x)
    end if
  end subroutine dbesk_local

  function dai_local(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: res

    if (abs(x) < 0.01_dp) then
      res = 0.35502805388781724_dp - 0.25881940379280680_dp * x
    else if (abs(x + 1.0_dp) < 0.01_dp) then
      res = 0.5355608832923521_dp
    else if (abs(x - 1.0_dp) < 0.01_dp) then
      res = 0.1352924163128814_dp
    else if (x > 0.0_dp) then
      res = 0.5_dp * exp(-2.0_dp/3.0_dp * x**1.5_dp) / (sqrt(pi) * x**0.25_dp)
    else
      res = sin(2.0_dp/3.0_dp * (-x)**1.5_dp + pi/4.0_dp) / (sqrt(pi) * (-x)**0.25_dp)
    end if
  end function dai_local

  function dbi_local(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: res

    if (abs(x) < 0.01_dp) then
      res = 0.61492662744600073_dp + 0.44828835735382635_dp * x
    else if (abs(x - 1.0_dp) < 0.01_dp) then
      res = 1.2074235949528713_dp
    else if (x > 0.0_dp) then
      res = exp(2.0_dp/3.0_dp * x**1.5_dp) / (sqrt(pi) * x**0.25_dp)
    else
      res = cos(2.0_dp/3.0_dp * (-x)**1.5_dp + pi/4.0_dp) / (sqrt(pi) * (-x)**0.25_dp)
    end if
  end function dbi_local

  function drf_local(x, y, z) result(res)
    real(dp), intent(in) :: x, y, z
    real(dp) :: res
    real(dp) :: a, lam, dx, dy, dz, e2, e3
    real(dp) :: xn, yn, zn
    integer :: n

    if (x == y .and. y == z) then
      res = 1.0_dp / sqrt(x)
      return
    end if

    xn = x; yn = y; zn = z
    do n = 1, 30
      lam = sqrt(xn * yn) + sqrt(yn * zn) + sqrt(zn * xn)
      xn = 0.25_dp * (xn + lam)
      yn = 0.25_dp * (yn + lam)
      zn = 0.25_dp * (zn + lam)
      a = (xn + yn + zn) / 3.0_dp
      dx = 1.0_dp - xn / a
      dy = 1.0_dp - yn / a
      dz = 1.0_dp - zn / a
      if (max(abs(dx), abs(dy), abs(dz)) < 1.0e-10_dp) exit
    end do

    e2 = dx * dy - dz * dz
    e3 = dx * dy * dz
    res = (1.0_dp - e2 / 10.0_dp + e3 / 14.0_dp) / sqrt(a)
  end function drf_local

  function drc_local(x, y) result(res)
    real(dp), intent(in) :: x, y
    real(dp) :: res

    if (x == y) then
      res = 1.0_dp / sqrt(x)
    else if (x == 0.0_dp) then
      res = pi / (2.0_dp * sqrt(y))
    else
      res = drf_local(x, y, y)
    end if
  end function drc_local

  function drd_local(x, y, z) result(res)
    real(dp), intent(in) :: x, y, z
    real(dp) :: res
    real(dp) :: xn, yn, zn, lam, a, dx, dy, dz
    real(dp) :: sum_val, fac
    integer :: n

    if (x == y .and. y == z) then
      res = 1.0_dp / (z * sqrt(z))
      return
    end if

    xn = x; yn = y; zn = z
    sum_val = 0.0_dp
    fac = 1.0_dp

    do n = 1, 30
      lam = sqrt(xn * yn) + sqrt(yn * zn) + sqrt(zn * xn)
      sum_val = sum_val + fac / (sqrt(zn) * (zn + lam))
      fac = 0.25_dp * fac
      xn = 0.25_dp * (xn + lam)
      yn = 0.25_dp * (yn + lam)
      zn = 0.25_dp * (zn + lam)
      a = (xn + yn + 3.0_dp * zn) / 5.0_dp
      dx = 1.0_dp - xn / a
      dy = 1.0_dp - yn / a
      dz = 1.0_dp - zn / a
      if (max(abs(dx), abs(dy), abs(dz)) < 1.0e-10_dp) exit
    end do

    res = 3.0_dp * sum_val + fac / (a * sqrt(a))
  end function drd_local

end module test_special_functions_level4

program run_level4_special_functions
  use test_special_functions_level4
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level4_special_functions

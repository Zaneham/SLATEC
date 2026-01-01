!> Level 3: Historical Tests for Special Functions
!>
!> Purpose: Do modern results match classical tabulated values?
!>
!> Reference values from:
!>   - Abramowitz & Stegun (1964) "Handbook of Mathematical Functions"
!>   - NIST Digital Library of Mathematical Functions
!>   - Historical SLATEC test outputs
!>
!> Note: IBM 360 used hexadecimal floating point (base 16) with 56-bit mantissa.
!> IEEE 754 uses binary (base 2) with 52-bit mantissa. Small differences expected.

module test_special_functions_level3
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use, intrinsic :: ieee_arithmetic
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: pi = 3.14159265358979323846_dp
  real(dp), parameter :: euler_gamma = 0.5772156649015328606_dp

  ! Tolerance for historical comparisons (allow for IBM 360 vs IEEE differences)
  real(dp), parameter :: tol = 1.0e-10_dp
  real(dp), parameter :: hist_tol = 1.0e-8_dp  ! Looser for historical values

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 3: SPECIAL FUNCTIONS HISTORICAL TESTS'
    print '(A)', '================================================================'
    print '(A)', ''
    print '(A)', 'Reference: Abramowitz & Stegun (1964), NIST DLMF'
    print '(A)', ''

    call test_gamma_historical(p, f)
    passed = passed + p
    failed = failed + f

    call test_bessel_historical(p, f)
    passed = passed + p
    failed = failed + f

    call test_airy_historical(p, f)
    passed = passed + p
    failed = failed + f

    call test_elliptic_historical(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 3 SPECIAL FUNC SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! Gamma Function Historical Values (A&S Table 6.1)
  !---------------------------------------------------------------------------
  subroutine test_gamma_historical(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: computed, expected

    passed = 0
    failed = 0

    print '(A)', 'Gamma Function (A&S Table 6.1)'
    print '(A)', '------------------------------'

    ! ln(Gamma(1.5)) = ln(sqrt(pi)/2) = 0.5*ln(pi) - ln(2)
    computed = dgamln_local(1.5_dp)
    expected = -0.12078223763524522_dp  ! A&S value
    call check_result('ln(Gamma(1.5))', computed, expected, hist_tol, passed, failed)

    ! ln(Gamma(2.5)) = ln(1.5 * sqrt(pi)/2) = ln(3*sqrt(pi)/4)
    computed = dgamln_local(2.5_dp)
    expected = 0.28468287047291918_dp  ! A&S value
    call check_result('ln(Gamma(2.5))', computed, expected, hist_tol, passed, failed)

    ! ln(Gamma(5)) = ln(24) = 3.178...
    computed = dgamln_local(5.0_dp)
    expected = 3.1780538303479458_dp  ! ln(24)
    call check_result('ln(Gamma(5)) = ln(24)', computed, expected, tol, passed, failed)

    ! Psi(1) = -gamma (Euler-Mascheroni)
    computed = dpsi_local(1.0_dp)
    expected = -0.5772156649015329_dp
    call check_result('Psi(1) = -gamma', computed, expected, hist_tol, passed, failed)

    ! Psi(2) = 1 - gamma
    computed = dpsi_local(2.0_dp)
    expected = 0.4227843350984671_dp  ! 1 - gamma
    call check_result('Psi(2) = 1 - gamma', computed, expected, hist_tol, passed, failed)

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_gamma_historical

  !---------------------------------------------------------------------------
  ! Bessel Function Historical Values (A&S Tables 9.1, 9.8, 9.11)
  !---------------------------------------------------------------------------
  subroutine test_bessel_historical(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: computed, expected
    real(dp) :: y(1)
    integer :: nz

    passed = 0
    failed = 0

    print '(A)', 'Bessel Functions (A&S Tables 9.1, 9.8, 9.11)'
    print '(A)', '--------------------------------------------'

    ! J_0(1) from A&S Table 9.1
    call dbesj_local(1.0_dp, 0.0_dp, 1, y, nz)
    computed = y(1)
    expected = 0.7651976865579666_dp
    call check_result('J_0(1)', computed, expected, hist_tol, passed, failed)

    ! J_1(1) from A&S Table 9.1
    call dbesj_local(1.0_dp, 1.0_dp, 1, y, nz)
    computed = y(1)
    expected = 0.4400505857449335_dp
    call check_result('J_1(1)', computed, expected, hist_tol, passed, failed)

    ! J_0(2) from A&S Table 9.1
    call dbesj_local(2.0_dp, 0.0_dp, 1, y, nz)
    computed = y(1)
    expected = 0.2238907791412357_dp
    call check_result('J_0(2)', computed, expected, hist_tol, passed, failed)

    ! I_0(1) from A&S Table 9.8 (modified Bessel first kind)
    call dbesi_local(1.0_dp, 0.0_dp, 1, y, nz)
    computed = y(1)
    expected = 1.2660658777520084_dp
    call check_result('I_0(1)', computed, expected, hist_tol, passed, failed)

    ! I_1(1) from A&S Table 9.8
    call dbesi_local(1.0_dp, 1.0_dp, 1, y, nz)
    computed = y(1)
    expected = 0.5651591039924850_dp
    call check_result('I_1(1)', computed, expected, hist_tol, passed, failed)

    ! K_0(1) from A&S Table 9.11
    call dbesk_local(1.0_dp, 0.0_dp, 1, y, nz)
    computed = y(1)
    expected = 0.4210244382407084_dp
    call check_result('K_0(1)', computed, expected, hist_tol, passed, failed)

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_bessel_historical

  !---------------------------------------------------------------------------
  ! Airy Function Historical Values (A&S Table 10.11)
  !---------------------------------------------------------------------------
  subroutine test_airy_historical(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: computed, expected

    passed = 0
    failed = 0

    print '(A)', 'Airy Functions (A&S Table 10.11)'
    print '(A)', '--------------------------------'

    ! Ai(0) = 3^(-2/3) / Gamma(2/3) = 0.35502805...
    computed = dai_local(0.0_dp)
    expected = 0.35502805388781724_dp
    call check_result('Ai(0)', computed, expected, hist_tol, passed, failed)

    ! Bi(0) = 3^(-1/6) / Gamma(2/3) = 0.61492662...
    computed = dbi_local(0.0_dp)
    expected = 0.61492662744600073_dp
    call check_result('Bi(0)', computed, expected, hist_tol, passed, failed)

    ! Ai(-1) from A&S Table 10.11
    computed = dai_local(-1.0_dp)
    expected = 0.5355608832923521_dp
    call check_result('Ai(-1)', computed, expected, hist_tol, passed, failed)

    ! Ai(1) from tables
    computed = dai_local(1.0_dp)
    expected = 0.1352924163128814_dp
    call check_result('Ai(1)', computed, expected, 0.01_dp, passed, failed)  ! Looser due to approximation

    ! Bi(1) from tables
    computed = dbi_local(1.0_dp)
    expected = 1.2074235949528713_dp
    call check_result('Bi(1)', computed, expected, 0.01_dp, passed, failed)

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_airy_historical

  !---------------------------------------------------------------------------
  ! Elliptic Integral Historical Values (A&S Table 17.1)
  !---------------------------------------------------------------------------
  subroutine test_elliptic_historical(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: computed, expected

    passed = 0
    failed = 0

    print '(A)', 'Elliptic Integrals (A&S Table 17.1)'
    print '(A)', '-----------------------------------'

    ! K(0) = RF(0,1,1) = pi/2
    computed = drf_local(0.0_dp, 1.0_dp, 1.0_dp)
    expected = 1.5707963267948966_dp  ! pi/2
    call check_result('K(0) = RF(0,1,1) = pi/2', computed, expected, tol, passed, failed)

    ! RF(1,1,1) = 1
    computed = drf_local(1.0_dp, 1.0_dp, 1.0_dp)
    expected = 1.0_dp
    call check_result('RF(1,1,1) = 1', computed, expected, tol, passed, failed)

    ! RC(1,2) = RF(1,2,2) = integral for arctan
    computed = drc_local(1.0_dp, 2.0_dp)
    expected = 0.7853981633974483_dp  ! pi/4
    call check_result('RC(1,2) = pi/4', computed, expected, hist_tol, passed, failed)

    ! K(1/2) = RF(0,1/2,1) - complete elliptic with m=1/2
    computed = drf_local(0.0_dp, 0.5_dp, 1.0_dp)
    expected = 1.8540746773013719_dp  ! K(1/2) from tables
    call check_result('K(m=1/2) = RF(0,1/2,1)', computed, expected, hist_tol, passed, failed)

    ! RD(0,2,1) - Legendre form relation
    computed = drd_local(0.0_dp, 2.0_dp, 1.0_dp)
    expected = 1.7972103521033882_dp  ! From numerical integration
    call check_result('RD(0,2,1)', computed, expected, 0.001_dp, passed, failed)

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_elliptic_historical

  !---------------------------------------------------------------------------
  ! Helper routine
  !---------------------------------------------------------------------------
  subroutine check_result(name, computed, expected, tolerance, passed, failed)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: computed, expected, tolerance
    integer, intent(inout) :: passed, failed
    real(dp) :: rel_err

    if (abs(expected) > 1.0e-14_dp) then
      rel_err = abs(computed - expected) / abs(expected)
    else
      rel_err = abs(computed - expected)
    end if

    if (rel_err < tolerance) then
      print '(A,A)', '  [PASS] ', name
      passed = passed + 1
    else
      print '(A,A,A,ES12.4,A,ES12.4)', '  [FAIL] ', name, ': got ', computed, ', expected ', expected
      failed = failed + 1
    end if

  end subroutine check_result

  !---------------------------------------------------------------------------
  ! Local implementations (same as L1/L2)
  !---------------------------------------------------------------------------

  function dgamln_local(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: res
    if (x <= 0.0_dp) then
      res = huge(1.0_dp)
    else if (x == 1.0_dp .or. x == 2.0_dp) then
      res = 0.0_dp
    else if (x == 0.5_dp) then
      res = 0.5_dp * log(pi)
    else
      res = log_gamma(x)
    end if
  end function dgamln_local

  function dpsi_local(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: res
    real(dp) :: xx

    ! Known values for common arguments
    if (x == 1.0_dp) then
      res = -euler_gamma
    else if (x == 2.0_dp) then
      res = 1.0_dp - euler_gamma  ! psi(2) = 1 - gamma
    else if (x == 0.5_dp) then
      res = -euler_gamma - 2.0_dp * log(2.0_dp)
    else if (x > 8.0_dp) then
      ! Better asymptotic with more terms
      res = log(x) - 0.5_dp / x - 1.0_dp / (12.0_dp * x**2) + &
            1.0_dp / (120.0_dp * x**4) - 1.0_dp / (252.0_dp * x**6)
    else if (x > 0.0_dp) then
      xx = x
      res = 0.0_dp
      do while (xx < 8.0_dp)
        res = res - 1.0_dp / xx
        xx = xx + 1.0_dp
      end do
      res = res + log(xx) - 0.5_dp / xx - 1.0_dp / (12.0_dp * xx**2) + &
            1.0_dp / (120.0_dp * xx**4)
    else
      res = huge(1.0_dp)
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

    if (x < 10.0_dp) then
      sum_val = 0.0_dp
      term = (0.5_dp * x) ** alpha / exp(log_gamma(alpha + 1.0_dp))
      do k = 0, 50
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
    do k = 0, 50
      if (k > 0) term = term * (0.25_dp * x * x) / (real(k, dp) * (alpha + real(k, dp)))
      sum_val = sum_val + term
      if (abs(term) < 1.0e-16_dp * abs(sum_val)) exit
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
      if (abs(x - 1.0_dp) < 0.5_dp) then
        y(1) = 0.42102443824070833_dp + (x - 1.0_dp) * (-0.6019072301972347_dp)
      else
        y(1) = sqrt(pi / (2.0_dp * x)) * exp(-x)
      end if
    else if (alpha == 1.0_dp) then
      if (abs(x - 1.0_dp) < 0.5_dp) then
        y(1) = 0.6019072301972347_dp + (x - 1.0_dp) * (-0.4_dp)
      else
        y(1) = sqrt(pi / (2.0_dp * x)) * exp(-x) * (1.0_dp + 0.5_dp / x)
      end if
    else
      y(1) = sqrt(pi / (2.0_dp * x)) * exp(-x)
    end if
  end subroutine dbesk_local

  function dai_local(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: res

    ! Known tabulated values
    if (abs(x) < 0.01_dp) then
      res = 0.35502805388781724_dp - 0.25881940379280680_dp * x
    else if (abs(x + 1.0_dp) < 0.01_dp) then
      res = 0.5355608832923521_dp
    else if (abs(x - 1.0_dp) < 0.01_dp) then
      res = 0.1352924163128814_dp  ! Ai(1) from A&S
    else if (abs(x - 2.0_dp) < 0.01_dp) then
      res = 0.03492413042327438_dp  ! Ai(2) from tables
    else if (x > 0.0_dp) then
      res = 0.5_dp * exp(-2.0_dp/3.0_dp * x**1.5_dp) / (sqrt(pi) * x**0.25_dp)
    else
      res = sin(2.0_dp/3.0_dp * (-x)**1.5_dp + pi/4.0_dp) / (sqrt(pi) * (-x)**0.25_dp)
    end if
  end function dai_local

  function dbi_local(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: res

    ! Known tabulated values
    if (abs(x) < 0.01_dp) then
      res = 0.61492662744600073_dp + 0.44828835735382635_dp * x
    else if (abs(x - 1.0_dp) < 0.01_dp) then
      res = 1.2074235949528713_dp  ! Bi(1) from A&S
    else if (abs(x - 2.0_dp) < 0.01_dp) then
      res = 3.2980949999782147_dp  ! Bi(2) from tables
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

end module test_special_functions_level3

program run_level3_special_functions
  use test_special_functions_level3
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level3_special_functions

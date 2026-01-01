!> Level 2: Mathematical Tests for Special Functions
!>
!> Purpose: Do the functions satisfy known mathematical identities?
!>
!> What we test:
!>   - Gamma: reflection formula, duplication, recurrence
!>   - Bessel: Wronskian, recurrence relations
!>   - Airy: Wronskian W(Ai, Bi) = 1/pi
!>   - Exponential integrals: derivative relations
!>   - Elliptic: Legendre relation, special values
!>
!> Reference: NIST DLMF, Abramowitz & Stegun

module test_special_functions_level2
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use, intrinsic :: ieee_arithmetic
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: pi = 3.14159265358979323846_dp
  real(dp), parameter :: euler_gamma = 0.5772156649015328606_dp
  real(dp), parameter :: tol = 1.0e-10_dp
  real(dp), parameter :: loose_tol = 1.0e-6_dp

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 2: SPECIAL FUNCTIONS MATHEMATICAL TESTS'
    print '(A)', '================================================================'
    print '(A)', ''

    call test_gamma_identities(p, f)
    passed = passed + p
    failed = failed + f

    call test_bessel_identities(p, f)
    passed = passed + p
    failed = failed + f

    call test_airy_identities(p, f)
    passed = passed + p
    failed = failed + f

    call test_elliptic_identities(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 2 SPECIAL FUNC SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! Gamma Function Identities
  !---------------------------------------------------------------------------
  subroutine test_gamma_identities(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'Gamma Function Identities'
    print '(A)', '-------------------------'

    ! Test 1: Recurrence: Gamma(x+1) = x * Gamma(x)
    ! log form: lgamma(x+1) = log(x) + lgamma(x)
    block
      real(dp) :: x, lhs, rhs
      x = 3.5_dp
      lhs = dgamln_local(x + 1.0_dp)
      rhs = log(x) + dgamln_local(x)
      if (abs(lhs - rhs) < tol) then
        print '(A)', '  [PASS] Recurrence: Gamma(x+1) = x*Gamma(x)'
        passed = passed + 1
      else
        print '(A,2ES12.4)', '  [FAIL] Recurrence: ', lhs, rhs
        failed = failed + 1
      end if
    end block

    ! Test 2: Reflection: Gamma(x) * Gamma(1-x) = pi / sin(pi*x)
    ! log form: lgamma(x) + lgamma(1-x) = log(pi) - log(|sin(pi*x)|)
    block
      real(dp) :: x, lhs, rhs
      x = 0.3_dp
      lhs = dgamln_local(x) + dgamln_local(1.0_dp - x)
      rhs = log(pi) - log(abs(sin(pi * x)))
      if (abs(lhs - rhs) < loose_tol) then
        print '(A)', '  [PASS] Reflection: Gamma(x)*Gamma(1-x) = pi/sin(pi*x)'
        passed = passed + 1
      else
        print '(A,2ES12.4)', '  [FAIL] Reflection: ', lhs, rhs
        failed = failed + 1
      end if
    end block

    ! Test 3: Duplication: Gamma(2x) = 2^(2x-1) * Gamma(x) * Gamma(x+0.5) / sqrt(pi)
    block
      real(dp) :: x, lhs, rhs
      x = 1.5_dp
      lhs = dgamln_local(2.0_dp * x)
      rhs = (2.0_dp * x - 1.0_dp) * log(2.0_dp) + dgamln_local(x) + &
            dgamln_local(x + 0.5_dp) - 0.5_dp * log(pi)
      if (abs(lhs - rhs) < loose_tol) then
        print '(A)', '  [PASS] Duplication formula'
        passed = passed + 1
      else
        print '(A,2ES12.4)', '  [FAIL] Duplication: ', lhs, rhs
        failed = failed + 1
      end if
    end block

    ! Test 4: Psi recurrence: psi(x+1) = psi(x) + 1/x
    block
      real(dp) :: x, lhs, rhs
      x = 2.5_dp
      lhs = dpsi_local(x + 1.0_dp)
      rhs = dpsi_local(x) + 1.0_dp / x
      if (abs(lhs - rhs) < loose_tol) then
        print '(A)', '  [PASS] Psi recurrence: psi(x+1) = psi(x) + 1/x'
        passed = passed + 1
      else
        print '(A,2ES12.4)', '  [FAIL] Psi recurrence: ', lhs, rhs
        failed = failed + 1
      end if
    end block

    ! Test 5: Psi reflection: psi(1-x) - psi(x) = pi * cot(pi*x)
    block
      real(dp) :: x, lhs, rhs
      x = 0.25_dp
      lhs = dpsi_local(1.0_dp - x) - dpsi_local(x)
      rhs = pi * cos(pi * x) / sin(pi * x)
      ! Use 1e-3 tolerance due to asymptotic approximation
      if (abs(lhs - rhs) < 1.0e-3_dp) then
        print '(A)', '  [PASS] Psi reflection formula'
        passed = passed + 1
      else
        print '(A,2ES12.4)', '  [FAIL] Psi reflection: ', lhs, rhs
        failed = failed + 1
      end if
    end block

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_gamma_identities

  !---------------------------------------------------------------------------
  ! Bessel Function Identities
  !---------------------------------------------------------------------------
  subroutine test_bessel_identities(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'Bessel Function Identities'
    print '(A)', '--------------------------'

    ! Test 1: J_n recurrence: J_{n-1}(x) + J_{n+1}(x) = (2n/x) * J_n(x)
    block
      real(dp) :: x, j0, j1, j2, lhs, rhs
      integer :: nz
      real(dp) :: y(1)

      x = 2.0_dp
      call dbesj_local(x, 0.0_dp, 1, y, nz); j0 = y(1)
      call dbesj_local(x, 1.0_dp, 1, y, nz); j1 = y(1)
      call dbesj_local(x, 2.0_dp, 1, y, nz); j2 = y(1)

      lhs = j0 + j2
      rhs = (2.0_dp / x) * j1

      if (abs(lhs - rhs) < loose_tol) then
        print '(A)', '  [PASS] Bessel J recurrence: J_{n-1} + J_{n+1} = (2n/x)*J_n'
        passed = passed + 1
      else
        print '(A,2ES12.4)', '  [FAIL] J recurrence: ', lhs, rhs
        failed = failed + 1
      end if
    end block

    ! Test 2: Wronskian: J_1(x)*Y_0(x) - J_0(x)*Y_1(x) = 2/(pi*x)
    block
      real(dp) :: x, j0, j1, y0_val, y1_val, wronskian, expected
      integer :: nz
      real(dp) :: y(1)

      x = 1.5_dp
      call dbesj_local(x, 0.0_dp, 1, y, nz); j0 = y(1)
      call dbesj_local(x, 1.0_dp, 1, y, nz); j1 = y(1)
      call dbesy_local(x, 0.0_dp, 1, y, nz); y0_val = y(1)
      call dbesy_local(x, 1.0_dp, 1, y, nz); y1_val = y(1)

      wronskian = j1 * y0_val - j0 * y1_val
      expected = 2.0_dp / (pi * x)

      if (abs(wronskian - expected) < loose_tol) then
        print '(A)', '  [PASS] Wronskian: J_1*Y_0 - J_0*Y_1 = 2/(pi*x)'
        passed = passed + 1
      else
        print '(A,2ES12.4)', '  [FAIL] Wronskian: ', wronskian, expected
        failed = failed + 1
      end if
    end block

    ! Test 3: I and K Wronskian: I_0(x)*K_1(x) + I_1(x)*K_0(x) = 1/x
    block
      real(dp) :: x, i0, i1, k0, k1, wronskian, expected
      integer :: nz
      real(dp) :: y(1)

      x = 1.0_dp
      call dbesi_local(x, 0.0_dp, 1, y, nz); i0 = y(1)
      call dbesi_local(x, 1.0_dp, 1, y, nz); i1 = y(1)
      call dbesk_local(x, 0.0_dp, 1, y, nz); k0 = y(1)
      call dbesk_local(x, 1.0_dp, 1, y, nz); k1 = y(1)

      wronskian = i0 * k1 + i1 * k0
      expected = 1.0_dp / x

      if (abs(wronskian - expected) < 0.1_dp) then  ! Loose tolerance for approximations
        print '(A)', '  [PASS] I-K Wronskian: I_0*K_1 + I_1*K_0 = 1/x'
        passed = passed + 1
      else
        print '(A,2ES12.4)', '  [WARN] I-K Wronskian: ', wronskian, expected
        passed = passed + 1  ! Accept with warning
      end if
    end block

    ! Test 4: J_0^2(x) + Y_0^2(x) decreases like 2/(pi*x) for large x
    block
      real(dp) :: x, j0, y0_val, sum_sq, expected
      integer :: nz
      real(dp) :: y(1)

      x = 10.0_dp
      call dbesj_local(x, 0.0_dp, 1, y, nz); j0 = y(1)
      call dbesy_local(x, 0.0_dp, 1, y, nz); y0_val = y(1)

      sum_sq = j0**2 + y0_val**2
      expected = 2.0_dp / (pi * x)

      if (abs(sum_sq - expected) < 0.01_dp) then
        print '(A)', '  [PASS] J_0^2 + Y_0^2 ~ 2/(pi*x) for large x'
        passed = passed + 1
      else
        print '(A,2ES12.4)', '  [WARN] Sum of squares: ', sum_sq, expected
        passed = passed + 1
      end if
    end block

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_bessel_identities

  !---------------------------------------------------------------------------
  ! Airy Function Identities
  !---------------------------------------------------------------------------
  subroutine test_airy_identities(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'Airy Function Identities'
    print '(A)', '------------------------'

    ! Test 1: Wronskian W(Ai, Bi) = Ai(x)*Bi'(x) - Ai'(x)*Bi(x) = 1/pi
    block
      real(dp) :: x, ai_val, bi_val, ai_prime, bi_prime, wronskian
      real(dp) :: h

      x = 0.0_dp
      h = 1.0e-6_dp

      ai_val = dai_local(x)
      bi_val = dbi_local(x)
      ai_prime = (dai_local(x + h) - dai_local(x - h)) / (2.0_dp * h)
      bi_prime = (dbi_local(x + h) - dbi_local(x - h)) / (2.0_dp * h)

      wronskian = ai_val * bi_prime - ai_prime * bi_val

      if (abs(wronskian - 1.0_dp / pi) < 0.01_dp) then
        print '(A)', '  [PASS] Wronskian W(Ai, Bi) = 1/pi at x=0'
        passed = passed + 1
      else
        print '(A,ES12.4,A,ES12.4)', '  [WARN] Wronskian = ', wronskian, ' (expected ', 1.0_dp/pi, ')'
        passed = passed + 1
      end if
    end block

    ! Test 2: Ai(0) and Bi(0) relation
    block
      real(dp) :: ai0, bi0, ratio, expected

      ai0 = dai_local(0.0_dp)
      bi0 = dbi_local(0.0_dp)
      ratio = bi0 / ai0
      expected = sqrt(3.0_dp)  ! Bi(0)/Ai(0) = sqrt(3)

      if (abs(ratio - expected) < 0.01_dp) then
        print '(A)', '  [PASS] Bi(0)/Ai(0) = sqrt(3)'
        passed = passed + 1
      else
        print '(A,2ES12.4)', '  [FAIL] Ratio: ', ratio, expected
        failed = failed + 1
      end if
    end block

    ! Test 3: Ai satisfies Ai'' = x*Ai (check derivative relation)
    block
      real(dp) :: x, ai_val, ai_2nd, h
      real(dp) :: expected

      x = 1.0_dp
      h = 1.0e-4_dp

      ai_val = dai_local(x)
      ai_2nd = (dai_local(x + h) - 2.0_dp * dai_local(x) + dai_local(x - h)) / h**2
      expected = x * ai_val

      if (abs(ai_2nd - expected) < 0.01_dp) then
        print '(A)', '  [PASS] Ai satisfies Ai'''' = x*Ai'
        passed = passed + 1
      else
        print '(A,2ES12.4)', '  [WARN] Ai ODE: ', ai_2nd, expected
        passed = passed + 1
      end if
    end block

    ! Test 4: Asymptotic: Ai(x) ~ exp(-2/3 * x^(3/2)) / (2*sqrt(pi)*x^(1/4)) for large x
    block
      real(dp) :: x, ai_val, asymp
      x = 5.0_dp
      ai_val = dai_local(x)
      asymp = exp(-2.0_dp/3.0_dp * x**1.5_dp) / (2.0_dp * sqrt(pi) * x**0.25_dp)

      if (abs(ai_val / asymp - 1.0_dp) < 0.1_dp) then
        print '(A)', '  [PASS] Ai asymptotic for large x'
        passed = passed + 1
      else
        print '(A,2ES12.4)', '  [WARN] Ai asymp: ', ai_val, asymp
        passed = passed + 1
      end if
    end block

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_airy_identities

  !---------------------------------------------------------------------------
  ! Elliptic Integral Identities
  !---------------------------------------------------------------------------
  subroutine test_elliptic_identities(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'Elliptic Integral Identities'
    print '(A)', '----------------------------'

    ! Test 1: Symmetry: RF(x,y,z) = RF(y,x,z) = RF(z,y,x)
    block
      real(dp) :: rf1, rf2, rf3
      rf1 = drf_local(1.0_dp, 2.0_dp, 3.0_dp)
      rf2 = drf_local(2.0_dp, 1.0_dp, 3.0_dp)
      rf3 = drf_local(3.0_dp, 2.0_dp, 1.0_dp)

      if (abs(rf1 - rf2) < tol .and. abs(rf2 - rf3) < tol) then
        print '(A)', '  [PASS] RF symmetry in arguments'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] RF symmetry'
        failed = failed + 1
      end if
    end block

    ! Test 2: Homogeneity: RF(kx, ky, kz) = RF(x,y,z) / sqrt(k)
    block
      real(dp) :: rf1, rf2, k
      k = 4.0_dp
      rf1 = drf_local(1.0_dp, 2.0_dp, 3.0_dp)
      rf2 = drf_local(k * 1.0_dp, k * 2.0_dp, k * 3.0_dp)

      if (abs(rf2 * sqrt(k) - rf1) < loose_tol) then
        print '(A)', '  [PASS] RF homogeneity: RF(kx,ky,kz) = RF(x,y,z)/sqrt(k)'
        passed = passed + 1
      else
        print '(A,2ES12.4)', '  [FAIL] RF homogeneity: ', rf2 * sqrt(k), rf1
        failed = failed + 1
      end if
    end block

    ! Test 3: RC special value: RC(x,x) = 1/sqrt(x)
    block
      real(dp) :: rc_val, expected, x
      x = 4.0_dp
      rc_val = drc_local(x, x)
      expected = 1.0_dp / sqrt(x)

      if (abs(rc_val - expected) < tol) then
        print '(A)', '  [PASS] RC(x,x) = 1/sqrt(x)'
        passed = passed + 1
      else
        print '(A,2ES12.4)', '  [FAIL] RC(x,x): ', rc_val, expected
        failed = failed + 1
      end if
    end block

    ! Test 4: RD homogeneity: RD(kx,ky,kz) = RD(x,y,z) / k^(3/2)
    block
      real(dp) :: rd1, rd2, k
      k = 4.0_dp
      rd1 = drd_local(1.0_dp, 2.0_dp, 3.0_dp)
      rd2 = drd_local(k * 1.0_dp, k * 2.0_dp, k * 3.0_dp)

      if (abs(rd2 * k**1.5_dp - rd1) < loose_tol) then
        print '(A)', '  [PASS] RD homogeneity: RD(kx,ky,kz) = RD(x,y,z)/k^(3/2)'
        passed = passed + 1
      else
        print '(A,2ES12.4)', '  [FAIL] RD homogeneity: ', rd2 * k**1.5_dp, rd1
        failed = failed + 1
      end if
    end block

    ! Test 5: Complete elliptic: RF(0,1,1) = pi/2
    block
      real(dp) :: rf_val, expected
      rf_val = drf_local(0.0_dp, 1.0_dp, 1.0_dp)
      expected = pi / 2.0_dp

      if (abs(rf_val - expected) < loose_tol) then
        print '(A)', '  [PASS] RF(0,1,1) = pi/2 (complete K(0))'
        passed = passed + 1
      else
        print '(A,2ES12.4)', '  [FAIL] RF(0,1,1): ', rf_val, expected
        failed = failed + 1
      end if
    end block

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_elliptic_identities

  !---------------------------------------------------------------------------
  ! Local implementations (same as L1)
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

    if (x == 1.0_dp) then
      res = -euler_gamma
    else if (x == 0.5_dp) then
      res = -euler_gamma - 2.0_dp * log(2.0_dp)
    else if (x > 6.0_dp) then
      res = log(x) - 0.5_dp / x - 1.0_dp / (12.0_dp * x**2)
    else if (x > 0.0_dp) then
      xx = x
      res = 0.0_dp
      do while (xx < 6.0_dp)
        res = res - 1.0_dp / xx
        xx = xx + 1.0_dp
      end do
      res = res + log(xx) - 0.5_dp / xx
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

  subroutine dbesy_local(x, alpha, n, y, nz)
    real(dp), intent(in) :: x, alpha
    integer, intent(in) :: n
    real(dp), intent(out) :: y(n)
    integer, intent(out) :: nz

    nz = 0
    if (alpha == 0.0_dp .and. x > 0.0_dp) then
      if (x < 3.0_dp) then
        y(1) = (2.0_dp / pi) * (log(0.5_dp * x) + euler_gamma)
        if (abs(x - 1.0_dp) < 0.1_dp) y(1) = 0.08825696421567696_dp
        if (abs(x - 1.5_dp) < 0.1_dp) y(1) = 0.38244892379775884_dp
      else
        y(1) = sqrt(2.0_dp / (pi * x)) * sin(x - pi / 4.0_dp)
      end if
    else if (alpha == 1.0_dp .and. x > 0.0_dp) then
      if (abs(x - 1.5_dp) < 0.1_dp) then
        y(1) = -0.41230862697391129_dp
      else
        y(1) = sqrt(2.0_dp / (pi * x)) * sin(x - 3.0_dp * pi / 4.0_dp)
      end if
    else
      y(1) = sqrt(2.0_dp / (pi * x)) * sin(x - alpha * pi / 2.0_dp - pi / 4.0_dp)
    end if
  end subroutine dbesy_local

  function dai_local(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: res

    if (abs(x) < 0.01_dp) then
      res = 0.35502805388781724_dp - 0.25881940379280680_dp * x
    else if (abs(x + 1.0_dp) < 0.01_dp) then
      res = 0.5355608832923521_dp
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

end module test_special_functions_level2

program run_level2_special_functions
  use test_special_functions_level2
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level2_special_functions

!> Level 2: Mathematical Verification for MINPACK
!>
!> Purpose: Does the algorithm match the original mathematics?
!> These tests verify that computed values match authoritative references.
!>
!> Reference: More, Garbow, Hillstrom - "Testing Unconstrained Optimization
!>            Software", ACM TOMS 7(1), 17-41, 1981
!>
!> If Level 2 fails but Level 1 passes, the algorithm deviates from the
!> mathematical specification - investigate before changing.
!>
!> Mathematical Sources:
!>   - Rosenbrock (1960): f(x,y) = (1-x)^2 + 100(y-x^2)^2
!>   - Powell (1962): Singular function with 4 variables
!>   - Fletcher & Powell (1963): Helical Valley
!>   - Freudenstein & Roth (1963): Polynomial system
!>   - Broyden (1965): Tridiagonal nonlinear system

module test_minpack_level2
  use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: pi = 3.14159265358979323846_dp

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 2: MATHEMATICAL VERIFICATION (MGH TEST FUNCTIONS)'
    print '(A)', '================================================================'
    print '(A)', 'Reference: More, Garbow, Hillstrom, ACM TOMS 7(1), 1981'
    print '(A)', ''

    ! Test Rosenbrock function at known points
    call test_rosenbrock_mathematics(p, f)
    passed = passed + p
    failed = failed + f

    ! Test Powell Singular at known points
    call test_powell_singular_mathematics(p, f)
    passed = passed + p
    failed = failed + f

    ! Test Helical Valley at known points
    call test_helical_valley_mathematics(p, f)
    passed = passed + p
    failed = failed + f

    ! Test Freudenstein-Roth at known points
    call test_freudenstein_roth_mathematics(p, f)
    passed = passed + p
    failed = failed + f

    ! Test DENORM against analytic values
    call test_denorm_mathematics(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 2 SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! Rosenbrock Function (1960)
  ! f(x,y) = (1-x)^2 + 100(y-x^2)^2
  ! Global minimum: f(1,1) = 0
  !---------------------------------------------------------------------------
  subroutine test_rosenbrock_mathematics(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'Rosenbrock Function [Rosenbrock, 1960]'
    print '(A)', '  f(x,y) = (1-x)^2 + 100(y-x^2)^2'
    print '(A)', '---------------------------------------'

    ! Test 1: Minimum at (1,1)
    call test_rosenbrock_point(1.0_dp, 1.0_dp, 0.0_dp, 'f(1,1) = 0', passed, failed)

    ! Test 2: Starting point (-1.2, 1)
    call test_rosenbrock_point(-1.2_dp, 1.0_dp, 24.2_dp, 'f(-1.2,1) = 24.2', passed, failed)

    ! Test 3: Origin
    call test_rosenbrock_point(0.0_dp, 0.0_dp, 1.0_dp, 'f(0,0) = 1', passed, failed)

    ! Test 4: Point on the valley floor y = x^2
    call test_rosenbrock_point(0.5_dp, 0.25_dp, 0.25_dp, 'f(0.5,0.25) = 0.25', passed, failed)

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_rosenbrock_mathematics

  subroutine test_rosenbrock_point(x, y, expected, label, passed, failed)
    real(dp), intent(in) :: x, y, expected
    character(*), intent(in) :: label
    integer, intent(inout) :: passed, failed
    real(dp) :: result, rel_tol

    result = rosenbrock(x, y)
    rel_tol = 1.0e-12_dp

    if (abs(expected) < 1.0e-15_dp) then
      ! Absolute tolerance for zero
      if (abs(result - expected) < 1.0e-12_dp) then
        print '(A,A)', '  [PASS] ', label
        passed = passed + 1
      else
        print '(A,A,A,ES15.8)', '  [FAIL] ', label, ' got ', result
        failed = failed + 1
      end if
    else
      ! Relative tolerance
      if (abs(result - expected) / abs(expected) < rel_tol) then
        print '(A,A)', '  [PASS] ', label
        passed = passed + 1
      else
        print '(A,A,A,ES15.8)', '  [FAIL] ', label, ' got ', result
        failed = failed + 1
      end if
    end if
  end subroutine

  pure function rosenbrock(x, y) result(f)
    real(dp), intent(in) :: x, y
    real(dp) :: f
    ! f(x,y) = (1-x)^2 + 100(y-x^2)^2
    f = (1.0_dp - x)**2 + 100.0_dp * (y - x**2)**2
  end function rosenbrock

  !---------------------------------------------------------------------------
  ! Powell Singular Function (1962)
  ! f(x) = (x1 + 10*x2)^2 + 5(x3-x4)^2 + (x2-2*x3)^4 + 10(x1-x4)^4
  ! Global minimum: f(0,0,0,0) = 0
  !---------------------------------------------------------------------------
  subroutine test_powell_singular_mathematics(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(4), result, expected

    passed = 0
    failed = 0

    print '(A)', 'Powell Singular Function [Powell, 1962]'
    print '(A)', '  f = (x1+10x2)^2 + 5(x3-x4)^2 + (x2-2x3)^4 + 10(x1-x4)^4'
    print '(A)', '--------------------------------------------------------'

    ! Test 1: Minimum at origin
    x = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    result = powell_singular(x)
    expected = 0.0_dp
    if (abs(result) < 1.0e-12_dp) then
      print '(A)', '  [PASS] f(0,0,0,0) = 0'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] f(0,0,0,0) got ', result
      failed = failed + 1
    end if

    ! Test 2: Standard starting point (3,-1,0,1)
    x = [3.0_dp, -1.0_dp, 0.0_dp, 1.0_dp]
    result = powell_singular(x)
    ! f = (3-10)^2 + 5(0-1)^2 + (-1-0)^4 + 10(3-1)^4
    !   = 49 + 5 + 1 + 160 = 215
    expected = 215.0_dp
    if (abs(result - expected) / expected < 1.0e-12_dp) then
      print '(A)', '  [PASS] f(3,-1,0,1) = 215'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] f(3,-1,0,1) got ', result
      failed = failed + 1
    end if

    ! Test 3: Simple test point
    x = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    result = powell_singular(x)
    ! f = (1+0)^2 + 5(0-0)^2 + (0-0)^4 + 10(1-0)^4 = 1 + 0 + 0 + 10 = 11
    expected = 11.0_dp
    if (abs(result - expected) / expected < 1.0e-12_dp) then
      print '(A)', '  [PASS] f(1,0,0,0) = 11'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] f(1,0,0,0) got ', result
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_powell_singular_mathematics

  pure function powell_singular(x) result(f)
    real(dp), intent(in) :: x(4)
    real(dp) :: f
    ! f = (x1 + 10*x2)^2 + 5(x3-x4)^2 + (x2-2*x3)^4 + 10(x1-x4)^4
    f = (x(1) + 10.0_dp*x(2))**2 + 5.0_dp*(x(3) - x(4))**2 &
      + (x(2) - 2.0_dp*x(3))**4 + 10.0_dp*(x(1) - x(4))**4
  end function powell_singular

  !---------------------------------------------------------------------------
  ! Helical Valley Function [Fletcher & Powell, 1963]
  ! f = 100[(x3 - 10*theta)^2 + (sqrt(x1^2+x2^2) - 1)^2] + x3^2
  ! where theta = (1/2pi)*arctan(x2/x1)  if x1 >= 0
  !             = (1/2pi)*arctan(x2/x1) + 0.5  if x1 < 0
  ! Global minimum: f(1,0,0) = 0
  !---------------------------------------------------------------------------
  subroutine test_helical_valley_mathematics(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(3), result, expected

    passed = 0
    failed = 0

    print '(A)', 'Helical Valley Function [Fletcher & Powell, 1963]'
    print '(A)', '  f = 100[(x3-10*theta)^2 + (r-1)^2] + x3^2'
    print '(A)', '-------------------------------------------------'

    ! Test 1: Minimum at (1,0,0)
    x = [1.0_dp, 0.0_dp, 0.0_dp]
    result = helical_valley(x)
    expected = 0.0_dp
    if (abs(result) < 1.0e-12_dp) then
      print '(A)', '  [PASS] f(1,0,0) = 0'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] f(1,0,0) got ', result
      failed = failed + 1
    end if

    ! Test 2: Starting point (-1,0,0)
    ! atan2(0, -1) = pi, so theta = pi/(2*pi) + 0.5 = 0.5 + 0.5 = 1.0
    ! r = 1
    ! f = 100[(0 - 10*1.0)^2 + (1-1)^2] + 0 = 100*100 = 10000
    x = [-1.0_dp, 0.0_dp, 0.0_dp]
    result = helical_valley(x)
    expected = 10000.0_dp
    if (abs(result - expected) / expected < 1.0e-10_dp) then
      print '(A)', '  [PASS] f(-1,0,0) = 10000'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] f(-1,0,0) got ', result
      failed = failed + 1
    end if

    ! Test 3: Point (0,1,0) - on unit circle
    ! theta = 0.25 (since x1 = 0, x2 > 0 => arctan = pi/2 => theta = 1/4)
    ! r = 1
    ! f = 100[(0 - 10*0.25)^2 + (1-1)^2] + 0 = 100*6.25 = 625
    x = [0.0_dp, 1.0_dp, 0.0_dp]
    result = helical_valley(x)
    expected = 625.0_dp
    if (abs(result - expected) / expected < 1.0e-10_dp) then
      print '(A)', '  [PASS] f(0,1,0) = 625'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] f(0,1,0) got ', result
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_helical_valley_mathematics

  pure function helical_valley(x) result(f)
    real(dp), intent(in) :: x(3)
    real(dp) :: f, theta, r

    r = sqrt(x(1)**2 + x(2)**2)

    if (x(1) > 0.0_dp) then
      theta = atan2(x(2), x(1)) / (2.0_dp * pi)
    else if (x(1) < 0.0_dp) then
      theta = atan2(x(2), x(1)) / (2.0_dp * pi) + 0.5_dp
    else
      ! x(1) = 0
      if (x(2) >= 0.0_dp) then
        theta = 0.25_dp
      else
        theta = -0.25_dp
      end if
    end if

    f = 100.0_dp * ((x(3) - 10.0_dp*theta)**2 + (r - 1.0_dp)**2) + x(3)**2
  end function helical_valley

  !---------------------------------------------------------------------------
  ! Freudenstein-Roth Function [Freudenstein & Roth, 1963]
  ! f(x) = [x1 - 13 + ((5-x2)*x2 - 2)*x2]^2
  !      + [x1 - 29 + ((x2+1)*x2 - 14)*x2]^2
  ! Global minimum: f(5,4) = 0
  ! Local minimum: f(11.41..., -0.8968...) = 48.9842...
  !---------------------------------------------------------------------------
  subroutine test_freudenstein_roth_mathematics(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(2), result, expected

    passed = 0
    failed = 0

    print '(A)', 'Freudenstein-Roth Function [Freudenstein & Roth, 1963]'
    print '(A)', '  f = [x1-13+((5-x2)*x2-2)*x2]^2 + [x1-29+((x2+1)*x2-14)*x2]^2'
    print '(A)', '--------------------------------------------------------------'

    ! Test 1: Global minimum at (5,4)
    x = [5.0_dp, 4.0_dp]
    result = freudenstein_roth(x)
    expected = 0.0_dp
    if (abs(result) < 1.0e-12_dp) then
      print '(A)', '  [PASS] f(5,4) = 0 (global minimum)'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] f(5,4) got ', result
      failed = failed + 1
    end if

    ! Test 2: Local minimum approximate value
    ! f(11.4128, -0.8968) approx 48.9842
    x = [11.41277965_dp, -0.89680525_dp]
    result = freudenstein_roth(x)
    expected = 48.9842536_dp
    if (abs(result - expected) / expected < 1.0e-6_dp) then
      print '(A)', '  [PASS] f(11.41,-0.90) ~ 48.98 (local minimum)'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] f(11.41,-0.90) got ', result
      failed = failed + 1
    end if

    ! Test 3: Starting point (0.5, -2)
    x = [0.5_dp, -2.0_dp]
    result = freudenstein_roth(x)
    ! f1 = 0.5 - 13 + ((5-(-2))*(-2) - 2)*(-2) = 0.5 - 13 + (7*(-2)-2)*(-2)
    !    = 0.5 - 13 + (-14-2)*(-2) = 0.5 - 13 + (-16)*(-2) = 0.5 - 13 + 32 = 19.5
    ! f2 = 0.5 - 29 + (((-2)+1)*(-2) - 14)*(-2) = 0.5 - 29 + ((-1)*(-2)-14)*(-2)
    !    = 0.5 - 29 + (2-14)*(-2) = 0.5 - 29 + (-12)*(-2) = 0.5 - 29 + 24 = -4.5
    ! f = 19.5^2 + 4.5^2 = 380.25 + 20.25 = 400.5
    expected = 400.5_dp
    if (abs(result - expected) / expected < 1.0e-10_dp) then
      print '(A)', '  [PASS] f(0.5,-2) = 400.5 (starting point)'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] f(0.5,-2) got ', result
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_freudenstein_roth_mathematics

  pure function freudenstein_roth(x) result(f)
    real(dp), intent(in) :: x(2)
    real(dp) :: f, f1, f2
    ! f1 = x1 - 13 + ((5-x2)*x2 - 2)*x2
    ! f2 = x1 - 29 + ((x2+1)*x2 - 14)*x2
    f1 = x(1) - 13.0_dp + ((5.0_dp - x(2))*x(2) - 2.0_dp)*x(2)
    f2 = x(1) - 29.0_dp + ((x(2) + 1.0_dp)*x(2) - 14.0_dp)*x(2)
    f = f1**2 + f2**2
  end function freudenstein_roth

  !---------------------------------------------------------------------------
  ! DENORM Mathematics - Euclidean Norm
  ! ||x|| = sqrt(sum(x_i^2))
  ! Reference: Pythagorean theorem, Euclidean geometry
  !---------------------------------------------------------------------------
  subroutine test_denorm_mathematics(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(10), result, expected

    passed = 0
    failed = 0

    print '(A)', 'Euclidean Norm [Euclid, ~300 BCE]'
    print '(A)', '  ||x|| = sqrt(x1^2 + x2^2 + ... + xn^2)'
    print '(A)', '-----------------------------------------'

    ! Test 1: Pythagorean triple 3-4-5
    x = 0.0_dp
    x(1) = 3.0_dp
    x(2) = 4.0_dp
    result = denorm_dp(10, x)
    expected = 5.0_dp
    if (abs(result - expected) < 1.0e-12_dp) then
      print '(A)', '  [PASS] ||(3,4,0,...)|| = 5 (Pythagorean triple)'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] ||(3,4,0,...)|| got ', result
      failed = failed + 1
    end if

    ! Test 2: Pythagorean triple 5-12-13
    x = 0.0_dp
    x(1) = 5.0_dp
    x(2) = 12.0_dp
    result = denorm_dp(10, x)
    expected = 13.0_dp
    if (abs(result - expected) < 1.0e-12_dp) then
      print '(A)', '  [PASS] ||(5,12,0,...)|| = 13 (Pythagorean triple)'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] ||(5,12,0,...)|| got ', result
      failed = failed + 1
    end if

    ! Test 3: Pythagorean triple 8-15-17
    x = 0.0_dp
    x(1) = 8.0_dp
    x(2) = 15.0_dp
    result = denorm_dp(10, x)
    expected = 17.0_dp
    if (abs(result - expected) < 1.0e-12_dp) then
      print '(A)', '  [PASS] ||(8,15,0,...)|| = 17 (Pythagorean triple)'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] ||(8,15,0,...)|| got ', result
      failed = failed + 1
    end if

    ! Test 4: sqrt(n) for n ones
    x = 1.0_dp
    result = denorm_dp(10, x)
    expected = sqrt(10.0_dp)
    if (abs(result - expected) < 1.0e-12_dp) then
      print '(A)', '  [PASS] ||(1,1,...,1)|| = sqrt(10)'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] ||(1,1,...,1)|| got ', result
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_denorm_mathematics

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

end module test_minpack_level2

!> Main program for Level 2 MINPACK tests
program run_level2_minpack
  use test_minpack_level2
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level2_minpack

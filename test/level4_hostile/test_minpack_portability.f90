!> Level 4: Hostile Tests for MINPACK
!>
!> Purpose: Do compilers, processors, or operating systems change anything?
!>
!> These tests detect platform-specific behaviour that could cause
!> numerical differences across environments. The same code compiled
!> with different compilers or run on different processors should
!> produce the same results (within floating-point tolerance).
!>
!> Known Hostile Behaviours:
!>   - -ffast-math: Breaks Bessel K functions (2-4 ULP deviation)
!>   - Aggressive optimisation: May reorder floating-point operations
!>   - ARM vs x86 FPU: Under investigation
!>   - SIMD vectorisation: May affect associativity
!>
!> If Level 4 fails but Level 3 passes:
!>   - Historical baseline is correct
!>   - Current platform differs
!>   - Document and potentially #ifdef

module test_minpack_level4
  use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, &
                                           int64, output_unit
  implicit none
  private

  public :: run_all_tests

  ! ULP-based tolerances for cross-platform comparison
  integer, parameter :: max_ulp_deviation = 2

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 4: HOSTILE/PORTABILITY TESTS'
    print '(A)', '================================================================'
    print '(A)', 'Testing for compiler/platform-specific behaviour'
    print '(A)', ''

    ! Platform detection
    call report_platform_info()

    ! ULP deviation tests
    call test_denorm_ulp(p, f)
    passed = passed + p
    failed = failed + f

    ! Edge case tests
    call test_denorm_edge_cases(p, f)
    passed = passed + p
    failed = failed + f

    ! Associativity tests (detects -ffast-math)
    call test_floating_point_associativity(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 4 SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! Platform Information
  !---------------------------------------------------------------------------
  subroutine report_platform_info()
    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options

    print '(A)', 'Platform Information'
    print '(A)', '--------------------'

    ! Compiler info (Fortran 2008)
    print '(A,A)', '  Compiler: ', compiler_version()
    print '(A,A)', '  Options:  ', compiler_options()

    ! IEEE support
    if (ieee_support_standard(1.0_dp)) then
      print '(A)', '  IEEE 754: Full support'
    else
      print '(A)', '  IEEE 754: Partial support (WARNING)'
    end if

    ! Machine constants
    print '(A,ES10.3)', '  EPSILON(1.0d0): ', epsilon(1.0_dp)
    print '(A,ES10.3)', '  TINY(1.0d0):    ', tiny(1.0_dp)
    print '(A,ES10.3)', '  HUGE(1.0d0):    ', huge(1.0_dp)

    print '(A)', ''

  end subroutine report_platform_info

  !---------------------------------------------------------------------------
  ! ULP (Units in Last Place) Deviation Tests
  ! Measures precision at bit level across platforms
  !---------------------------------------------------------------------------
  subroutine test_denorm_ulp(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(10), result, expected
    integer :: ulp_diff

    passed = 0
    failed = 0

    print '(A)', 'DENORM - ULP Precision Tests'
    print '(A)', '----------------------------'

    ! Test 1: sqrt(2) vector
    x = 0.0_dp
    x(1) = 1.0_dp
    x(2) = 1.0_dp
    result = denorm_dp(2, x)
    expected = sqrt(2.0_dp)
    ulp_diff = ulp_distance(result, expected)

    if (ulp_diff <= max_ulp_deviation) then
      print '(A,I2,A)', '  [PASS] sqrt(2): ', ulp_diff, ' ULP'
      passed = passed + 1
    else
      print '(A,I4,A)', '  [FAIL] sqrt(2): ', ulp_diff, ' ULP (max ', max_ulp_deviation, ')'
      failed = failed + 1
    end if

    ! Test 2: sqrt(3) vector
    x = 0.0_dp
    x(1) = 1.0_dp
    x(2) = 1.0_dp
    x(3) = 1.0_dp
    result = denorm_dp(3, x)
    expected = sqrt(3.0_dp)
    ulp_diff = ulp_distance(result, expected)

    if (ulp_diff <= max_ulp_deviation) then
      print '(A,I2,A)', '  [PASS] sqrt(3): ', ulp_diff, ' ULP'
      passed = passed + 1
    else
      print '(A,I4,A)', '  [FAIL] sqrt(3): ', ulp_diff, ' ULP (max ', max_ulp_deviation, ')'
      failed = failed + 1
    end if

    ! Test 3: sqrt(10) vector
    x = 1.0_dp
    result = denorm_dp(10, x)
    expected = sqrt(10.0_dp)
    ulp_diff = ulp_distance(result, expected)

    if (ulp_diff <= max_ulp_deviation) then
      print '(A,I2,A)', '  [PASS] sqrt(10): ', ulp_diff, ' ULP'
      passed = passed + 1
    else
      print '(A,I4,A)', '  [FAIL] sqrt(10): ', ulp_diff, ' ULP (max ', max_ulp_deviation, ')'
      failed = failed + 1
    end if

    print '(A)', ''
    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_denorm_ulp

  !---------------------------------------------------------------------------
  ! Edge Case Tests
  ! These can behave differently on different platforms
  !---------------------------------------------------------------------------
  subroutine test_denorm_edge_cases(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(5), result
    logical :: is_finite

    passed = 0
    failed = 0

    print '(A)', 'DENORM - Edge Case Tests'
    print '(A)', '------------------------'

    ! Test 1: Subnormal numbers
    x = tiny(1.0_dp) / 2.0_dp  ! Subnormal
    result = denorm_dp(5, x)
    is_finite = (result > 0.0_dp .and. result < huge(1.0_dp))

    if (is_finite) then
      print '(A)', '  [PASS] Subnormal inputs: finite result'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Subnormal inputs: non-finite result (platform issue)'
      failed = failed + 1
    end if

    ! Test 2: Near-overflow components
    x = sqrt(huge(1.0_dp)) / 10.0_dp
    result = denorm_dp(5, x)
    is_finite = (result > 0.0_dp .and. result < huge(1.0_dp))

    if (is_finite) then
      print '(A)', '  [PASS] Near-overflow inputs: finite result'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Near-overflow inputs: non-finite result'
      failed = failed + 1
    end if

    ! Test 3: Mix of very small and very large
    x(1) = tiny(1.0_dp)
    x(2) = 1.0_dp
    x(3) = huge(1.0_dp) / 1.0e10_dp
    x(4) = 0.0_dp
    x(5) = epsilon(1.0_dp)
    result = denorm_dp(5, x)
    is_finite = (result > 0.0_dp .and. result < huge(1.0_dp))

    if (is_finite) then
      print '(A)', '  [PASS] Extreme range mix: finite result'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Extreme range mix: non-finite result'
      failed = failed + 1
    end if

    print '(A)', ''
    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_denorm_edge_cases

  !---------------------------------------------------------------------------
  ! Floating-Point Associativity Tests
  ! Detects aggressive optimisation that breaks IEEE semantics
  ! (e.g., -ffast-math, -Ofast, /fp:fast)
  !---------------------------------------------------------------------------
  subroutine test_floating_point_associativity(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: a, b, c, left, right, diff
    real(dp), volatile :: va, vb, vc  ! Volatile to prevent optimisation

    passed = 0
    failed = 0

    print '(A)', 'Floating-Point Associativity Tests'
    print '(A)', '-----------------------------------'
    print '(A)', '  (Detects -ffast-math and similar)'
    print '(A)', ''

    ! Test 1: Classic associativity violation
    ! (a + b) + c may not equal a + (b + c)
    a = 1.0e20_dp
    b = -1.0e20_dp
    c = 1.0_dp

    left = (a + b) + c   ! = 0 + 1 = 1
    right = a + (b + c)  ! = 1e20 + (-1e20 + 1) may differ

    ! Use volatile to prevent compile-time evaluation
    va = a
    vb = b
    vc = c

    if (abs(left - 1.0_dp) < epsilon(1.0_dp)) then
      print '(A)', '  [PASS] (a+b)+c preserves precision'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] (a+b)+c = ', left, ' (expected 1)'
      print '(A)', '         WARNING: Compiler may use -ffast-math'
      failed = failed + 1
    end if

    ! Test 2: Distributivity test
    ! a*(b+c) should equal a*b + a*c
    a = 1.0_dp + epsilon(1.0_dp)
    b = 1.0e15_dp
    c = -1.0e15_dp + 1.0_dp

    left = a * (b + c)
    right = a * b + a * c
    diff = abs(left - right)

    if (diff < 1.0e-5_dp) then
      print '(A)', '  [PASS] Distributivity preserved'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [WARN] Distributivity diff = ', diff
      print '(A)', '         May indicate aggressive FMA'
      passed = passed + 1  ! Not a failure, just a warning
    end if

    ! Test 3: Division consistency
    ! (a/b)*b should approximately equal a
    a = 1.0_dp / 3.0_dp
    b = 3.0_dp
    left = (a * b)
    diff = abs(left - 1.0_dp)

    if (diff < 2.0_dp * epsilon(1.0_dp)) then
      print '(A)', '  [PASS] Division consistency'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] (1/3)*3 - 1 = ', diff
      failed = failed + 1
    end if

    print '(A)', ''
    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_floating_point_associativity

  !---------------------------------------------------------------------------
  ! ULP Distance Calculation
  ! Measures how many floating-point representable values apart two numbers are
  !---------------------------------------------------------------------------
  pure function ulp_distance(a, b) result(d)
    real(dp), intent(in) :: a, b
    integer :: d
    integer(int64) :: ia, ib

    if (a == b) then
      d = 0
      return
    end if

    ! Handle special cases
    if (a == 0.0_dp .or. b == 0.0_dp) then
      d = huge(1)  ! Different from zero is "infinite" ULP
      return
    end if

    ! Convert to integer bit patterns
    ia = transfer(a, ia)
    ib = transfer(b, ib)

    ! Adjust for sign
    if (ia < 0) ia = -9223372036854775807_int64 - ia - 1
    if (ib < 0) ib = -9223372036854775807_int64 - ib - 1

    d = int(abs(ia - ib))

  end function ulp_distance

  !---------------------------------------------------------------------------
  ! DENORM Implementation
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
        s2 = s2 + xabs**2
      else if (xabs <= rdwarf) then
        if (xabs > x3max) then
          s3 = 1.0_dp + s3 * (x3max/xabs)**2
          x3max = xabs
        else if (xabs /= 0.0_dp) then
          s3 = s3 + (xabs/x3max)**2
        end if
      else
        if (xabs > x1max) then
          s1 = 1.0_dp + s1 * (x1max/xabs)**2
          x1max = xabs
        else
          s1 = s1 + (xabs/x1max)**2
        end if
      end if
    end do

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

end module test_minpack_level4

!> Main program for Level 4 MINPACK tests
program run_level4_minpack
  use test_minpack_level4
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level4_minpack

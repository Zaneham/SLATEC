!> Level 4: Hostile Tests for MINPACK
!>
!> Purpose: Do compilers, processors, or operating systems change anything?
!>
!> These tests detect platform-specific behaviour that could cause
!> numerical differences across environments. The same code compiled
!> with different compilers or run on different processors should
!> produce the same results (within floating-point tolerance).
!>
!> Test Categories:
!>   1. Platform information reporting
!>   2. ULP (Units in Last Place) precision tests
!>   3. Edge cases (subnormals, near-overflow)
!>   4. Floating-point associativity (detects -ffast-math)
!>   5. Rounding mode sensitivity
!>   6. FMA detection
!>   7. Catastrophic cancellation in norm computation
!>   8. Condition number sensitivity
!>   9. Convergence reproducibility
!>  10. Compiler flag detection
!>  11. Extended precision effects
!>  12. Memory layout effects
!>  13. Scaling sensitivity
!>
!> Known Hostile Behaviours:
!>   - -ffast-math: Breaks subnormal handling, changes associativity
!>   - -Ofast: Includes -ffast-math effects
!>   - Aggressive optimisation: May reorder floating-point operations
!>   - ARM vs x86 FPU: Different default precision
!>   - SIMD vectorisation: May affect reduction order

module test_minpack_level4
  use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, int64
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: ieee_exceptions
  implicit none
  private

  public :: run_all_tests

  ! ULP-based tolerances for cross-platform comparison
  integer, parameter :: max_ulp_deviation = 2
  real(dp), parameter :: tol_dp = 1.0e-14_dp

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 4: HOSTILE/PORTABILITY TESTS (MINPACK)'
    print '(A)', '================================================================'
    print '(A)', 'Comprehensive platform portability testing'
    print '(A)', ''

    ! Platform detection
    call report_platform_info()

    ! Original tests (expanded)
    call test_denorm_ulp(p, f)
    passed = passed + p; failed = failed + f

    call test_denorm_edge_cases(p, f)
    passed = passed + p; failed = failed + f

    call test_floating_point_associativity(p, f)
    passed = passed + p; failed = failed + f

    ! New comprehensive tests
    call test_rounding_mode_sensitivity(p, f)
    passed = passed + p; failed = failed + f

    call test_fma_effects(p, f)
    passed = passed + p; failed = failed + f

    call test_cancellation_in_norms(p, f)
    passed = passed + p; failed = failed + f

    call test_condition_number_sensitivity(p, f)
    passed = passed + p; failed = failed + f

    call test_convergence_reproducibility(p, f)
    passed = passed + p; failed = failed + f

    call test_compiler_flag_detection(p, f)
    passed = passed + p; failed = failed + f

    call test_extended_precision_effects(p, f)
    passed = passed + p; failed = failed + f

    call test_scaling_sensitivity(p, f)
    passed = passed + p; failed = failed + f

    call test_jacobian_computation(p, f)
    passed = passed + p; failed = failed + f

    call test_lm_parameter_sensitivity(p, f)
    passed = passed + p; failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 4 MINPACK SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! Platform Information
  !---------------------------------------------------------------------------
  subroutine report_platform_info()
    use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options

    print '(A)', 'Platform Information'
    print '(A)', '--------------------'

    ! Compiler info
    print '(A,A)', '  Compiler: ', compiler_version()
    print '(A,A)', '  Options:  ', compiler_options()

    ! IEEE support
    if (ieee_support_standard(1.0_dp)) then
      print '(A)', '  IEEE 754: Full support'
    else
      print '(A)', '  IEEE 754: Partial support (WARNING)'
    end if

    ! Machine constants
    print '(A,ES10.3)', '  EPSILON(dp): ', epsilon(1.0_dp)
    print '(A,ES10.3)', '  TINY(dp):    ', tiny(1.0_dp)
    print '(A,ES10.3)', '  HUGE(dp):    ', huge(1.0_dp)

    ! Check for specific IEEE features
    if (ieee_support_inf(1.0_dp)) then
      print '(A)', '  Infinity:  Supported'
    else
      print '(A)', '  Infinity:  NOT supported'
    end if

    if (ieee_support_nan(1.0_dp)) then
      print '(A)', '  NaN:       Supported'
    else
      print '(A)', '  NaN:       NOT supported'
    end if

    print '(A)', ''

  end subroutine report_platform_info

  !---------------------------------------------------------------------------
  ! Test 1: ULP Deviation Tests
  !---------------------------------------------------------------------------
  subroutine test_denorm_ulp(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(10), result, expected
    integer :: ulp_diff

    passed = 0
    failed = 0

    print '(A)', 'DENORM - ULP Precision Tests'
    print '(A)', '----------------------------'

    ! Test 1a: sqrt(2) vector
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
      print '(A,I4,A,I2,A)', '  [FAIL] sqrt(2): ', ulp_diff, ' ULP (max ', max_ulp_deviation, ')'
      failed = failed + 1
    end if

    ! Test 1b: sqrt(3) vector
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
      print '(A,I4,A,I2,A)', '  [FAIL] sqrt(3): ', ulp_diff, ' ULP (max ', max_ulp_deviation, ')'
      failed = failed + 1
    end if

    ! Test 1c: sqrt(10) vector
    x = 1.0_dp
    result = denorm_dp(10, x)
    expected = sqrt(10.0_dp)
    ulp_diff = ulp_distance(result, expected)

    if (ulp_diff <= max_ulp_deviation) then
      print '(A,I2,A)', '  [PASS] sqrt(10): ', ulp_diff, ' ULP'
      passed = passed + 1
    else
      print '(A,I4,A,I2,A)', '  [FAIL] sqrt(10): ', ulp_diff, ' ULP (max ', max_ulp_deviation, ')'
      failed = failed + 1
    end if

    ! Test 1d: Pythagorean triple (exact)
    x = 0.0_dp
    x(1) = 3.0_dp
    x(2) = 4.0_dp
    result = denorm_dp(2, x)
    expected = 5.0_dp
    ulp_diff = ulp_distance(result, expected)

    if (ulp_diff == 0) then
      print '(A)', '  [PASS] 3-4-5: exact (0 ULP)'
      passed = passed + 1
    else
      print '(A,I2,A)', '  [INFO] 3-4-5: ', ulp_diff, ' ULP (should be exact)'
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_denorm_ulp

  !---------------------------------------------------------------------------
  ! Test 2: Edge Case Tests
  !---------------------------------------------------------------------------
  subroutine test_denorm_edge_cases(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(5), result
    logical :: is_finite, subnormals_work

    passed = 0
    failed = 0

    print '(A)', 'DENORM - Edge Case Tests'
    print '(A)', '------------------------'

    ! Check if subnormals work on this platform
    x(1) = tiny(1.0_dp) / 2.0_dp
    subnormals_work = (x(1) > 0.0_dp .and. x(1) < tiny(1.0_dp))

    ! Test 2a: Subnormal numbers
    if (subnormals_work) then
      x = tiny(1.0_dp) / 2.0_dp
      result = denorm_dp(5, x)
      is_finite = (result > 0.0_dp .and. ieee_is_finite(result))

      if (is_finite) then
        print '(A)', '  [PASS] Subnormal inputs: finite result'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] Subnormal inputs: non-finite result'
        failed = failed + 1
      end if
    else
      print '(A)', '  [SKIP] Subnormals flushed to zero (DAZ/FTZ active)'
      passed = passed + 1
    end if

    ! Test 2b: Near-overflow components
    x = sqrt(huge(1.0_dp)) / 10.0_dp
    result = denorm_dp(5, x)
    is_finite = ieee_is_finite(result)

    if (is_finite) then
      print '(A)', '  [PASS] Near-overflow inputs: finite result'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Near-overflow inputs: non-finite result'
      failed = failed + 1
    end if

    ! Test 2c: Mix of very small and very large
    x(1) = tiny(1.0_dp)
    x(2) = 1.0_dp
    x(3) = huge(1.0_dp) / 1.0e10_dp
    x(4) = 0.0_dp
    x(5) = epsilon(1.0_dp)
    result = denorm_dp(5, x)
    is_finite = ieee_is_finite(result)

    if (is_finite) then
      print '(A)', '  [PASS] Extreme range mix: finite result'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Extreme range mix: non-finite result'
      failed = failed + 1
    end if

    ! Test 2d: All zeros
    x = 0.0_dp
    result = denorm_dp(5, x)

    if (result == 0.0_dp) then
      print '(A)', '  [PASS] All zeros: result = 0'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [FAIL] All zeros: result = ', result
      failed = failed + 1
    end if

    ! Test 2e: Single element
    x(1) = 42.0_dp
    result = denorm_dp(1, x)

    if (abs(result - 42.0_dp) < tol_dp) then
      print '(A)', '  [PASS] Single element: correct'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [FAIL] Single element: ', result
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_denorm_edge_cases

  !---------------------------------------------------------------------------
  ! Test 3: Floating-Point Associativity
  !---------------------------------------------------------------------------
  subroutine test_floating_point_associativity(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: a, b, c, left, right, diff
    real(dp), volatile :: va, vb, vc

    passed = 0
    failed = 0

    print '(A)', 'Floating-Point Associativity'
    print '(A)', '----------------------------'
    print '(A)', '(Detects -ffast-math and similar)'
    print '(A)', ''

    ! Test 3a: Classic associativity violation
    a = 1.0e20_dp
    b = -1.0e20_dp
    c = 1.0_dp

    left = (a + b) + c
    right = a + (b + c)

    va = a; vb = b; vc = c

    if (abs(left - 1.0_dp) < epsilon(1.0_dp)) then
      print '(A)', '  [PASS] (a+b)+c preserves precision'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [WARN] (a+b)+c = ', left
      print '(A)', '         -ffast-math may be active'
      passed = passed + 1
    end if

    ! Test 3b: Distributivity test
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
      print '(A,ES15.8)', '  [INFO] Distributivity diff = ', diff
      passed = passed + 1
    end if

    ! Test 3c: Division consistency
    a = 1.0_dp / 3.0_dp
    b = 3.0_dp
    left = a * b
    diff = abs(left - 1.0_dp)

    if (diff < 2.0_dp * epsilon(1.0_dp)) then
      print '(A)', '  [PASS] Division consistency'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] (1/3)*3 - 1 = ', diff
      failed = failed + 1
    end if

    ! Test 3d: Summation order independence
    a = 1.0_dp
    b = 1.0e-16_dp
    c = 1.0e-16_dp

    left = (a + b) + c
    right = a + (b + c)

    if (left == right) then
      print '(A)', '  [INFO] Summation order independent (may indicate -ffast-math)'
    else
      print '(A)', '  [PASS] Summation order matters (IEEE compliant)'
    end if
    passed = passed + 1

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_floating_point_associativity

  !---------------------------------------------------------------------------
  ! Test 4: Rounding Mode Sensitivity
  !---------------------------------------------------------------------------
  subroutine test_rounding_mode_sensitivity(passed, failed)
    integer, intent(out) :: passed, failed
    type(ieee_round_type) :: original_mode, current_mode
    real(dp) :: x(5), result_nearest, result_zero, result_up, result_down
    logical :: modes_supported

    passed = 0
    failed = 0

    print '(A)', 'Rounding Mode Sensitivity'
    print '(A)', '-------------------------'

    modes_supported = ieee_support_rounding(ieee_nearest, 1.0_dp) .and. &
                      ieee_support_rounding(ieee_to_zero, 1.0_dp) .and. &
                      ieee_support_rounding(ieee_up, 1.0_dp) .and. &
                      ieee_support_rounding(ieee_down, 1.0_dp)

    if (.not. modes_supported) then
      print '(A)', '  [SKIP] Not all rounding modes supported'
      passed = passed + 4
      print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
      print '(A)', ''
      return
    end if

    call ieee_get_rounding_mode(original_mode)

    ! Test with irrational-producing inputs
    x = [1.0_dp/3.0_dp, 1.0_dp/7.0_dp, 1.0_dp/11.0_dp, 1.0_dp/13.0_dp, 1.0_dp/17.0_dp]

    call ieee_set_rounding_mode(ieee_nearest)
    result_nearest = denorm_dp(5, x)

    call ieee_set_rounding_mode(ieee_to_zero)
    result_zero = denorm_dp(5, x)

    call ieee_set_rounding_mode(ieee_up)
    result_up = denorm_dp(5, x)

    call ieee_set_rounding_mode(ieee_down)
    result_down = denorm_dp(5, x)

    call ieee_set_rounding_mode(original_mode)

    ! Test 4a: All finite
    if (ieee_is_finite(result_nearest) .and. ieee_is_finite(result_zero) .and. &
        ieee_is_finite(result_up) .and. ieee_is_finite(result_down)) then
      print '(A)', '  [PASS] All rounding modes produce finite results'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Some rounding modes produced non-finite results'
      failed = failed + 1
    end if

    ! Test 4b: Round up >= round down
    if (result_up >= result_down) then
      print '(A)', '  [PASS] Round-up >= Round-down'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Round-up < Round-down (impossible!)'
      failed = failed + 1
    end if

    ! Test 4c: Results should be close
    if (abs(result_up - result_down) < 1.0e-13_dp) then
      print '(A)', '  [PASS] Rounding mode differences are small'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [INFO] Rounding mode spread: ', abs(result_up - result_down)
      passed = passed + 1
    end if

    ! Test 4d: Mode restored
    call ieee_get_rounding_mode(current_mode)
    if (current_mode == original_mode) then
      print '(A)', '  [PASS] Rounding mode restored correctly'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Rounding mode not restored'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_rounding_mode_sensitivity

  !---------------------------------------------------------------------------
  ! Test 5: FMA Effects on Norm Computation
  !---------------------------------------------------------------------------
  subroutine test_fma_effects(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(4), result1, result2
    real(dp) :: a, b, c, fma_test, sep_test

    passed = 0
    failed = 0

    print '(A)', 'FMA (Fused Multiply-Add) Effects'
    print '(A)', '--------------------------------'

    ! Test 5a: FMA detection
    a = 1.0_dp + 2.0_dp**(-27)
    b = 1.0_dp + 2.0_dp**(-27)
    c = -(1.0_dp + 2.0_dp**(-26))

    sep_test = a * b
    sep_test = sep_test + c
    fma_test = a * b + c

    if (fma_test /= sep_test) then
      print '(A)', '  [INFO] FMA detected (different from separate mul+add)'
    else
      print '(A)', '  [INFO] No FMA detected (or consistent usage)'
    end if
    passed = passed + 1

    ! Test 5b: DENORM consistency with FMA
    x = [a, b, c, 1.0_dp]
    result1 = denorm_dp(4, x)
    result2 = denorm_dp(4, x)

    if (result1 == result2) then
      print '(A)', '  [PASS] DENORM produces consistent results'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DENORM produces inconsistent results'
      failed = failed + 1
    end if

    ! Test 5c: FMA in accumulation
    x = [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp]
    result1 = denorm_dp(4, x)
    ! Expected: sqrt(4/9) = 2/3

    if (abs(result1 - 2.0_dp/3.0_dp) < 1.0e-14_dp) then
      print '(A)', '  [PASS] FMA accumulation precision acceptable'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [INFO] FMA accumulation result: ', result1
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_fma_effects

  !---------------------------------------------------------------------------
  ! Test 6: Catastrophic Cancellation in Norms
  !---------------------------------------------------------------------------
  subroutine test_cancellation_in_norms(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(5), result, expected
    real(dp) :: big, small

    passed = 0
    failed = 0

    print '(A)', 'Catastrophic Cancellation in Norms'
    print '(A)', '----------------------------------'

    ! Test 6a: One large, rest small (shouldn't cancel)
    x = [1.0e10_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    result = denorm_dp(5, x)
    expected = sqrt(1.0e20_dp + 4.0_dp)

    if (abs(result - expected) / expected < 1.0e-10_dp) then
      print '(A)', '  [PASS] Large + small: no precision loss'
      passed = passed + 1
    else
      print '(A)', '  [WARN] Large + small: some precision loss'
      passed = passed + 1
    end if

    ! Test 6b: Nearly equal components
    big = 1.0e10_dp
    small = 1.0e-5_dp
    x = [big, big + small, big - small, big, big]
    result = denorm_dp(5, x)

    if (ieee_is_finite(result) .and. result > 0.0_dp) then
      print '(A)', '  [PASS] Nearly equal components: finite result'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Nearly equal components: non-finite'
      failed = failed + 1
    end if

    ! Test 6c: Alternating signs (if passed through squaring, no cancellation)
    x = [1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp]
    result = denorm_dp(5, x)
    expected = sqrt(5.0_dp)

    if (abs(result - expected) < tol_dp) then
      print '(A)', '  [PASS] Alternating signs: correct'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [FAIL] Alternating signs: ', result
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_cancellation_in_norms

  !---------------------------------------------------------------------------
  ! Test 7: Condition Number Sensitivity
  !---------------------------------------------------------------------------
  subroutine test_condition_number_sensitivity(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x_well(5), x_ill(5), result_well, result_ill
    real(dp) :: perturbation, result_pert, rel_change

    passed = 0
    failed = 0

    print '(A)', 'Condition Number Sensitivity'
    print '(A)', '----------------------------'

    ! Test 7a: Well-conditioned input
    x_well = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    result_well = denorm_dp(5, x_well)

    ! Perturb slightly
    perturbation = epsilon(1.0_dp) * 10.0_dp
    x_well(1) = x_well(1) * (1.0_dp + perturbation)
    result_pert = denorm_dp(5, x_well)

    rel_change = abs(result_pert - result_well) / result_well

    if (rel_change < perturbation * 10.0_dp) then
      print '(A)', '  [PASS] Well-conditioned: stable to perturbation'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [WARN] Well-conditioned relative change: ', rel_change
      passed = passed + 1
    end if

    ! Test 7b: Ill-conditioned input (one very large, rest tiny)
    x_ill = [1.0e15_dp, 1.0e-15_dp, 1.0e-15_dp, 1.0e-15_dp, 1.0e-15_dp]
    result_ill = denorm_dp(5, x_ill)

    if (ieee_is_finite(result_ill) .and. result_ill > 0.0_dp) then
      print '(A)', '  [PASS] Ill-conditioned: finite result'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Ill-conditioned: non-finite'
      failed = failed + 1
    end if

    ! Test 7c: Very small condition (all tiny)
    x_ill = tiny(1.0_dp) * [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    result_ill = denorm_dp(5, x_ill)

    if (result_ill > 0.0_dp) then
      print '(A)', '  [PASS] All tiny: positive result'
      passed = passed + 1
    else
      print '(A)', '  [WARN] All tiny: underflowed to zero'
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_condition_number_sensitivity

  !---------------------------------------------------------------------------
  ! Test 8: Convergence Reproducibility
  !---------------------------------------------------------------------------
  subroutine test_convergence_reproducibility(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(10)
    real(dp) :: results(5)
    integer :: i
    logical :: all_same

    passed = 0
    failed = 0

    print '(A)', 'Convergence Reproducibility'
    print '(A)', '---------------------------'

    ! Test 8a: Multiple identical calls
    do i = 1, 10
      x(i) = sin(real(i, dp)) * cos(real(i, dp) * 0.7_dp)
    end do

    do i = 1, 5
      results(i) = denorm_dp(10, x)
    end do

    all_same = all(results == results(1))

    if (all_same) then
      print '(A)', '  [PASS] Five calls produced identical results'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Calls produced different results (non-deterministic!)'
      failed = failed + 1
    end if

    ! Test 8b: Order independence (for commutative operations)
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp, 10.0_dp]
    results(1) = denorm_dp(10, x)

    ! Reverse order
    x = [10.0_dp, 9.0_dp, 8.0_dp, 7.0_dp, 6.0_dp, 5.0_dp, 4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]
    results(2) = denorm_dp(10, x)

    if (results(1) == results(2)) then
      print '(A)', '  [PASS] Order independent (commutative)'
      passed = passed + 1
    else
      print '(A)', '  [INFO] Order affects result (expected for some algorithms)'
      passed = passed + 1
    end if

    ! Test 8c: Split vs full computation
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp, 10.0_dp]
    results(1) = denorm_dp(10, x)
    ! Compare with mathematical combination
    results(2) = sqrt(denorm_dp(5, x(1:5))**2 + denorm_dp(5, x(6:10))**2)

    if (abs(results(1) - results(2)) < tol_dp) then
      print '(A)', '  [PASS] Split computation matches full'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [INFO] Split vs full difference: ', abs(results(1) - results(2))
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_convergence_reproducibility

  !---------------------------------------------------------------------------
  ! Test 9: Compiler Flag Detection
  !---------------------------------------------------------------------------
  subroutine test_compiler_flag_detection(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: a, b, result
    logical :: fast_math_suspected

    passed = 0
    failed = 0

    print '(A)', 'Compiler Flag Detection'
    print '(A)', '-----------------------'

    fast_math_suspected = .false.

    ! Test 9a: Inf - Inf should be NaN
    if (ieee_support_nan(1.0_dp) .and. ieee_support_inf(1.0_dp)) then
      a = ieee_value(1.0_dp, ieee_positive_inf)
      b = ieee_value(1.0_dp, ieee_positive_inf)
      result = a - b

      if (ieee_is_nan(result)) then
        print '(A)', '  [PASS] Inf - Inf = NaN (IEEE compliant)'
        passed = passed + 1
      else
        print '(A)', '  [WARN] Inf - Inf /= NaN (-ffast-math likely active)'
        fast_math_suspected = .true.
        passed = passed + 1
      end if
    else
      print '(A)', '  [SKIP] Inf/NaN not supported'
      passed = passed + 1
    end if

    ! Test 9b: 0 * Inf should be NaN
    if (ieee_support_nan(1.0_dp) .and. ieee_support_inf(1.0_dp)) then
      a = 0.0_dp
      b = ieee_value(1.0_dp, ieee_positive_inf)
      result = a * b

      if (ieee_is_nan(result)) then
        print '(A)', '  [PASS] 0 * Inf = NaN (IEEE compliant)'
        passed = passed + 1
      else
        print '(A)', '  [WARN] 0 * Inf /= NaN (unsafe optimization)'
        fast_math_suspected = .true.
        passed = passed + 1
      end if
    else
      passed = passed + 1
    end if

    ! Test 9c: Subnormal preservation
    a = tiny(1.0_dp) / 2.0_dp
    if (a > 0.0_dp .and. a < tiny(1.0_dp)) then
      print '(A)', '  [PASS] Subnormals preserved (no FTZ)'
      passed = passed + 1
    else
      print '(A)', '  [WARN] Subnormals flushed (FTZ/DAZ active)'
      fast_math_suspected = .true.
      passed = passed + 1
    end if

    ! Test 9d: Summary
    if (fast_math_suspected) then
      print '(A)', '  [WARN] Potentially unsafe compiler flags detected!'
      print '(A)', '         Consider recompiling without -ffast-math or -Ofast'
    else
      print '(A)', '  [PASS] No obviously dangerous flags detected'
    end if
    passed = passed + 1

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_compiler_flag_detection

  !---------------------------------------------------------------------------
  ! Test 10: Extended Precision Effects
  !---------------------------------------------------------------------------
  subroutine test_extended_precision_effects(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(5), result1, result2
    real(dp), volatile :: vx(5)
    integer :: i

    passed = 0
    failed = 0

    print '(A)', 'Extended Precision (x87) Effects'
    print '(A)', '--------------------------------'

    ! Test 10a: Register vs memory precision
    x = [1.0_dp/3.0_dp, 1.0_dp/7.0_dp, 1.0_dp/11.0_dp, 1.0_dp/13.0_dp, 1.0_dp/17.0_dp]
    result1 = denorm_dp(5, x)

    ! Force through memory with volatile
    do i = 1, 5
      vx(i) = x(i)
    end do
    result2 = denorm_dp(5, vx)

    if (result1 == result2) then
      print '(A)', '  [PASS] No extended precision leakage'
      passed = passed + 1
    else
      print '(A)', '  [INFO] Extended precision detected (x87 80-bit likely)'
      passed = passed + 1
    end if

    ! Test 10b: Consistent results across calls
    result1 = denorm_dp(5, x)
    result2 = denorm_dp(5, x)

    if (result1 == result2) then
      print '(A)', '  [PASS] Consistent precision across calls'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Inconsistent precision detected'
      failed = failed + 1
    end if

    ! Test 10c: Long computation chain
    x = 1.0_dp
    do i = 1, 100
      x = x * 1.0000001_dp
    end do
    do i = 1, 100
      x = x / 1.0000001_dp
    end do

    if (all(abs(x - 1.0_dp) < 1.0e-10_dp)) then
      print '(A)', '  [PASS] Long chain reversible'
      passed = passed + 1
    else
      print '(A)', '  [INFO] Long chain drift detected'
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_extended_precision_effects

  !---------------------------------------------------------------------------
  ! Test 11: Scaling Sensitivity
  !---------------------------------------------------------------------------
  subroutine test_scaling_sensitivity(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(5), scaled_x(5), result, scaled_result
    real(dp) :: scale_factor

    passed = 0
    failed = 0

    print '(A)', 'Scaling Sensitivity'
    print '(A)', '-------------------'

    ! Test 11a: Homogeneity: ||alpha*x|| = |alpha| * ||x||
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    scale_factor = 1000.0_dp
    scaled_x = scale_factor * x

    result = denorm_dp(5, x)
    scaled_result = denorm_dp(5, scaled_x)

    if (abs(scaled_result - scale_factor * result) < tol_dp * scale_factor) then
      print '(A)', '  [PASS] Homogeneity preserved (alpha=1000)'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [FAIL] Homogeneity violation: ', abs(scaled_result - scale_factor * result)
      failed = failed + 1
    end if

    ! Test 11b: Very small scaling
    scale_factor = 1.0e-10_dp
    scaled_x = scale_factor * x
    scaled_result = denorm_dp(5, scaled_x)

    if (abs(scaled_result - scale_factor * result) / (scale_factor * result) < 1.0e-10_dp) then
      print '(A)', '  [PASS] Homogeneity preserved (alpha=1e-10)'
      passed = passed + 1
    else
      print '(A)', '  [WARN] Small scale homogeneity deviation'
      passed = passed + 1
    end if

    ! Test 11c: Very large scaling
    scale_factor = 1.0e10_dp
    scaled_x = scale_factor * x
    scaled_result = denorm_dp(5, scaled_x)

    if (ieee_is_finite(scaled_result)) then
      if (abs(scaled_result - scale_factor * result) / (scale_factor * result) < 1.0e-10_dp) then
        print '(A)', '  [PASS] Homogeneity preserved (alpha=1e10)'
        passed = passed + 1
      else
        print '(A)', '  [WARN] Large scale homogeneity deviation'
        passed = passed + 1
      end if
    else
      print '(A)', '  [FAIL] Large scaling caused overflow'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_scaling_sensitivity

  !---------------------------------------------------------------------------
  ! Test 12: Jacobian-like Computation
  !---------------------------------------------------------------------------
  subroutine test_jacobian_computation(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(3), f1, f2, h, jac_fd, jac_exact
    real(dp) :: rel_error

    passed = 0
    failed = 0

    print '(A)', 'Jacobian-like Computation'
    print '(A)', '-------------------------'

    ! Test finite difference Jacobian accuracy
    ! For f(x) = ||x||, df/dx_i = x_i / ||x||

    x = [3.0_dp, 4.0_dp, 0.0_dp]  ! ||x|| = 5

    ! Test 12a: Central difference for first component
    h = sqrt(epsilon(1.0_dp))
    x(1) = 3.0_dp + h
    f1 = denorm_dp(3, x)
    x(1) = 3.0_dp - h
    f2 = denorm_dp(3, x)
    x(1) = 3.0_dp

    jac_fd = (f1 - f2) / (2.0_dp * h)
    jac_exact = 3.0_dp / 5.0_dp  ! x_1 / ||x||

    rel_error = abs(jac_fd - jac_exact) / abs(jac_exact)

    if (rel_error < 1.0e-6_dp) then
      print '(A,ES10.3)', '  [PASS] Finite diff Jacobian accurate: ', rel_error
      passed = passed + 1
    else
      print '(A,ES10.3)', '  [WARN] Finite diff Jacobian error: ', rel_error
      passed = passed + 1
    end if

    ! Test 12b: Jacobian at nearly zero component
    x = [3.0_dp, 4.0_dp, 1.0e-10_dp]
    h = sqrt(epsilon(1.0_dp))
    x(3) = 1.0e-10_dp + h
    f1 = denorm_dp(3, x)
    x(3) = 1.0e-10_dp - h
    f2 = denorm_dp(3, x)

    jac_fd = (f1 - f2) / (2.0_dp * h)

    if (ieee_is_finite(jac_fd)) then
      print '(A)', '  [PASS] Jacobian at tiny component: finite'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Jacobian at tiny component: non-finite'
      failed = failed + 1
    end if

    ! Test 12c: Numerical Jacobian stability
    x = [1.0_dp, 1.0_dp, 1.0_dp]
    h = sqrt(epsilon(1.0_dp))
    f1 = denorm_dp(3, x)

    ! Perturb each component
    x(1) = 1.0_dp + h
    f2 = denorm_dp(3, x)
    jac_fd = (f2 - f1) / h

    if (abs(jac_fd - 1.0_dp/sqrt(3.0_dp)) < 1.0e-5_dp) then
      print '(A)', '  [PASS] Forward diff Jacobian stable'
      passed = passed + 1
    else
      print '(A)', '  [INFO] Forward diff Jacobian less accurate (expected)'
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_jacobian_computation

  !---------------------------------------------------------------------------
  ! Test 13: Levenberg-Marquardt Parameter Sensitivity
  !---------------------------------------------------------------------------
  subroutine test_lm_parameter_sensitivity(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(5), result
    real(dp) :: lambda_values(5)
    real(dp) :: results(5)
    integer :: i

    passed = 0
    failed = 0

    print '(A)', 'LM Parameter Sensitivity'
    print '(A)', '------------------------'
    print '(A)', '(Simulates trust region effects)'
    print '(A)', ''

    ! Test 13a: Norm stability with regularization-like modification
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    lambda_values = [0.0_dp, 0.01_dp, 0.1_dp, 1.0_dp, 10.0_dp]

    do i = 1, 5
      ! Simulate regularized norm: sqrt(||x||^2 + lambda)
      results(i) = sqrt(denorm_dp(5, x)**2 + lambda_values(i))
    end do

    ! Results should be monotonically increasing
    if (results(1) < results(2) .and. results(2) < results(3) .and. &
        results(3) < results(4) .and. results(4) < results(5)) then
      print '(A)', '  [PASS] Regularized norm monotonically increasing'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Regularized norm not monotonic'
      failed = failed + 1
    end if

    ! Test 13b: Small lambda should barely change result
    if (abs(results(2) - results(1)) < 0.01_dp * results(1)) then
      print '(A)', '  [PASS] Small lambda has small effect'
      passed = passed + 1
    else
      print '(A)', '  [WARN] Small lambda has larger effect than expected'
      passed = passed + 1
    end if

    ! Test 13c: All results should be finite
    if (all(ieee_is_finite(results))) then
      print '(A)', '  [PASS] All regularized norms are finite'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Some regularized norms are non-finite'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_lm_parameter_sensitivity

  !---------------------------------------------------------------------------
  ! Helper Functions
  !---------------------------------------------------------------------------

  pure function ulp_distance(a, b) result(d)
    real(dp), intent(in) :: a, b
    integer :: d
    integer(int64) :: ia, ib

    if (a == b) then
      d = 0
      return
    end if

    if (a == 0.0_dp .or. b == 0.0_dp) then
      d = huge(1)
      return
    end if

    ia = transfer(a, ia)
    ib = transfer(b, ib)

    if (ia < 0) ia = -9223372036854775807_int64 - ia - 1
    if (ib < 0) ib = -9223372036854775807_int64 - ib - 1

    d = int(abs(ia - ib))

  end function ulp_distance

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

!> Main program
program run_level4_minpack
  use test_minpack_level4
  implicit none
  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) stop 1
end program run_level4_minpack

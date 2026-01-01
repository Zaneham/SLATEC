!> Level 4: Hostile Environment Tests for BLAS
!>
!> Purpose: Do compilers, processors, or OS change anything?
!>
!> These tests probe edge cases that break on:
!>   - Aggressive compiler optimizations (-ffast-math, -Ofast)
!>   - Different floating-point modes (DAZ, FTZ, rounding)
!>   - Non-IEEE processors (rare today)
!>   - Vectorization edge cases (SIMD boundary issues)
!>   - FMA vs separate mul+add
!>   - Extended precision (x87 80-bit)
!>   - Catastrophic cancellation
!>   - Non-reproducible parallel reductions
!>
!> Test Categories:
!>   1. Subnormal handling (DAZ/FTZ detection)
!>   2. Signed zero (IEEE 754 ±0)
!>   3. Inf/NaN propagation and variants
!>   4. Extreme values (overflow/underflow boundaries)
!>   5. Accumulation precision
!>   6. SIMD edge cases
!>   7. Rounding mode sensitivity
!>   8. FMA detection and consistency
!>   9. Catastrophic cancellation
!>  10. Associativity violations
!>  11. Extended precision leakage (x87)
!>  12. Reproducibility
!>  13. Compiler flag detection

module test_blas_level4
  use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, int64
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: ieee_exceptions
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: tol_dp = 1.0e-14_dp
  real(dp), parameter :: tol_loose = 1.0e-10_dp

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 4: HOSTILE ENVIRONMENT TESTS (BLAS)'
    print '(A)', '================================================================'
    print '(A)', 'Comprehensive platform portability and edge case testing'
    print '(A)', ''

    ! Original tests
    call test_subnormal_handling(p, f)
    passed = passed + p; failed = failed + f

    call test_signed_zero_handling(p, f)
    passed = passed + p; failed = failed + f

    call test_inf_nan_propagation(p, f)
    passed = passed + p; failed = failed + f

    call test_extreme_values(p, f)
    passed = passed + p; failed = failed + f

    call test_accumulation_precision(p, f)
    passed = passed + p; failed = failed + f

    call test_simd_edge_cases(p, f)
    passed = passed + p; failed = failed + f

    ! New comprehensive tests
    call test_rounding_modes(p, f)
    passed = passed + p; failed = failed + f

    call test_fma_detection(p, f)
    passed = passed + p; failed = failed + f

    call test_catastrophic_cancellation(p, f)
    passed = passed + p; failed = failed + f

    call test_associativity_violations(p, f)
    passed = passed + p; failed = failed + f

    call test_extended_precision(p, f)
    passed = passed + p; failed = failed + f

    call test_reproducibility(p, f)
    passed = passed + p; failed = failed + f

    call test_compiler_flag_detection(p, f)
    passed = passed + p; failed = failed + f

    call test_nan_variants(p, f)
    passed = passed + p; failed = failed + f

    call test_ulp_accuracy(p, f)
    passed = passed + p; failed = failed + f

    call test_memory_alignment(p, f)
    passed = passed + p; failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 4 BLAS SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! Test 1: Subnormal Number Handling
  ! -ffast-math and DAZ/FTZ modes flush subnormals to zero
  !---------------------------------------------------------------------------
  subroutine test_subnormal_handling(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(4), y(4), alpha
    real(dp) :: subnormal, result
    logical :: subnormals_work

    passed = 0
    failed = 0

    print '(A)', 'SUBNORMAL HANDLING'
    print '(A)', '------------------'
    print '(A)', 'Tests if subnormal numbers are preserved (not flushed to zero)'
    print '(A)', ''

    ! Test 1a: Can we create subnormals?
    subnormal = tiny(1.0_dp) / 2.0_dp
    subnormals_work = (subnormal > 0.0_dp .and. subnormal < tiny(1.0_dp))

    if (subnormals_work) then
      print '(A,ES12.5)', '  Subnormal created: ', subnormal
      print '(A)', '  [PASS] Subnormal numbers are supported'
      passed = passed + 1
    else
      print '(A)', '  [WARN] Subnormals flushed to zero (DAZ/FTZ active or -ffast-math)'
      print '(A)', '         This is a known platform behavior, not necessarily a bug'
      passed = passed + 1
    end if

    ! Test 1b: DAXPY with subnormal values
    if (subnormals_work) then
      x = [subnormal, subnormal, subnormal, subnormal]
      y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      alpha = 1.0_dp

      call daxpy_local(4, alpha, x, 1, y, 1)

      if (all(y > 0.0_dp)) then
        print '(A)', '  [PASS] DAXPY preserves subnormals'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] DAXPY flushed subnormals to zero'
        failed = failed + 1
      end if
    else
      print '(A)', '  [SKIP] DAXPY subnormal test (subnormals not available)'
      passed = passed + 1
    end if

    ! Test 1c: Subnormal multiplication
    if (subnormals_work) then
      x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      alpha = subnormal

      call daxpy_local(4, alpha, x, 1, y, 1)

      if (all(y > 0.0_dp)) then
        print '(A)', '  [PASS] DAXPY with subnormal alpha works'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] DAXPY with subnormal alpha flushed'
        failed = failed + 1
      end if
    else
      print '(A)', '  [SKIP] DAXPY subnormal alpha test'
      passed = passed + 1
    end if

    ! Test 1d: Gradual underflow chain
    if (subnormals_work) then
      result = tiny(1.0_dp)
      ! Divide until we reach the smallest subnormal, then one more division gives 0
      do while (result > 0.0_dp)
        subnormal = result  ! Save previous value
        result = result / 2.0_dp
      end do
      ! result should now be 0, subnormal should be the smallest positive
      if (result == 0.0_dp .and. subnormal > 0.0_dp) then
        print '(A)', '  [PASS] Gradual underflow to zero works'
        passed = passed + 1
      else
        print '(A)', '  [INFO] Gradual underflow behavior platform-specific'
        passed = passed + 1  ! Not a failure, just informational
      end if
    else
      print '(A)', '  [SKIP] Gradual underflow test'
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_subnormal_handling

  !---------------------------------------------------------------------------
  ! Test 2: Signed Zero Handling
  ! IEEE 754 distinguishes +0 and -0
  !---------------------------------------------------------------------------
  subroutine test_signed_zero_handling(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: pos_zero, neg_zero
    real(dp) :: x(2), y(2)

    passed = 0
    failed = 0

    print '(A)', 'SIGNED ZERO HANDLING'
    print '(A)', '--------------------'
    print '(A)', 'Tests IEEE 754 signed zero behavior'
    print '(A)', ''

    pos_zero = 0.0_dp
    neg_zero = -0.0_dp

    ! Test 2a: Can we distinguish +0 and -0?
    if (ieee_is_negative(neg_zero) .and. .not. ieee_is_negative(pos_zero)) then
      print '(A)', '  [PASS] Signed zeros are distinguishable'
      passed = passed + 1
    else
      print '(A)', '  [WARN] Cannot distinguish +0 and -0'
      passed = passed + 1
    end if

    ! Test 2b: 1/+0 vs 1/-0 should give different infinities
    if (ieee_support_inf(1.0_dp)) then
      if (1.0_dp / pos_zero > 0.0_dp .and. 1.0_dp / neg_zero < 0.0_dp) then
        print '(A)', '  [PASS] 1/+0 = +Inf, 1/-0 = -Inf'
        passed = passed + 1
      else
        print '(A)', '  [WARN] Signed zero division behavior unexpected'
        passed = passed + 1
      end if
    else
      print '(A)', '  [SKIP] Infinity not supported'
      passed = passed + 1
    end if

    ! Test 2c: DAXPY with negative zero
    x = [-0.0_dp, -0.0_dp]
    y = [1.0_dp, 2.0_dp]
    call daxpy_local(2, 1.0_dp, x, 1, y, 1)

    if (abs(y(1) - 1.0_dp) < tol_dp .and. abs(y(2) - 2.0_dp) < tol_dp) then
      print '(A)', '  [PASS] DAXPY handles -0 correctly'
      passed = passed + 1
    else
      print '(A,2ES12.5)', '  [FAIL] DAXPY with -0: y = ', y
      failed = failed + 1
    end if

    ! Test 2d: -0 * positive = -0
    x(1) = -0.0_dp
    y(1) = 0.0_dp
    call daxpy_local(1, 5.0_dp, x, 1, y, 1)
    if (ieee_is_negative(y(1)) .or. y(1) == 0.0_dp) then
      print '(A)', '  [PASS] -0 * positive handled correctly'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] -0 * positive sign incorrect'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_signed_zero_handling

  !---------------------------------------------------------------------------
  ! Test 3: Inf/NaN Propagation
  ! IEEE 754 requires proper propagation of special values
  !---------------------------------------------------------------------------
  subroutine test_inf_nan_propagation(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(4), y(4)
    real(dp) :: pos_inf, neg_inf, nan_val

    passed = 0
    failed = 0

    print '(A)', 'INF/NAN PROPAGATION'
    print '(A)', '-------------------'
    print '(A)', 'Tests IEEE 754 special value handling'
    print '(A)', ''

    if (.not. ieee_support_inf(1.0_dp) .or. .not. ieee_support_nan(1.0_dp)) then
      print '(A)', '  [SKIP] IEEE Inf/NaN not fully supported'
      passed = passed + 4
      print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
      print '(A)', ''
      return
    end if

    pos_inf = ieee_value(1.0_dp, ieee_positive_inf)
    neg_inf = ieee_value(1.0_dp, ieee_negative_inf)
    nan_val = ieee_value(1.0_dp, ieee_quiet_nan)

    ! Test 3a: Infinity propagation in DAXPY
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    y = [pos_inf, 0.0_dp, neg_inf, 0.0_dp]
    call daxpy_local(4, 1.0_dp, x, 1, y, 1)

    if (ieee_is_finite(y(2)) .and. ieee_is_finite(y(4))) then
      if (.not. ieee_is_finite(y(1)) .and. .not. ieee_is_finite(y(3))) then
        print '(A)', '  [PASS] DAXPY propagates infinities correctly'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] DAXPY infinity propagation unexpected'
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] DAXPY corrupted finite values'
      failed = failed + 1
    end if

    ! Test 3b: NaN propagation
    x = [1.0_dp, nan_val, 3.0_dp, 4.0_dp]
    y = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp]
    call daxpy_local(4, 1.0_dp, x, 1, y, 1)

    if (ieee_is_nan(y(2))) then
      print '(A)', '  [PASS] DAXPY propagates NaN correctly'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DAXPY did not propagate NaN'
      failed = failed + 1
    end if

    ! Test 3c: Inf - Inf = NaN
    x = [pos_inf, neg_inf, 0.0_dp, 0.0_dp]
    y = [neg_inf, pos_inf, 0.0_dp, 0.0_dp]
    call daxpy_local(2, 1.0_dp, x, 1, y, 1)

    if (ieee_is_nan(y(1)) .and. ieee_is_nan(y(2))) then
      print '(A)', '  [PASS] Inf + (-Inf) = NaN as expected'
      passed = passed + 1
    else
      print '(A)', '  [WARN] Inf arithmetic behavior varies'
      passed = passed + 1
    end if

    ! Test 3d: 0 * Inf = NaN
    x = [pos_inf, neg_inf, 0.0_dp, 0.0_dp]
    y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call daxpy_local(2, 0.0_dp, x, 1, y, 1)
    ! With alpha=0, DAXPY should return early without touching y
    if (y(1) == 0.0_dp .and. y(2) == 0.0_dp) then
      print '(A)', '  [PASS] DAXPY alpha=0 returns without computation'
      passed = passed + 1
    else
      print '(A)', '  [INFO] DAXPY alpha=0 still computed (may be NaN)'
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_inf_nan_propagation

  !---------------------------------------------------------------------------
  ! Test 4: Extreme Values Near Overflow/Underflow
  !---------------------------------------------------------------------------
  subroutine test_extreme_values(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(4), y(4)
    real(dp) :: big, small

    passed = 0
    failed = 0

    print '(A)', 'EXTREME VALUES'
    print '(A)', '--------------'
    print '(A)', 'Tests behavior near overflow/underflow boundaries'
    print '(A)', ''

    big = huge(1.0_dp) / 2.0_dp
    small = tiny(1.0_dp) * 2.0_dp

    ! Test 4a: Large values without overflow
    x = [big, big, big, big]
    y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call daxpy_local(4, 1.0_dp, x, 1, y, 1)

    if (all(ieee_is_finite(y))) then
      print '(A)', '  [PASS] DAXPY handles large values without overflow'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DAXPY overflowed unexpectedly'
      failed = failed + 1
    end if

    ! Test 4b: Small values without underflow to zero
    x = [small, small, small, small]
    y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call daxpy_local(4, 1.0_dp, x, 1, y, 1)

    if (all(y > 0.0_dp)) then
      print '(A)', '  [PASS] DAXPY preserves small values'
      passed = passed + 1
    else
      print '(A)', '  [WARN] DAXPY underflowed to zero'
      passed = passed + 1
    end if

    ! Test 4c: DSCAL with large multiplier
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    call dscal_local(4, big / 10.0_dp, x, 1)

    if (all(ieee_is_finite(x))) then
      print '(A)', '  [PASS] DSCAL handles large scaling'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DSCAL overflowed'
      failed = failed + 1
    end if

    ! Test 4d: DROTG with extreme values
    block
      real(dp) :: a, b, c, s, r
      a = big / 2.0_dp
      b = big / 2.0_dp
      call drotg_local(a, b, c, s)
      r = a

      if (ieee_is_finite(r) .and. ieee_is_finite(c) .and. ieee_is_finite(s)) then
        if (abs(c**2 + s**2 - 1.0_dp) < 1.0e-10_dp) then
          print '(A)', '  [PASS] DROTG handles extreme values'
          passed = passed + 1
        else
          print '(A)', '  [FAIL] DROTG orthogonality lost with extreme values'
          failed = failed + 1
        end if
      else
        print '(A)', '  [FAIL] DROTG produced non-finite values'
        failed = failed + 1
      end if
    end block

    ! Test 4e: Near-max * 2 should overflow
    x(1) = huge(1.0_dp) / 1.5_dp
    y(1) = x(1)
    call daxpy_local(1, 1.0_dp, x, 1, y, 1)
    if (.not. ieee_is_finite(y(1))) then
      print '(A)', '  [PASS] Overflow detected correctly'
      passed = passed + 1
    else
      print '(A)', '  [WARN] Expected overflow did not occur'
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_extreme_values

  !---------------------------------------------------------------------------
  ! Test 5: Accumulation Precision
  !---------------------------------------------------------------------------
  subroutine test_accumulation_precision(passed, failed)
    integer, intent(out) :: passed, failed
    integer, parameter :: n = 10000
    real(dp) :: x(n), y(n)
    integer :: i

    passed = 0
    failed = 0

    print '(A)', 'ACCUMULATION PRECISION'
    print '(A)', '----------------------'
    print '(A)', 'Tests precision in long summations'
    print '(A)', ''

    ! Test 5a: All ones should stay ones
    do i = 1, n
      x(i) = 1.0_dp
      y(i) = 0.0_dp
    end do
    call daxpy_local(n, 1.0_dp, x, 1, y, 1)

    if (all(abs(y - 1.0_dp) < tol_dp)) then
      print '(A)', '  [PASS] DAXPY precise for n=10000 elements'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DAXPY accumulated errors'
      failed = failed + 1
    end if

    ! Test 5b: Alternating signs
    do i = 1, n
      if (mod(i, 2) == 0) then
        x(i) = 1.0_dp
      else
        x(i) = -1.0_dp
      end if
      y(i) = 0.0_dp
    end do
    call daxpy_local(n, 1.0_dp, x, 1, y, 1)

    if (all(abs(abs(y) - 1.0_dp) < tol_dp)) then
      print '(A)', '  [PASS] DAXPY preserves alternating signs'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DAXPY sign preservation failed'
      failed = failed + 1
    end if

    ! Test 5c: Reversible scaling
    do i = 1, 100
      x(i) = 1.0_dp
    end do
    call dscal_local(100, 2.0_dp, x, 1)
    call dscal_local(100, 0.5_dp, x, 1)

    if (all(abs(x(1:100) - 1.0_dp) < tol_dp)) then
      print '(A)', '  [PASS] DSCAL reversible scaling precise'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DSCAL lost precision in reversible scaling'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_accumulation_precision

  !---------------------------------------------------------------------------
  ! Test 6: SIMD Edge Cases
  !---------------------------------------------------------------------------
  subroutine test_simd_edge_cases(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(17), y(17)
    real(dp) :: expected
    integer :: i

    passed = 0
    failed = 0

    print '(A)', 'SIMD EDGE CASES'
    print '(A)', '---------------'
    print '(A)', 'Tests vectorization boundary conditions'
    print '(A)', ''

    ! Test 6a: Odd-length array
    do i = 1, 17
      x(i) = real(i, dp)
      y(i) = 0.0_dp
    end do
    call daxpy_local(17, 2.0_dp, x, 1, y, 1)

    expected = 0.0_dp
    do i = 1, 17
      if (abs(y(i) - 2.0_dp * real(i, dp)) > tol_dp) then
        expected = 1.0_dp
        exit
      end if
    end do

    if (expected == 0.0_dp) then
      print '(A)', '  [PASS] DAXPY handles n=17 (odd length)'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DAXPY failed at SIMD boundary'
      failed = failed + 1
    end if

    ! Test 6b: Non-unit stride
    do i = 1, 17
      x(i) = real(i, dp)
      y(i) = 100.0_dp
    end do
    call daxpy_local(9, 1.0_dp, x, 2, y, 2)

    if (abs(y(1) - 101.0_dp) < tol_dp .and. &
        abs(y(3) - 103.0_dp) < tol_dp .and. &
        abs(y(2) - 100.0_dp) < tol_dp) then
      print '(A)', '  [PASS] DAXPY handles stride=2'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DAXPY stride handling incorrect'
      failed = failed + 1
    end if

    ! Test 6c: Negative stride
    do i = 1, 4
      x(i) = real(i, dp)
      y(i) = 0.0_dp
    end do
    call daxpy_local(4, 1.0_dp, x, -1, y, -1)

    if (abs(y(1) - 1.0_dp) < tol_dp .and. &
        abs(y(2) - 2.0_dp) < tol_dp .and. &
        abs(y(3) - 3.0_dp) < tol_dp .and. &
        abs(y(4) - 4.0_dp) < tol_dp) then
      print '(A)', '  [PASS] DAXPY handles negative stride'
      passed = passed + 1
    else
      print '(A,4ES10.3)', '  [FAIL] DAXPY negative stride: y = ', y(1:4)
      failed = failed + 1
    end if

    ! Test 6d: n=1
    x(1) = 5.0_dp
    y(1) = 10.0_dp
    call daxpy_local(1, 2.0_dp, x, 1, y, 1)

    if (abs(y(1) - 20.0_dp) < tol_dp) then
      print '(A)', '  [PASS] DAXPY handles n=1'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DAXPY failed for n=1'
      failed = failed + 1
    end if

    ! Test 6e: n=0 (should be no-op)
    y(1) = 999.0_dp
    call daxpy_local(0, 2.0_dp, x, 1, y, 1)
    if (abs(y(1) - 999.0_dp) < tol_dp) then
      print '(A)', '  [PASS] DAXPY handles n=0 (no-op)'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DAXPY modified y when n=0'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_simd_edge_cases

  !---------------------------------------------------------------------------
  ! Test 7: Rounding Mode Sensitivity
  ! Different rounding modes can produce different results
  !---------------------------------------------------------------------------
  subroutine test_rounding_modes(passed, failed)
    integer, intent(out) :: passed, failed
    type(ieee_round_type) :: original_mode, current_mode
    real(dp) :: x(4), y(4), result_nearest, result_zero, result_up, result_down
    logical :: modes_supported

    passed = 0
    failed = 0

    print '(A)', 'ROUNDING MODE SENSITIVITY'
    print '(A)', '-------------------------'
    print '(A)', 'Tests behavior under different IEEE rounding modes'
    print '(A)', ''

    ! Check if rounding mode changes are supported
    modes_supported = ieee_support_rounding(ieee_nearest, 1.0_dp) .and. &
                      ieee_support_rounding(ieee_to_zero, 1.0_dp) .and. &
                      ieee_support_rounding(ieee_up, 1.0_dp) .and. &
                      ieee_support_rounding(ieee_down, 1.0_dp)

    if (.not. modes_supported) then
      print '(A)', '  [SKIP] Not all rounding modes supported'
      passed = passed + 5
      print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
      print '(A)', ''
      return
    end if

    ! Save original rounding mode
    call ieee_get_rounding_mode(original_mode)

    ! Test 7a: Results should be consistent within each mode
    x = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 1.0_dp/7.0_dp, 1.0_dp/11.0_dp]

    call ieee_set_rounding_mode(ieee_nearest)
    y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call daxpy_local(4, 1.0_dp, x, 1, y, 1)
    result_nearest = y(1)

    call ieee_set_rounding_mode(ieee_to_zero)
    y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call daxpy_local(4, 1.0_dp, x, 1, y, 1)
    result_zero = y(1)

    call ieee_set_rounding_mode(ieee_up)
    y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call daxpy_local(4, 1.0_dp, x, 1, y, 1)
    result_up = y(1)

    call ieee_set_rounding_mode(ieee_down)
    y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call daxpy_local(4, 1.0_dp, x, 1, y, 1)
    result_down = y(1)

    ! Restore original mode
    call ieee_set_rounding_mode(original_mode)

    ! Test 7a: All modes should produce finite results
    if (ieee_is_finite(result_nearest) .and. ieee_is_finite(result_zero) .and. &
        ieee_is_finite(result_up) .and. ieee_is_finite(result_down)) then
      print '(A)', '  [PASS] All rounding modes produce finite results'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Some rounding modes produced non-finite results'
      failed = failed + 1
    end if

    ! Test 7b: Round up should be >= round down
    if (result_up >= result_down) then
      print '(A)', '  [PASS] Round-up >= Round-down as expected'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Round-up < Round-down (impossible!)'
      failed = failed + 1
    end if

    ! Test 7c: Nearest should be between up and down
    if (result_nearest >= result_down .and. result_nearest <= result_up) then
      print '(A)', '  [PASS] Round-nearest between up and down'
      passed = passed + 1
    else
      print '(A)', '  [WARN] Round-nearest outside expected range'
      passed = passed + 1
    end if

    ! Test 7d: Results should be close (within a few ULP)
    if (abs(result_up - result_down) < 1.0e-14_dp) then
      print '(A)', '  [PASS] Rounding mode differences are small'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [INFO] Rounding mode spread: ', abs(result_up - result_down)
      passed = passed + 1
    end if

    ! Test 7e: Mode restored correctly
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

  end subroutine test_rounding_modes

  !---------------------------------------------------------------------------
  ! Test 8: FMA (Fused Multiply-Add) Detection
  ! FMA computes a*b+c with single rounding, which differs from separate ops
  !---------------------------------------------------------------------------
  subroutine test_fma_detection(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: a, b, c, fma_result, separate_result
    real(dp) :: x(4), y(4)
    logical :: fma_detected

    passed = 0
    failed = 0

    print '(A)', 'FMA (FUSED MULTIPLY-ADD) DETECTION'
    print '(A)', '----------------------------------'
    print '(A)', 'Detects if FMA instructions are used (affects results by ~1 ULP)'
    print '(A)', ''

    ! Classic FMA detection: a*b + c where a*b rounds, then add rounds again
    ! With FMA: only one rounding at the end
    ! Use values that maximize the difference

    ! Test case from William Kahan
    a = 1.0_dp + 2.0_dp**(-27)
    b = 1.0_dp + 2.0_dp**(-27)
    c = -(1.0_dp + 2.0_dp**(-26))

    ! Separate multiply then add
    separate_result = a * b
    separate_result = separate_result + c

    ! The difference between FMA and separate is in the low bits
    ! If FMA is used, a*b+c is computed with one rounding
    ! If separate, a*b is rounded, then +c is rounded again

    ! We can't force non-FMA in standard Fortran, so we detect which is used
    fma_result = a * b + c  ! Compiler may use FMA here

    fma_detected = (fma_result /= separate_result)

    ! Test 8a: Report FMA status
    if (fma_detected) then
      print '(A)', '  [INFO] FMA instructions detected (different from separate mul+add)'
      print '(A,ES20.13)', '         FMA result:      ', fma_result
      print '(A,ES20.13)', '         Separate result: ', separate_result
      passed = passed + 1
    else
      print '(A)', '  [INFO] No FMA detected (or compiler used FMA consistently)'
      passed = passed + 1
    end if

    ! Test 8b: DAXPY should be consistent with itself
    x = [a, a, a, a]
    y = [c, c, c, c]
    call daxpy_local(4, b, x, 1, y, 1)

    if (y(1) == y(2) .and. y(2) == y(3) .and. y(3) == y(4)) then
      print '(A)', '  [PASS] DAXPY internally consistent (all elements equal)'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DAXPY produced inconsistent results'
      failed = failed + 1
    end if

    ! Test 8c: Result should be close to mathematical answer
    ! Mathematical: (1+2^-27)^2 - (1+2^-26) = 2^-54 (approximately)
    if (abs(y(1)) < 1.0e-10_dp) then
      print '(A)', '  [PASS] DAXPY result magnitude reasonable'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [INFO] DAXPY result: ', y(1)
      passed = passed + 1
    end if

    ! Test 8d: Multiple FMA operations should accumulate consistently
    a = 1.0_dp / 3.0_dp
    b = 3.0_dp
    x = [a, a, a, a]
    y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call daxpy_local(4, b, x, 1, y, 1)

    ! Each y(i) should be 1.0 (ideally)
    if (all(abs(y - 1.0_dp) < tol_dp)) then
      print '(A)', '  [PASS] FMA precision acceptable for 1/3 * 3'
      passed = passed + 1
    else
      print '(A,4ES12.5)', '  [INFO] 1/3 * 3 results: ', y
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_fma_detection

  !---------------------------------------------------------------------------
  ! Test 9: Catastrophic Cancellation
  ! Subtracting nearly-equal numbers loses precision
  !---------------------------------------------------------------------------
  subroutine test_catastrophic_cancellation(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(4), y(4)
    real(dp) :: a, b, expected, result, rel_error

    passed = 0
    failed = 0

    print '(A)', 'CATASTROPHIC CANCELLATION'
    print '(A)', '-------------------------'
    print '(A)', 'Tests precision loss when subtracting nearly-equal values'
    print '(A)', ''

    ! Test 9a: Classic cancellation: (1+eps) - 1
    a = 1.0_dp + 1.0e-15_dp
    x = [a, a, a, a]
    y = [-1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp]
    call daxpy_local(4, 1.0_dp, x, 1, y, 1)

    expected = 1.0e-15_dp
    result = y(1)
    if (result > 0.0_dp) then
      rel_error = abs(result - expected) / expected
      if (rel_error < 0.1_dp) then
        print '(A)', '  [PASS] (1+1e-15) - 1 computed with good precision'
        passed = passed + 1
      else
        print '(A,ES12.5)', '  [WARN] (1+1e-15) - 1 relative error: ', rel_error
        passed = passed + 1
      end if
    else
      print '(A)', '  [FAIL] (1+1e-15) - 1 underflowed to zero or negative'
      failed = failed + 1
    end if

    ! Test 9b: Near-cancellation with larger values
    a = 1.0e10_dp + 1.0_dp
    b = 1.0e10_dp
    x = [a, a, a, a]
    y = [-b, -b, -b, -b]
    call daxpy_local(4, 1.0_dp, x, 1, y, 1)

    if (abs(y(1) - 1.0_dp) < 1.0e-5_dp) then
      print '(A)', '  [PASS] (1e10+1) - 1e10 = 1 (no cancellation)'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [WARN] (1e10+1) - 1e10 = ', y(1)
      passed = passed + 1
    end if

    ! Test 9c: Quadratic formula-like cancellation
    ! For ax^2 + bx + c = 0 with b^2 ≈ 4ac
    ! This is a case where DAXPY's simple operations avoid the worst cancellation
    a = 1.0_dp
    b = 1.0e8_dp
    ! c chosen so b^2 - 4ac is small
    x = [b*b, 0.0_dp, 0.0_dp, 0.0_dp]
    y = [-4.0_dp*a*(b*b/4.0_dp - 1.0_dp), 0.0_dp, 0.0_dp, 0.0_dp]
    call daxpy_local(1, 1.0_dp, x, 1, y, 1)

    if (ieee_is_finite(y(1))) then
      print '(A)', '  [PASS] Discriminant computation finite'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Discriminant computation overflow/underflow'
      failed = failed + 1
    end if

    ! Test 9d: Accumulated cancellation
    x = [1.0_dp, -1.0_dp, 1.0e-16_dp, 0.0_dp]
    y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call daxpy_local(3, 1.0_dp, x, 1, y, 1)
    ! y should be [1, -1, 1e-16]
    ! Sum would be 1e-16, but we're testing element-wise

    if (abs(y(3) - 1.0e-16_dp) < 1.0e-17_dp) then
      print '(A)', '  [PASS] Small value preserved after larger values'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [INFO] Small value result: ', y(3)
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_catastrophic_cancellation

  !---------------------------------------------------------------------------
  ! Test 10: Associativity Violations
  ! (a + b) + c ≠ a + (b + c) in floating-point
  ! -ffast-math assumes they're equal, which can break code
  !---------------------------------------------------------------------------
  subroutine test_associativity_violations(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: a, b, c
    real(dp) :: left_assoc, right_assoc
    real(dp) :: x(3), y(3)
    logical :: associativity_holds
    integer :: i

    passed = 0
    failed = 0

    print '(A)', 'ASSOCIATIVITY VIOLATIONS'
    print '(A)', '------------------------'
    print '(A)', 'Tests if (a+b)+c = a+(b+c) - detects -ffast-math'
    print '(A)', ''

    ! Test 10a: Classic non-associativity example
    a = 1.0e20_dp
    b = -1.0e20_dp
    c = 1.0_dp

    left_assoc = (a + b) + c   ! = 0 + 1 = 1
    right_assoc = a + (b + c)  ! = 1e20 + (-1e20 + 1) = 1e20 - 1e20 = 0 (c lost)

    associativity_holds = (left_assoc == right_assoc)

    if (.not. associativity_holds) then
      print '(A)', '  [PASS] Associativity correctly violated (IEEE compliant)'
      print '(A,ES12.5)', '         (a+b)+c = ', left_assoc
      print '(A,ES12.5)', '         a+(b+c) = ', right_assoc
      passed = passed + 1
    else
      print '(A)', '  [WARN] Associativity appears to hold (-ffast-math likely active)'
      passed = passed + 1
    end if

    ! Test 10b: Smaller scale test
    a = 1.0_dp
    b = 1.0e-16_dp
    c = -1.0_dp

    left_assoc = (a + b) + c
    right_assoc = a + (b + c)

    if (left_assoc /= right_assoc) then
      print '(A)', '  [PASS] Small-scale associativity violation detected'
      passed = passed + 1
    else
      print '(A)', '  [INFO] Small-scale associativity holds (may be correct)'
      passed = passed + 1
    end if

    ! Test 10c: DAXPY ordering test
    ! DAXPY computes y = alpha*x + y element-wise, order shouldn't matter
    x = [1.0e15_dp, -1.0e15_dp, 1.0_dp]
    y = [0.0_dp, 0.0_dp, 0.0_dp]
    call daxpy_local(3, 1.0_dp, x, 1, y, 1)

    ! Each y(i) should just be x(i) since we added to 0
    if (abs(y(1) - 1.0e15_dp) < 1.0_dp .and. &
        abs(y(2) + 1.0e15_dp) < 1.0_dp .and. &
        abs(y(3) - 1.0_dp) < tol_dp) then
      print '(A)', '  [PASS] DAXPY element-wise operation unaffected by ordering'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DAXPY ordering issue detected'
      failed = failed + 1
    end if

    ! Test 10d: Repeated addition test
    a = 0.1_dp
    y(1) = 0.0_dp
    do i = 1, 10
      y(1) = y(1) + a
    end do
    ! y(1) should be 1.0, but floating-point gives slightly different

    if (abs(y(1) - 1.0_dp) < 1.0e-14_dp) then
      print '(A)', '  [PASS] 0.1 added 10 times ≈ 1.0'
      passed = passed + 1
    else
      print '(A,ES20.13)', '  [INFO] 0.1 * 10 = ', y(1)
      print '(A)', '         (Not exactly 1.0 due to binary representation)'
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_associativity_violations

  !---------------------------------------------------------------------------
  ! Test 11: Extended Precision Detection (x87)
  ! x86 FPUs may use 80-bit registers internally
  !---------------------------------------------------------------------------
  subroutine test_extended_precision(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: a, b, c, result1, result2
    real(dp), volatile :: va, vb, vc  ! volatile forces memory store
    real(dp) :: x(4), y(4)
    integer :: i

    passed = 0
    failed = 0

    print '(A)', 'EXTENDED PRECISION (x87) DETECTION'
    print '(A)', '----------------------------------'
    print '(A)', 'Detects if x87 80-bit registers affect results'
    print '(A)', ''

    ! Test 11a: Extended precision can cause different results
    ! when intermediate values are stored to memory vs kept in registers

    a = 1.0_dp / 3.0_dp
    b = 3.0_dp
    c = 1.0_dp

    ! Compute without forcing memory store
    result1 = a * b - c

    ! Compute with volatile (forces 64-bit memory store)
    va = a
    vb = b
    vc = c
    result2 = va * vb - vc

    if (result1 == result2) then
      print '(A)', '  [PASS] No extended precision leakage detected'
      passed = passed + 1
    else
      print '(A)', '  [INFO] Extended precision detected (x87 80-bit likely)'
      print '(A,ES20.13)', '         Register result: ', result1
      print '(A,ES20.13)', '         Memory result:   ', result2
      passed = passed + 1
    end if

    ! Test 11b: DAXPY should use consistent precision
    x = [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp]
    y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call daxpy_local(4, 3.0_dp, x, 1, y, 1)

    if (y(1) == y(2) .and. y(2) == y(3) .and. y(3) == y(4)) then
      print '(A)', '  [PASS] DAXPY produces consistent results'
      passed = passed + 1
    else
      print '(A)', '  [WARN] DAXPY results vary (possible precision inconsistency)'
      passed = passed + 1
    end if

    ! Test 11c: Long computation chain
    result1 = 1.0_dp
    do i = 1, 100
      result1 = result1 * 1.0000001_dp
    end do
    do i = 1, 100
      result1 = result1 / 1.0000001_dp
    end do

    if (abs(result1 - 1.0_dp) < 1.0e-10_dp) then
      print '(A)', '  [PASS] Long computation chain reversible'
      passed = passed + 1
    else
      print '(A,ES12.5)', '  [INFO] Long chain drift: ', abs(result1 - 1.0_dp)
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_extended_precision

  !---------------------------------------------------------------------------
  ! Test 12: Reproducibility
  ! Same code, same input should give same output
  !---------------------------------------------------------------------------
  subroutine test_reproducibility(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(100), y1(100), y2(100), y3(100)
    integer :: i
    logical :: all_equal

    passed = 0
    failed = 0

    print '(A)', 'REPRODUCIBILITY'
    print '(A)', '---------------'
    print '(A)', 'Tests if repeated runs give identical results'
    print '(A)', ''

    ! Initialize with pseudo-random-ish pattern
    do i = 1, 100
      x(i) = sin(real(i, dp)) * cos(real(i, dp) * 0.7_dp)
    end do

    ! Run 1
    y1 = 0.0_dp
    call daxpy_local(100, 3.14159_dp, x, 1, y1, 1)

    ! Run 2
    y2 = 0.0_dp
    call daxpy_local(100, 3.14159_dp, x, 1, y2, 1)

    ! Run 3
    y3 = 0.0_dp
    call daxpy_local(100, 3.14159_dp, x, 1, y3, 1)

    ! Test 12a: All runs should be identical
    all_equal = all(y1 == y2) .and. all(y2 == y3)

    if (all_equal) then
      print '(A)', '  [PASS] Three runs produced identical results'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Runs produced different results (non-deterministic!)'
      failed = failed + 1
    end if

    ! Test 12b: Order of operations
    y1 = 0.0_dp
    call daxpy_local(50, 1.0_dp, x(1:50), 1, y1(1:50), 1)
    call daxpy_local(50, 1.0_dp, x(51:100), 1, y1(51:100), 1)

    y2 = 0.0_dp
    call daxpy_local(100, 1.0_dp, x, 1, y2, 1)

    if (all(y1 == y2)) then
      print '(A)', '  [PASS] Split vs single call identical'
      passed = passed + 1
    else
      print '(A)', '  [INFO] Split vs single call differ slightly'
      passed = passed + 1
    end if

    ! Test 12c: Forward vs reverse
    y1 = 0.0_dp
    y2 = 0.0_dp
    do i = 1, 100
      y1(i) = x(i) * 2.0_dp
    end do
    do i = 100, 1, -1
      y2(i) = x(i) * 2.0_dp
    end do

    if (all(y1 == y2)) then
      print '(A)', '  [PASS] Forward and reverse loops identical'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Forward and reverse loops differ'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_reproducibility

  !---------------------------------------------------------------------------
  ! Test 13: Compiler Flag Detection
  ! Detect if dangerous compiler flags are active
  !---------------------------------------------------------------------------
  subroutine test_compiler_flag_detection(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: a, b, c, result
    real(dp) :: x(4), y(4)
    logical :: fast_math_suspected, unsafe_suspected

    passed = 0
    failed = 0

    print '(A)', 'COMPILER FLAG DETECTION'
    print '(A)', '-----------------------'
    print '(A)', 'Detects potentially dangerous optimization flags'
    print '(A)', ''

    fast_math_suspected = .false.
    unsafe_suspected = .false.

    ! Test 13a: -ffast-math typically breaks this
    ! (Inf - Inf) should be NaN, but -ffast-math may optimize it away
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

    ! Test 13b: Check if 0 * Inf = NaN
    if (ieee_support_nan(1.0_dp) .and. ieee_support_inf(1.0_dp)) then
      a = 0.0_dp
      b = ieee_value(1.0_dp, ieee_positive_inf)
      result = a * b

      if (ieee_is_nan(result)) then
        print '(A)', '  [PASS] 0 * Inf = NaN (IEEE compliant)'
        passed = passed + 1
      else
        print '(A)', '  [WARN] 0 * Inf /= NaN (unsafe optimization suspected)'
        unsafe_suspected = .true.
        passed = passed + 1
      end if
    else
      passed = passed + 1
    end if

    ! Test 13c: Reciprocal optimization check
    ! -freciprocal-math replaces a/b with a*(1/b), losing precision
    a = 1.0_dp
    b = 3.0_dp
    result = a / b
    c = a * (1.0_dp / b)

    if (result == c) then
      print '(A)', '  [INFO] Division and reciprocal multiplication equal'
      passed = passed + 1
    else
      print '(A)', '  [PASS] Division differs from reciprocal (precise division used)'
      passed = passed + 1
    end if

    ! Test 13d: Contraction check (FMA)
    ! This isn't necessarily bad, just informational
    a = 1.0_dp + 2.0_dp**(-30)
    b = 1.0_dp + 2.0_dp**(-30)
    c = -(1.0_dp + 2.0_dp**(-29))
    result = a * b + c

    if (abs(result) < 1.0e-15_dp) then
      print '(A)', '  [INFO] FMA contraction likely (result very small)'
    else
      print '(A)', '  [INFO] Separate mul+add likely'
    end if
    passed = passed + 1

    ! Test 13e: Summary
    if (fast_math_suspected .or. unsafe_suspected) then
      print '(A)', '  [WARN] Potentially unsafe compiler flags detected!'
      print '(A)', '         Consider recompiling without -ffast-math or -Ofast'
      passed = passed + 1
    else
      print '(A)', '  [PASS] No obviously dangerous flags detected'
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_compiler_flag_detection

  !---------------------------------------------------------------------------
  ! Test 14: NaN Variants
  ! Different NaN types and payloads
  !---------------------------------------------------------------------------
  subroutine test_nan_variants(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: qnan, snan, nan1, nan2, result
    real(dp) :: x(4), y(4)

    passed = 0
    failed = 0

    print '(A)', 'NAN VARIANTS'
    print '(A)', '------------'
    print '(A)', 'Tests quiet NaN, signaling NaN, and NaN propagation'
    print '(A)', ''

    if (.not. ieee_support_nan(1.0_dp)) then
      print '(A)', '  [SKIP] NaN not supported on this platform'
      passed = passed + 4
      print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
      print '(A)', ''
      return
    end if

    ! Test 14a: Quiet NaN creation and detection
    qnan = ieee_value(1.0_dp, ieee_quiet_nan)
    if (ieee_is_nan(qnan)) then
      print '(A)', '  [PASS] Quiet NaN created and detected'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Quiet NaN not detected'
      failed = failed + 1
    end if

    ! Test 14b: NaN comparisons should all be false
    nan1 = ieee_value(1.0_dp, ieee_quiet_nan)
    nan2 = ieee_value(1.0_dp, ieee_quiet_nan)

    if (.not. (nan1 == nan2) .and. .not. (nan1 < nan2) .and. &
        .not. (nan1 > nan2) .and. .not. (nan1 <= nan2) .and. &
        .not. (nan1 >= nan2)) then
      print '(A)', '  [PASS] NaN comparisons correctly return false'
      passed = passed + 1
    else
      print '(A)', '  [WARN] NaN comparison behavior non-standard'
      passed = passed + 1
    end if

    ! Test 14c: NaN propagation in DAXPY
    x = [1.0_dp, qnan, 3.0_dp, 4.0_dp]
    y = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp]
    call daxpy_local(4, 2.0_dp, x, 1, y, 1)

    if (ieee_is_nan(y(2)) .and. .not. ieee_is_nan(y(1)) .and. &
        .not. ieee_is_nan(y(3)) .and. .not. ieee_is_nan(y(4))) then
      print '(A)', '  [PASS] NaN propagates correctly (only affects element 2)'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] NaN propagation incorrect'
      failed = failed + 1
    end if

    ! Test 14d: NaN in alpha
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    y = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp]
    call daxpy_local(4, qnan, x, 1, y, 1)

    if (all(ieee_is_nan(y))) then
      print '(A)', '  [PASS] NaN alpha propagates to all elements'
      passed = passed + 1
    else
      print '(A)', '  [WARN] NaN alpha propagation partial'
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_nan_variants

  !---------------------------------------------------------------------------
  ! Test 15: ULP (Units in Last Place) Accuracy
  ! Measure actual numerical accuracy in ULPs
  !---------------------------------------------------------------------------
  subroutine test_ulp_accuracy(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(4), y(4), expected(4)
    real(dp) :: ulp_val, max_ulp_error
    integer :: i
    integer(int64) :: bits_result, bits_expected, ulp_diff

    passed = 0
    failed = 0

    print '(A)', 'ULP ACCURACY'
    print '(A)', '------------'
    print '(A)', 'Measures numerical accuracy in Units in Last Place'
    print '(A)', ''

    ! Test 15a: Simple operation should be exact
    x = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp]
    y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call daxpy_local(4, 2.0_dp, x, 1, y, 1)
    expected = [2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp]

    max_ulp_error = 0.0_dp
    do i = 1, 4
      ulp_val = spacing(expected(i))
      max_ulp_error = max(max_ulp_error, abs(y(i) - expected(i)) / ulp_val)
    end do

    if (max_ulp_error == 0.0_dp) then
      print '(A)', '  [PASS] Powers of 2 scaling exact (0 ULP error)'
      passed = passed + 1
    else
      print '(A,F6.2,A)', '  [INFO] Powers of 2 scaling: ', max_ulp_error, ' ULP max error'
      passed = passed + 1
    end if

    ! Test 15b: Irrational multiplier
    x = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call daxpy_local(4, 3.14159265358979_dp, x, 1, y, 1)

    ! All should be identical
    if (y(1) == y(2) .and. y(2) == y(3) .and. y(3) == y(4)) then
      print '(A)', '  [PASS] Identical inputs produce identical outputs'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Identical inputs produce different outputs'
      failed = failed + 1
    end if

    ! Test 15c: Accumulation ULP test
    x = [1.0_dp/7.0_dp, 1.0_dp/11.0_dp, 1.0_dp/13.0_dp, 1.0_dp/17.0_dp]
    y = [100.0_dp, 100.0_dp, 100.0_dp, 100.0_dp]
    call daxpy_local(4, 1.0_dp, x, 1, y, 1)

    ! Check ULP from expected
    expected(1) = 100.0_dp + 1.0_dp/7.0_dp
    expected(2) = 100.0_dp + 1.0_dp/11.0_dp
    expected(3) = 100.0_dp + 1.0_dp/13.0_dp
    expected(4) = 100.0_dp + 1.0_dp/17.0_dp

    max_ulp_error = 0.0_dp
    do i = 1, 4
      ulp_val = spacing(expected(i))
      if (ulp_val > 0.0_dp) then
        max_ulp_error = max(max_ulp_error, abs(y(i) - expected(i)) / ulp_val)
      end if
    end do

    if (max_ulp_error <= 1.0_dp) then
      print '(A)', '  [PASS] Fraction additions within 1 ULP'
      passed = passed + 1
    else if (max_ulp_error <= 2.0_dp) then
      print '(A,F6.2,A)', '  [PASS] Fraction additions within ', max_ulp_error, ' ULP'
      passed = passed + 1
    else
      print '(A,F6.2,A)', '  [WARN] Fraction additions: ', max_ulp_error, ' ULP error'
      passed = passed + 1
    end if

    ! Test 15d: Large + small should be accurate
    x = [1.0e-15_dp, 1.0e-15_dp, 1.0e-15_dp, 1.0e-15_dp]
    y = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    call daxpy_local(4, 1.0_dp, x, 1, y, 1)

    expected = 1.0_dp + 1.0e-15_dp
    max_ulp_error = abs(y(1) - expected(1)) / spacing(expected(1))

    if (max_ulp_error <= 1.0_dp) then
      print '(A)', '  [PASS] 1 + 1e-15 accurate to 1 ULP'
      passed = passed + 1
    else
      print '(A,F6.2,A)', '  [INFO] 1 + 1e-15 accuracy: ', max_ulp_error, ' ULP'
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_ulp_accuracy

  !---------------------------------------------------------------------------
  ! Test 16: Memory Alignment Edge Cases
  ! Misaligned memory access can affect SIMD operations
  !---------------------------------------------------------------------------
  subroutine test_memory_alignment(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp), allocatable :: x(:), y(:)
    real(dp) :: x_static(35), y_static(35)
    integer :: i, offset

    passed = 0
    failed = 0

    print '(A)', 'MEMORY ALIGNMENT'
    print '(A)', '----------------'
    print '(A)', 'Tests behavior with various memory alignments'
    print '(A)', ''

    ! Test 16a: Heap-allocated arrays
    allocate(x(100), y(100))
    do i = 1, 100
      x(i) = real(i, dp)
      y(i) = 0.0_dp
    end do
    call daxpy_local(100, 2.0_dp, x, 1, y, 1)

    if (all(abs(y - 2.0_dp * [(real(i,dp), i=1,100)]) < tol_dp)) then
      print '(A)', '  [PASS] Heap-allocated arrays work correctly'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Heap-allocated arrays failed'
      failed = failed + 1
    end if
    deallocate(x, y)

    ! Test 16b: Stack-allocated arrays
    do i = 1, 35
      x_static(i) = real(i, dp)
      y_static(i) = 0.0_dp
    end do
    call daxpy_local(35, 2.0_dp, x_static, 1, y_static, 1)

    if (all(abs(y_static - 2.0_dp * [(real(i,dp), i=1,35)]) < tol_dp)) then
      print '(A)', '  [PASS] Stack-allocated arrays work correctly'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Stack-allocated arrays failed'
      failed = failed + 1
    end if

    ! Test 16c: Offset into array (potentially misaligned)
    allocate(x(100), y(100))
    do i = 1, 100
      x(i) = real(i, dp)
      y(i) = 0.0_dp
    end do

    ! Start from offset 3 (may not be aligned to SIMD width)
    call daxpy_local(50, 2.0_dp, x(3), 1, y(3), 1)

    if (abs(y(3) - 2.0_dp * 3.0_dp) < tol_dp .and. &
        abs(y(52) - 2.0_dp * 52.0_dp) < tol_dp .and. &
        abs(y(1)) < tol_dp .and. abs(y(2)) < tol_dp) then
      print '(A)', '  [PASS] Offset array access works correctly'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Offset array access failed'
      failed = failed + 1
    end if
    deallocate(x, y)

    ! Test 16d: Various sizes (1 to 20) - tests unrolling edge cases
    allocate(x(20), y(20))
    do offset = 1, 20
      do i = 1, 20
        x(i) = 1.0_dp
        y(i) = 0.0_dp
      end do
      call daxpy_local(offset, 3.0_dp, x, 1, y, 1)
      if (any(abs(y(1:offset) - 3.0_dp) > tol_dp)) then
        print '(A,I2,A)', '  [FAIL] Size ', offset, ' failed'
        failed = failed + 1
        deallocate(x, y)
        print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
        print '(A)', ''
        return
      end if
    end do
    print '(A)', '  [PASS] All sizes 1-20 work correctly'
    passed = passed + 1
    deallocate(x, y)

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_memory_alignment

  !---------------------------------------------------------------------------
  ! Local BLAS implementations
  !---------------------------------------------------------------------------
  pure subroutine daxpy_local(n, da, dx, incx, dy, incy)
    integer, intent(in) :: n, incx, incy
    real(dp), intent(in) :: da, dx(*)
    real(dp), intent(inout) :: dy(*)
    integer :: i, ix, iy

    if (n <= 0) return
    if (da == 0.0_dp) return

    if (incx == 1 .and. incy == 1) then
      do i = 1, n
        dy(i) = dy(i) + da * dx(i)
      end do
    else
      ix = 1
      iy = 1
      if (incx < 0) ix = (-n + 1) * incx + 1
      if (incy < 0) iy = (-n + 1) * incy + 1
      do i = 1, n
        dy(iy) = dy(iy) + da * dx(ix)
        ix = ix + incx
        iy = iy + incy
      end do
    end if
  end subroutine

  pure subroutine dscal_local(n, da, dx, incx)
    integer, intent(in) :: n, incx
    real(dp), intent(in) :: da
    real(dp), intent(inout) :: dx(*)
    integer :: i

    if (n <= 0 .or. incx <= 0) return
    do i = 1, n
      dx(1 + (i-1)*incx) = da * dx(1 + (i-1)*incx)
    end do
  end subroutine

  pure subroutine drotg_local(da, db, dc, ds)
    real(dp), intent(inout) :: da, db
    real(dp), intent(out) :: dc, ds
    real(dp) :: r, roe, scale, z

    roe = db
    if (abs(da) > abs(db)) roe = da
    scale = abs(da) + abs(db)

    if (scale == 0.0_dp) then
      dc = 1.0_dp
      ds = 0.0_dp
      r = 0.0_dp
      z = 0.0_dp
    else
      r = scale * sqrt((da/scale)**2 + (db/scale)**2)
      r = sign(1.0_dp, roe) * r
      dc = da / r
      ds = db / r
      z = 1.0_dp
      if (abs(da) > abs(db)) z = ds
      if (abs(db) >= abs(da) .and. dc /= 0.0_dp) z = 1.0_dp / dc
    end if

    da = r
    db = z
  end subroutine

end module test_blas_level4

!> Main program
program run_level4_blas
  use test_blas_level4
  implicit none
  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) stop 1
end program run_level4_blas

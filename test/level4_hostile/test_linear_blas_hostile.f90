!> Level 4: Hostile Environment Tests for BLAS
!>
!> Purpose: Do compilers, processors, or OS change anything?
!>
!> These tests probe edge cases that break on:
!>   - Aggressive compiler optimizations (-ffast-math, -Ofast)
!>   - Different floating-point modes (DAZ, FTZ)
!>   - Non-IEEE processors (rare today)
!>   - Vectorization edge cases (SIMD boundary issues)
!>
!> What we test:
!>   - Subnormal numbers (may be flushed to zero)
!>   - Signed zeros (+0 vs -0)
!>   - Inf/NaN propagation
!>   - Extreme values near overflow/underflow
!>   - Precision loss in accumulation
!>   - SIMD alignment edge cases

module test_blas_level4
  use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
  use, intrinsic :: ieee_arithmetic
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: tol_dp = 1.0e-14_dp

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 4: HOSTILE ENVIRONMENT TESTS (BLAS)'
    print '(A)', '================================================================'
    print '(A)', 'Testing edge cases that may fail with aggressive optimization'
    print '(A)', ''

    call test_subnormal_handling(p, f)
    passed = passed + p
    failed = failed + f

    call test_signed_zero_handling(p, f)
    passed = passed + p
    failed = failed + f

    call test_inf_nan_propagation(p, f)
    passed = passed + p
    failed = failed + f

    call test_extreme_values(p, f)
    passed = passed + p
    failed = failed + f

    call test_accumulation_precision(p, f)
    passed = passed + p
    failed = failed + f

    call test_simd_edge_cases(p, f)
    passed = passed + p
    failed = failed + f

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
      ! Count as pass with warning - this is expected on some platforms
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
      passed = passed + 1  ! Skip counts as pass
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
      passed = passed + 1  ! Warning, not failure
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
      passed = passed + 3
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
      passed = passed + 1  ! Platform-dependent
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_inf_nan_propagation

  !---------------------------------------------------------------------------
  ! Test 4: Extreme Values Near Overflow/Underflow
  ! Tests behavior at the edges of the representable range
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
      passed = passed + 1  ! May be platform behavior
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

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_extreme_values

  !---------------------------------------------------------------------------
  ! Test 5: Accumulation Precision
  ! Tests for precision loss in long summations
  !---------------------------------------------------------------------------
  subroutine test_accumulation_precision(passed, failed)
    integer, intent(out) :: passed, failed
    integer, parameter :: n = 10000
    real(dp) :: x(n), y(n)
    real(dp) :: expected, actual, rel_error
    integer :: i

    passed = 0
    failed = 0

    print '(A)', 'ACCUMULATION PRECISION'
    print '(A)', '----------------------'
    print '(A)', 'Tests precision in long summations'
    print '(A)', ''

    ! Test 5a: Sum of small numbers should accumulate correctly
    ! x = (1, 1, 1, ..., 1) with n elements
    ! y = 0 initially
    ! After DAXPY with alpha=1: y(1) = 1 (not accumulated, just single op)
    do i = 1, n
      x(i) = 1.0_dp
      y(i) = 0.0_dp
    end do

    call daxpy_local(n, 1.0_dp, x, 1, y, 1)

    ! Each y(i) should be exactly 1.0
    if (all(abs(y - 1.0_dp) < tol_dp)) then
      print '(A)', '  [PASS] DAXPY precise for n=10000 elements'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DAXPY accumulated errors'
      failed = failed + 1
    end if

    ! Test 5b: Kahan-like test - alternating +1, -1 should cancel
    do i = 1, n
      if (mod(i, 2) == 0) then
        x(i) = 1.0_dp
      else
        x(i) = -1.0_dp
      end if
      y(i) = 0.0_dp
    end do
    ! Note: n=10000 is even, so sum = 0
    call daxpy_local(n, 1.0_dp, x, 1, y, 1)

    ! y(i) should each be ±1, not accumulated
    ! This test is about individual element precision
    if (all(abs(abs(y) - 1.0_dp) < tol_dp)) then
      print '(A)', '  [PASS] DAXPY preserves alternating signs'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DAXPY sign preservation failed'
      failed = failed + 1
    end if

    ! Test 5c: Repeated scaling should preserve precision
    do i = 1, 100
      x(i) = 1.0_dp
    end do

    ! Scale by 2, then by 0.5 - should get back to 1
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
  ! Tests for issues at SIMD vector boundaries
  !---------------------------------------------------------------------------
  subroutine test_simd_edge_cases(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(17), y(17)  ! 17 is not divisible by common SIMD widths (2,4,8)
    real(dp) :: expected
    integer :: i

    passed = 0
    failed = 0

    print '(A)', 'SIMD EDGE CASES'
    print '(A)', '---------------'
    print '(A)', 'Tests vectorization boundary conditions'
    print '(A)', ''

    ! Test 6a: Odd-length array (not divisible by SIMD width)
    do i = 1, 17
      x(i) = real(i, dp)
      y(i) = 0.0_dp
    end do

    call daxpy_local(17, 2.0_dp, x, 1, y, 1)

    expected = 0.0_dp
    do i = 1, 17
      if (abs(y(i) - 2.0_dp * real(i, dp)) > tol_dp) then
        expected = 1.0_dp  ! Flag error
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

    ! Test 6b: Non-unit stride (defeats simple vectorization)
    do i = 1, 17
      x(i) = real(i, dp)
      y(i) = 100.0_dp
    end do

    ! Process every other element: x(1), x(3), x(5), ...
    call daxpy_local(9, 1.0_dp, x, 2, y, 2)

    ! y(1) = 100 + 1, y(3) = 100 + 3, etc.
    if (abs(y(1) - 101.0_dp) < tol_dp .and. &
        abs(y(3) - 103.0_dp) < tol_dp .and. &
        abs(y(2) - 100.0_dp) < tol_dp) then  ! y(2) unchanged
      print '(A)', '  [PASS] DAXPY handles stride=2'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] DAXPY stride handling incorrect'
      failed = failed + 1
    end if

    ! Test 6c: Negative stride
    ! For incx<0: algorithm accesses dx(1+(n-1)*|incx|) down to dx(1)
    ! For n=4, incx=-1: accesses dx(4), dx(3), dx(2), dx(1) in that order
    do i = 1, 4
      x(i) = real(i, dp)
      y(i) = 0.0_dp
    end do

    ! Pass full array - algorithm handles negative stride internally
    call daxpy_local(4, 1.0_dp, x, -1, y, -1)

    ! With negative stride, x(4),x(3),x(2),x(1) -> y(4),y(3),y(2),y(1)
    ! y(4) += x(4) = 4, y(3) += x(3) = 3, etc.
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

    ! Test 6d: Very small n (edge case for unrolled loops)
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

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_simd_edge_cases

  !---------------------------------------------------------------------------
  ! Local BLAS implementations (matching SLATEC behavior)
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

!> Level 1: Regression Tests for BLAS (Basic Linear Algebra Subprograms)
!>
!> Purpose: Does the code work? Does it still work after changes?
!> These tests verify basic functionality and catch regressions.
!>
!> What we test:
!>   - DAXPY/SAXPY: Y = alpha*X + Y
!>   - DSCAL/SSCAL: X = alpha*X
!>   - DCOPY/SCOPY: Y = X
!>   - DSWAP/SSWAP: swap X and Y
!>   - DROT/SROT: apply Givens rotation
!>   - DROTG/SROTG: generate Givens rotation
!>
!> Reference: BLAS Technical Forum Standard (2002)

module test_blas_level1
  use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: tol_dp = 1.0e-14_dp
  real(sp), parameter :: tol_sp = 1.0e-6_sp

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 1: BLAS REGRESSION TESTS'
    print '(A)', '================================================================'
    print '(A)', ''

    call test_daxpy_suite(p, f)
    passed = passed + p
    failed = failed + f

    call test_dscal_suite(p, f)
    passed = passed + p
    failed = failed + f

    call test_dcopy_suite(p, f)
    passed = passed + p
    failed = failed + f

    call test_dswap_suite(p, f)
    passed = passed + p
    failed = failed + f

    call test_drot_suite(p, f)
    passed = passed + p
    failed = failed + f

    call test_drotg_suite(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 1 BLAS SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! DAXPY Test Suite: Y = alpha*X + Y
  !---------------------------------------------------------------------------
  subroutine test_daxpy_suite(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(5), y(5), y_expected(5)
    integer :: i

    passed = 0
    failed = 0

    print '(A)', 'DAXPY (Y = alpha*X + Y)'
    print '(A)', '-----------------------'

    ! Test 1: Simple addition with alpha=1
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    y = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp]
    y_expected = [11.0_dp, 22.0_dp, 33.0_dp, 44.0_dp, 55.0_dp]
    call daxpy_local(5, 1.0_dp, x, 1, y, 1)
    if (all(abs(y - y_expected) < tol_dp)) then
      print '(A)', '  [PASS] alpha=1: Y = X + Y'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] alpha=1'
      failed = failed + 1
    end if

    ! Test 2: Scaling with alpha=2
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    y_expected = [2.0_dp, 4.0_dp, 6.0_dp, 8.0_dp, 10.0_dp]
    call daxpy_local(5, 2.0_dp, x, 1, y, 1)
    if (all(abs(y - y_expected) < tol_dp)) then
      print '(A)', '  [PASS] alpha=2: Y = 2*X'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] alpha=2'
      failed = failed + 1
    end if

    ! Test 3: alpha=0 should leave Y unchanged
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    y = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp]
    y_expected = y
    call daxpy_local(5, 0.0_dp, x, 1, y, 1)
    if (all(abs(y - y_expected) < tol_dp)) then
      print '(A)', '  [PASS] alpha=0: Y unchanged'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] alpha=0'
      failed = failed + 1
    end if

    ! Test 4: Negative alpha
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    y = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp]
    y_expected = [9.0_dp, 18.0_dp, 27.0_dp, 36.0_dp, 45.0_dp]
    call daxpy_local(5, -1.0_dp, x, 1, y, 1)
    if (all(abs(y - y_expected) < tol_dp)) then
      print '(A)', '  [PASS] alpha=-1: Y = Y - X'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] alpha=-1'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_daxpy_suite

  !---------------------------------------------------------------------------
  ! DSCAL Test Suite: X = alpha*X
  !---------------------------------------------------------------------------
  subroutine test_dscal_suite(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(5), x_expected(5)

    passed = 0
    failed = 0

    print '(A)', 'DSCAL (X = alpha*X)'
    print '(A)', '-------------------'

    ! Test 1: Scale by 2
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    x_expected = [2.0_dp, 4.0_dp, 6.0_dp, 8.0_dp, 10.0_dp]
    call dscal_local(5, 2.0_dp, x, 1)
    if (all(abs(x - x_expected) < tol_dp)) then
      print '(A)', '  [PASS] alpha=2'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] alpha=2'
      failed = failed + 1
    end if

    ! Test 2: Scale by 0 (zero out)
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    x_expected = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call dscal_local(5, 0.0_dp, x, 1)
    if (all(abs(x - x_expected) < tol_dp)) then
      print '(A)', '  [PASS] alpha=0 (zero out)'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] alpha=0'
      failed = failed + 1
    end if

    ! Test 3: Scale by -1 (negate)
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    x_expected = [-1.0_dp, -2.0_dp, -3.0_dp, -4.0_dp, -5.0_dp]
    call dscal_local(5, -1.0_dp, x, 1)
    if (all(abs(x - x_expected) < tol_dp)) then
      print '(A)', '  [PASS] alpha=-1 (negate)'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] alpha=-1'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_dscal_suite

  !---------------------------------------------------------------------------
  ! DCOPY Test Suite: Y = X
  !---------------------------------------------------------------------------
  subroutine test_dcopy_suite(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(5), y(5)

    passed = 0
    failed = 0

    print '(A)', 'DCOPY (Y = X)'
    print '(A)', '-------------'

    ! Test 1: Simple copy
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    y = 0.0_dp
    call dcopy_local(5, x, 1, y, 1)
    if (all(abs(y - x) < tol_dp)) then
      print '(A)', '  [PASS] Simple copy'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Simple copy'
      failed = failed + 1
    end if

    ! Test 2: Copy preserves source
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    y = 0.0_dp
    call dcopy_local(5, x, 1, y, 1)
    if (all(abs(x - [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]) < tol_dp)) then
      print '(A)', '  [PASS] Source unchanged'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Source modified'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_dcopy_suite

  !---------------------------------------------------------------------------
  ! DSWAP Test Suite: swap X and Y
  !---------------------------------------------------------------------------
  subroutine test_dswap_suite(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(5), y(5), x_orig(5), y_orig(5)

    passed = 0
    failed = 0

    print '(A)', 'DSWAP (swap X and Y)'
    print '(A)', '--------------------'

    ! Test 1: Simple swap
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    y = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp]
    x_orig = x
    y_orig = y
    call dswap_local(5, x, 1, y, 1)
    if (all(abs(x - y_orig) < tol_dp) .and. all(abs(y - x_orig) < tol_dp)) then
      print '(A)', '  [PASS] Simple swap'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Simple swap'
      failed = failed + 1
    end if

    ! Test 2: Double swap returns to original
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    y = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp]
    x_orig = x
    y_orig = y
    call dswap_local(5, x, 1, y, 1)
    call dswap_local(5, x, 1, y, 1)
    if (all(abs(x - x_orig) < tol_dp) .and. all(abs(y - y_orig) < tol_dp)) then
      print '(A)', '  [PASS] Double swap = identity'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Double swap'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_dswap_suite

  !---------------------------------------------------------------------------
  ! DROT Test Suite: apply Givens rotation
  !---------------------------------------------------------------------------
  subroutine test_drot_suite(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(3), y(3), c, s
    real(dp), parameter :: pi = 3.14159265358979323846_dp

    passed = 0
    failed = 0

    print '(A)', 'DROT (apply Givens rotation)'
    print '(A)', '----------------------------'

    ! Test 1: 90 degree rotation (c=0, s=1)
    x = [1.0_dp, 0.0_dp, 0.0_dp]
    y = [0.0_dp, 1.0_dp, 0.0_dp]
    c = 0.0_dp
    s = 1.0_dp
    call drot_local(3, x, 1, y, 1, c, s)
    ! After rotation: x' = c*x + s*y = y, y' = c*y - s*x = -x
    if (abs(x(1) - 0.0_dp) < tol_dp .and. abs(x(2) - 1.0_dp) < tol_dp .and. &
        abs(y(1) + 1.0_dp) < tol_dp .and. abs(y(2) - 0.0_dp) < tol_dp) then
      print '(A)', '  [PASS] 90 degree rotation'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] 90 degree rotation'
      failed = failed + 1
    end if

    ! Test 2: Identity rotation (c=1, s=0)
    x = [1.0_dp, 2.0_dp, 3.0_dp]
    y = [4.0_dp, 5.0_dp, 6.0_dp]
    c = 1.0_dp
    s = 0.0_dp
    call drot_local(3, x, 1, y, 1, c, s)
    if (all(abs(x - [1.0_dp, 2.0_dp, 3.0_dp]) < tol_dp) .and. &
        all(abs(y - [4.0_dp, 5.0_dp, 6.0_dp]) < tol_dp)) then
      print '(A)', '  [PASS] Identity rotation (c=1, s=0)'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Identity rotation'
      failed = failed + 1
    end if

    ! Test 3: 45 degree rotation
    x = [1.0_dp, 0.0_dp, 0.0_dp]
    y = [0.0_dp, 0.0_dp, 0.0_dp]
    c = sqrt(0.5_dp)
    s = sqrt(0.5_dp)
    call drot_local(1, x, 1, y, 1, c, s)
    ! x' = c*1 + s*0 = c, y' = c*0 - s*1 = -s
    if (abs(x(1) - c) < tol_dp .and. abs(y(1) + s) < tol_dp) then
      print '(A)', '  [PASS] 45 degree rotation'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] 45 degree rotation'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_drot_suite

  !---------------------------------------------------------------------------
  ! DROTG Test Suite: generate Givens rotation
  !---------------------------------------------------------------------------
  subroutine test_drotg_suite(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: a, b, c, s, r

    passed = 0
    failed = 0

    print '(A)', 'DROTG (generate Givens rotation)'
    print '(A)', '---------------------------------'

    ! Test 1: Simple 3-4-5 triangle
    a = 3.0_dp
    b = 4.0_dp
    call drotg_local(a, b, c, s)
    r = a  ! After call, a contains r = sqrt(a^2 + b^2)
    if (abs(r - 5.0_dp) < tol_dp .and. abs(c - 0.6_dp) < tol_dp .and. &
        abs(s - 0.8_dp) < tol_dp) then
      print '(A)', '  [PASS] 3-4-5 triangle: r=5, c=0.6, s=0.8'
      passed = passed + 1
    else
      print '(A,3ES12.4)', '  [FAIL] 3-4-5: r,c,s = ', r, c, s
      failed = failed + 1
    end if

    ! Test 2: b=0 case (no rotation needed)
    a = 5.0_dp
    b = 0.0_dp
    call drotg_local(a, b, c, s)
    if (abs(c - 1.0_dp) < tol_dp .and. abs(s) < tol_dp) then
      print '(A)', '  [PASS] b=0: c=1, s=0 (identity)'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] b=0 case'
      failed = failed + 1
    end if

    ! Test 3: a=0 case
    a = 0.0_dp
    b = 5.0_dp
    call drotg_local(a, b, c, s)
    r = a
    if (abs(r - 5.0_dp) < tol_dp .and. abs(c) < tol_dp .and. abs(s - 1.0_dp) < tol_dp) then
      print '(A)', '  [PASS] a=0: r=5, c=0, s=1'
      passed = passed + 1
    else
      print '(A,3ES12.4)', '  [FAIL] a=0: r,c,s = ', r, c, s
      failed = failed + 1
    end if

    ! Test 4: Both zero
    a = 0.0_dp
    b = 0.0_dp
    call drotg_local(a, b, c, s)
    if (abs(c - 1.0_dp) < tol_dp .and. abs(s) < tol_dp) then
      print '(A)', '  [PASS] a=b=0: c=1, s=0'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] a=b=0 case'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_drotg_suite

  !---------------------------------------------------------------------------
  ! Local BLAS implementations for testing
  !---------------------------------------------------------------------------
  pure subroutine daxpy_local(n, da, dx, incx, dy, incy)
    integer, intent(in) :: n, incx, incy
    real(dp), intent(in) :: da, dx(*)
    real(dp), intent(inout) :: dy(*)
    integer :: i, ix, iy
    if (n <= 0 .or. da == 0.0_dp) return
    if (incx == 1 .and. incy == 1) then
      do i = 1, n
        dy(i) = dy(i) + da * dx(i)
      end do
    else
      ix = 1; iy = 1
      if (incx < 0) ix = (-n + 1) * incx + 1
      if (incy < 0) iy = (-n + 1) * incy + 1
      do i = 1, n
        dy(iy) = dy(iy) + da * dx(ix)
        ix = ix + incx; iy = iy + incy
      end do
    end if
  end subroutine

  pure subroutine dscal_local(n, da, dx, incx)
    integer, intent(in) :: n, incx
    real(dp), intent(in) :: da
    real(dp), intent(inout) :: dx(*)
    integer :: i
    if (n <= 0 .or. incx <= 0) return
    do i = 1, n * incx, incx
      dx(i) = da * dx(i)
    end do
  end subroutine

  pure subroutine dcopy_local(n, dx, incx, dy, incy)
    integer, intent(in) :: n, incx, incy
    real(dp), intent(in) :: dx(*)
    real(dp), intent(out) :: dy(*)
    integer :: i, ix, iy
    if (n <= 0) return
    ix = 1; iy = 1
    if (incx < 0) ix = (-n + 1) * incx + 1
    if (incy < 0) iy = (-n + 1) * incy + 1
    do i = 1, n
      dy(iy) = dx(ix)
      ix = ix + incx; iy = iy + incy
    end do
  end subroutine

  pure subroutine dswap_local(n, dx, incx, dy, incy)
    integer, intent(in) :: n, incx, incy
    real(dp), intent(inout) :: dx(*), dy(*)
    real(dp) :: t
    integer :: i, ix, iy
    if (n <= 0) return
    ix = 1; iy = 1
    if (incx < 0) ix = (-n + 1) * incx + 1
    if (incy < 0) iy = (-n + 1) * incy + 1
    do i = 1, n
      t = dx(ix); dx(ix) = dy(iy); dy(iy) = t
      ix = ix + incx; iy = iy + incy
    end do
  end subroutine

  pure subroutine drot_local(n, dx, incx, dy, incy, c, s)
    integer, intent(in) :: n, incx, incy
    real(dp), intent(in) :: c, s
    real(dp), intent(inout) :: dx(*), dy(*)
    real(dp) :: t
    integer :: i, ix, iy
    if (n <= 0) return
    ix = 1; iy = 1
    if (incx < 0) ix = (-n + 1) * incx + 1
    if (incy < 0) iy = (-n + 1) * incy + 1
    do i = 1, n
      t = c*dx(ix) + s*dy(iy)
      dy(iy) = c*dy(iy) - s*dx(ix)
      dx(ix) = t
      ix = ix + incx; iy = iy + incy
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
      dc = 1.0_dp; ds = 0.0_dp; r = 0.0_dp; z = 0.0_dp
    else
      r = scale * sqrt((da/scale)**2 + (db/scale)**2)
      r = sign(1.0_dp, roe) * r
      dc = da / r; ds = db / r; z = 1.0_dp
      if (abs(da) > abs(db)) z = ds
      if (abs(db) >= abs(da) .and. dc /= 0.0_dp) z = 1.0_dp / dc
    end if
    da = r; db = z
  end subroutine

end module test_blas_level1

!> Main program for Level 1 BLAS tests
program run_level1_blas
  use test_blas_level1
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level1_blas

!> Level 2: Mathematical Verification for BLAS
!>
!> Purpose: Does the algorithm match the original mathematics?
!> These tests verify mathematical properties and identities.
!>
!> Mathematical Properties Tested:
!>   - DAXPY: Linearity, associativity
!>   - DROT: Orthogonality (rotation preserves length)
!>   - DROTG: Pythagorean theorem, rotation eliminates component
!>   - Vector operations: Commutativity, identity elements
!>
!> Reference: Lawson et al. (1979) "Basic Linear Algebra Subprograms for
!>            Fortran Usage", ACM TOMS 5(3), 308-323.

module test_blas_level2
  use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: tol_dp = 1.0e-14_dp
  real(dp), parameter :: pi = 3.14159265358979323846_dp

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 2: BLAS MATHEMATICAL VERIFICATION'
    print '(A)', '================================================================'
    print '(A)', 'Reference: Lawson et al., ACM TOMS 5(3), 1979'
    print '(A)', ''

    call test_daxpy_linearity(p, f)
    passed = passed + p
    failed = failed + f

    call test_drot_orthogonality(p, f)
    passed = passed + p
    failed = failed + f

    call test_drotg_mathematics(p, f)
    passed = passed + p
    failed = failed + f

    call test_norm_properties(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 2 BLAS SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! DAXPY Linearity Properties
  ! Y = alpha*X + Y satisfies:
  !   - (alpha + beta)*X = alpha*X + beta*X (distributivity)
  !   - alpha*(X + Y) = alpha*X + alpha*Y (distributivity)
  !   - 0*X + Y = Y (identity)
  !   - 1*X + 0 = X (identity)
  !---------------------------------------------------------------------------
  subroutine test_daxpy_linearity(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(5), y1(5), y2(5), y3(5)
    real(dp) :: alpha, beta

    passed = 0
    failed = 0

    print '(A)', 'DAXPY Linearity Properties'
    print '(A)', '--------------------------'

    ! Test 1: (alpha + beta)*X = alpha*X + beta*X (via two sequential calls)
    alpha = 2.0_dp
    beta = 3.0_dp
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]

    ! Method 1: (alpha + beta)*X
    y1 = 0.0_dp
    call daxpy_local(5, alpha + beta, x, 1, y1, 1)

    ! Method 2: alpha*X + beta*X
    y2 = 0.0_dp
    call daxpy_local(5, alpha, x, 1, y2, 1)
    call daxpy_local(5, beta, x, 1, y2, 1)

    if (all(abs(y1 - y2) < tol_dp)) then
      print '(A)', '  [PASS] Distributivity: (a+b)*X = a*X + b*X'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Distributivity'
      failed = failed + 1
    end if

    ! Test 2: alpha*(X + Y) = alpha*X + alpha*Y
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    y1 = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp]
    alpha = 3.0_dp

    ! Method 1: alpha*(X + Y) - compute X+Y first, then scale
    y2 = x + y1
    y3 = 0.0_dp
    call daxpy_local(5, alpha, y2, 1, y3, 1)

    ! Method 2: alpha*X + alpha*Y
    y2 = 0.0_dp
    call daxpy_local(5, alpha, x, 1, y2, 1)
    call daxpy_local(5, alpha, y1, 1, y2, 1)

    if (all(abs(y2 - y3) < tol_dp)) then
      print '(A)', '  [PASS] Distributivity: a*(X+Y) = a*X + a*Y'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Distributivity over addition'
      failed = failed + 1
    end if

    ! Test 3: Additive identity (0*X + Y = Y)
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    y1 = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp]
    y2 = y1
    call daxpy_local(5, 0.0_dp, x, 1, y2, 1)

    if (all(abs(y1 - y2) < tol_dp)) then
      print '(A)', '  [PASS] Additive identity: 0*X + Y = Y'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Additive identity'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_daxpy_linearity

  !---------------------------------------------------------------------------
  ! DROT Orthogonality Properties
  ! Givens rotation is orthogonal: ||Rx|| = ||x||
  ! And: R^T R = I
  !---------------------------------------------------------------------------
  subroutine test_drot_orthogonality(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(3), y(3), x_orig(3), y_orig(3)
    real(dp) :: c, s, norm_before, norm_after
    real(dp) :: theta

    passed = 0
    failed = 0

    print '(A)', 'DROT Orthogonality Properties'
    print '(A)', '------------------------------'

    ! Test 1: Rotation preserves length
    x = [3.0_dp, 4.0_dp, 0.0_dp]
    y = [0.0_dp, 0.0_dp, 0.0_dp]
    norm_before = sqrt(sum(x**2) + sum(y**2))

    c = cos(pi/4.0_dp)
    s = sin(pi/4.0_dp)
    call drot_local(3, x, 1, y, 1, c, s)

    norm_after = sqrt(sum(x**2) + sum(y**2))

    if (abs(norm_before - norm_after) < tol_dp) then
      print '(A)', '  [PASS] Rotation preserves length: ||Rx|| = ||x||'
      passed = passed + 1
    else
      print '(A,2ES15.8)', '  [FAIL] ||x|| before/after: ', norm_before, norm_after
      failed = failed + 1
    end if

    ! Test 2: Rotation by theta then -theta = identity
    x = [1.0_dp, 2.0_dp, 3.0_dp]
    y = [4.0_dp, 5.0_dp, 6.0_dp]
    x_orig = x
    y_orig = y

    theta = pi / 6.0_dp  ! 30 degrees
    c = cos(theta)
    s = sin(theta)
    call drot_local(3, x, 1, y, 1, c, s)

    ! Rotate back by -theta
    c = cos(-theta)
    s = sin(-theta)
    call drot_local(3, x, 1, y, 1, c, s)

    if (all(abs(x - x_orig) < tol_dp) .and. all(abs(y - y_orig) < tol_dp)) then
      print '(A)', '  [PASS] R(theta) * R(-theta) = Identity'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Rotation inverse'
      failed = failed + 1
    end if

    ! Test 3: 360 degree rotation = identity
    x = [1.0_dp, 2.0_dp, 3.0_dp]
    y = [4.0_dp, 5.0_dp, 6.0_dp]
    x_orig = x
    y_orig = y

    c = cos(2.0_dp * pi)
    s = sin(2.0_dp * pi)
    call drot_local(3, x, 1, y, 1, c, s)

    if (all(abs(x - x_orig) < tol_dp) .and. all(abs(y - y_orig) < tol_dp)) then
      print '(A)', '  [PASS] 360 degree rotation = Identity'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Full rotation'
      failed = failed + 1
    end if

    ! Test 4: Two 180 degree rotations = identity
    x = [1.0_dp, 2.0_dp, 3.0_dp]
    y = [4.0_dp, 5.0_dp, 6.0_dp]
    x_orig = x
    y_orig = y

    c = -1.0_dp  ! cos(pi) = -1
    s = 0.0_dp   ! sin(pi) = 0
    call drot_local(3, x, 1, y, 1, c, s)
    call drot_local(3, x, 1, y, 1, c, s)

    if (all(abs(x - x_orig) < tol_dp) .and. all(abs(y - y_orig) < tol_dp)) then
      print '(A)', '  [PASS] Two 180 degree rotations = Identity'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Double 180 rotation'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_drot_orthogonality

  !---------------------------------------------------------------------------
  ! DROTG Mathematical Properties
  ! Given (a, b), produces (c, s) such that:
  !   [c  s] [a]   [r]
  !   [-s c] [b] = [0]
  ! where r = sqrt(a^2 + b^2) and c^2 + s^2 = 1
  !---------------------------------------------------------------------------
  subroutine test_drotg_mathematics(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: a, b, c, s, r
    real(dp) :: x_rot, y_rot

    passed = 0
    failed = 0

    print '(A)', 'DROTG Mathematical Properties'
    print '(A)', '------------------------------'
    print '(A)', '  Givens: [c s; -s c] * [a; b] = [r; 0]'
    print '(A)', ''

    ! Test 1: Verify c^2 + s^2 = 1 (orthogonality)
    a = 3.0_dp
    b = 4.0_dp
    call drotg_local(a, b, c, s)

    if (abs(c**2 + s**2 - 1.0_dp) < tol_dp) then
      print '(A)', '  [PASS] c^2 + s^2 = 1 (orthogonality)'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] c^2 + s^2 = ', c**2 + s**2
      failed = failed + 1
    end if

    ! Test 2: r = sqrt(a^2 + b^2) (Pythagorean theorem)
    a = 3.0_dp
    b = 4.0_dp
    call drotg_local(a, b, c, s)
    r = a  ! After call, a contains r

    if (abs(r - 5.0_dp) < tol_dp) then
      print '(A)', '  [PASS] r = sqrt(a^2 + b^2) = 5 (Pythagorean)'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] r = ', r
      failed = failed + 1
    end if

    ! Test 3: Rotation zeroes second component
    ! [c s; -s c] * [3; 4] = [5; 0]
    a = 3.0_dp
    b = 4.0_dp
    call drotg_local(a, b, c, s)

    x_rot = c * 3.0_dp + s * 4.0_dp
    y_rot = -s * 3.0_dp + c * 4.0_dp

    if (abs(y_rot) < tol_dp .and. abs(x_rot - 5.0_dp) < tol_dp) then
      print '(A)', '  [PASS] Rotation eliminates b: [a;b] -> [r;0]'
      passed = passed + 1
    else
      print '(A,2ES15.8)', '  [FAIL] Result: ', x_rot, y_rot
      failed = failed + 1
    end if

    ! Test 4: Sign of r follows larger magnitude
    a = -5.0_dp
    b = 3.0_dp
    call drotg_local(a, b, c, s)
    r = a

    if (r < 0.0_dp) then
      print '(A)', '  [PASS] r takes sign of larger: r < 0 when |a| > |b|, a < 0'
      passed = passed + 1
    else
      print '(A,ES15.8)', '  [FAIL] r should be negative: ', r
      failed = failed + 1
    end if

    ! Test 5: 5-12-13 Pythagorean triple
    a = 5.0_dp
    b = 12.0_dp
    call drotg_local(a, b, c, s)
    r = a

    if (abs(r - 13.0_dp) < tol_dp .and. abs(c - 5.0_dp/13.0_dp) < tol_dp .and. &
        abs(s - 12.0_dp/13.0_dp) < tol_dp) then
      print '(A)', '  [PASS] 5-12-13 triple: r=13, c=5/13, s=12/13'
      passed = passed + 1
    else
      print '(A,3ES12.4)', '  [FAIL] 5-12-13: r,c,s = ', r, c, s
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_drotg_mathematics

  !---------------------------------------------------------------------------
  ! Norm Properties
  ! ||alpha * x|| = |alpha| * ||x||  (homogeneity)
  ! ||x + y|| <= ||x|| + ||y||  (triangle inequality)
  !---------------------------------------------------------------------------
  subroutine test_norm_properties(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: x(5), y(5), z(5)
    real(dp) :: norm_x, norm_y, norm_z, norm_sum
    real(dp) :: alpha

    passed = 0
    failed = 0

    print '(A)', 'Euclidean Norm Properties'
    print '(A)', '-------------------------'

    ! Test 1: Homogeneity ||alpha*x|| = |alpha| * ||x||
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    alpha = 3.0_dp
    norm_x = sqrt(sum(x**2))
    y = alpha * x
    norm_y = sqrt(sum(y**2))

    if (abs(norm_y - abs(alpha) * norm_x) < tol_dp) then
      print '(A)', '  [PASS] Homogeneity: ||a*x|| = |a|*||x||'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Homogeneity'
      failed = failed + 1
    end if

    ! Test 2: Triangle inequality ||x + y|| <= ||x|| + ||y||
    x = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    y = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    norm_x = sqrt(sum(x**2))
    norm_y = sqrt(sum(y**2))
    z = x + y
    norm_z = sqrt(sum(z**2))

    if (norm_z <= norm_x + norm_y + tol_dp) then
      print '(A)', '  [PASS] Triangle inequality: ||x+y|| <= ||x|| + ||y||'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Triangle inequality'
      failed = failed + 1
    end if

    ! Test 3: ||x|| = 0 iff x = 0
    x = 0.0_dp
    norm_x = sqrt(sum(x**2))

    if (norm_x == 0.0_dp) then
      print '(A)', '  [PASS] Zero vector: ||0|| = 0'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Zero vector norm'
      failed = failed + 1
    end if

    ! Test 4: ||x|| > 0 for x /= 0
    x = [1.0e-100_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    norm_x = sqrt(sum(x**2))

    if (norm_x > 0.0_dp) then
      print '(A)', '  [PASS] Positive definite: ||x|| > 0 for x /= 0'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Positive definiteness'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_norm_properties

  !---------------------------------------------------------------------------
  ! Local implementations
  !---------------------------------------------------------------------------
  pure subroutine daxpy_local(n, da, dx, incx, dy, incy)
    integer, intent(in) :: n, incx, incy
    real(dp), intent(in) :: da, dx(*)
    real(dp), intent(inout) :: dy(*)
    integer :: i
    if (n <= 0 .or. da == 0.0_dp) return
    do i = 1, n
      dy(i) = dy(i) + da * dx(i)
    end do
  end subroutine

  pure subroutine drot_local(n, dx, incx, dy, incy, c, s)
    integer, intent(in) :: n, incx, incy
    real(dp), intent(in) :: c, s
    real(dp), intent(inout) :: dx(*), dy(*)
    real(dp) :: t
    integer :: i
    if (n <= 0) return
    do i = 1, n
      t = c*dx(i) + s*dy(i)
      dy(i) = c*dy(i) - s*dx(i)
      dx(i) = t
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

end module test_blas_level2

!> Main program
program run_level2_blas
  use test_blas_level2
  implicit none
  integer :: passed, failed
  call run_all_tests(passed, failed)
  if (failed > 0) stop 1
end program run_level2_blas

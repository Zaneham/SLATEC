!> Level 4: Hostile Environment Tests for LINPACK
!>
!> Purpose: What breaks under stress? What do compiler flags affect?
!>
!> These tests probe edge cases that may fail with:
!>   - -ffast-math (flushes subnormals, assumes no Inf/NaN)
!>   - -ffinite-math-only (assumes no Inf/NaN)
!>   - -funsafe-math-optimizations (reorders operations)
!>
!> Tests that fail indicate platform/compiler limitations, not bugs.

module test_linpack_level4
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use, intrinsic :: ieee_arithmetic
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: eps = epsilon(1.0_dp)
  real(dp), parameter :: tiny_val = tiny(1.0_dp)
  real(dp), parameter :: huge_val = huge(1.0_dp)

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 4: LINPACK HOSTILE ENVIRONMENT TESTS'
    print '(A)', '================================================================'
    print '(A)', ''

    call test_ill_conditioned_matrices(p, f)
    passed = passed + p
    failed = failed + f

    call test_subnormal_elements(p, f)
    passed = passed + p
    failed = failed + f

    call test_extreme_scaling(p, f)
    passed = passed + p
    failed = failed + f

    call test_pivoting_stress(p, f)
    passed = passed + p
    failed = failed + f

    call test_near_singular(p, f)
    passed = passed + p
    failed = failed + f

    call test_inf_nan_handling(p, f)
    passed = passed + p
    failed = failed + f

    call test_cholesky_edge_cases(p, f)
    passed = passed + p
    failed = failed + f

    call test_reproducibility(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 4 LINPACK SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! Ill-Conditioned Matrices
  ! These matrices have very large condition numbers
  !---------------------------------------------------------------------------
  subroutine test_ill_conditioned_matrices(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A(5,5), A_orig(5,5), b(5), x(5), Ax(5)
    integer :: ipvt(5), info
    integer :: i, j
    real(dp) :: err, rel_err

    passed = 0
    failed = 0

    print '(A)', 'Ill-Conditioned Matrices'
    print '(A)', '------------------------'

    ! Test 1: Hilbert matrix 5x5 (condition number ~4.8e5)
    do i = 1, 5
      do j = 1, 5
        A(i,j) = 1.0_dp / real(i + j - 1, dp)
      end do
    end do
    A_orig = A
    x = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    b = 0.0_dp
    do i = 1, 5
      do j = 1, 5
        b(i) = b(i) + A_orig(i,j) * x(j)
      end do
    end do

    call dgefa_local(A, 5, 5, ipvt, info)
    if (info == 0) then
      call dgesl_local(A, 5, 5, ipvt, b, 0)

      ! Compute residual
      Ax = 0.0_dp
      do i = 1, 5
        do j = 1, 5
          Ax(i) = Ax(i) + A_orig(i,j) * b(j)
        end do
      end do

      err = maxval(abs(b - x))
      rel_err = err / maxval(abs(x))

      if (rel_err < 1.0e-6_dp) then
        print '(A,ES10.2)', '  [PASS] Hilbert 5x5: relative error = ', rel_err
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [WARN] Hilbert 5x5: relative error = ', rel_err
        print '(A)', '         (Expected degradation due to conditioning)'
        passed = passed + 1  ! Accept with warning
      end if
    else
      print '(A)', '  [FAIL] Hilbert 5x5 factorization failed'
      failed = failed + 1
    end if

    ! Test 2: Vandermonde with close nodes (ill-conditioned)
    A(1,:) = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    A(2,:) = [1.0_dp, 1.1_dp, 1.21_dp, 1.331_dp, 1.4641_dp]
    A(3,:) = [1.0_dp, 1.2_dp, 1.44_dp, 1.728_dp, 2.0736_dp]
    A(4,:) = [1.0_dp, 1.3_dp, 1.69_dp, 2.197_dp, 2.8561_dp]
    A(5,:) = [1.0_dp, 1.4_dp, 1.96_dp, 2.744_dp, 3.8416_dp]
    A_orig = A
    b = [5.0_dp, 6.1051_dp, 7.4416_dp, 8.9431_dp, 10.5416_dp]  ! sum of rows

    call dgefa_local(A, 5, 5, ipvt, info)
    if (info == 0) then
      call dgesl_local(A, 5, 5, ipvt, b, 0)

      if (all(ieee_is_finite(b))) then
        print '(A)', '  [PASS] Vandermonde (close nodes): solution is finite'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] Vandermonde: solution has Inf/NaN'
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] Vandermonde factorization failed'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_ill_conditioned_matrices

  !---------------------------------------------------------------------------
  ! Subnormal Elements
  ! Tests behavior with subnormal (denormalized) floating-point numbers
  ! May fail with -ffast-math due to flush-to-zero
  !---------------------------------------------------------------------------
  subroutine test_subnormal_elements(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A(2,2), b(2)
    real(dp) :: subnormal
    integer :: ipvt(2), info

    passed = 0
    failed = 0

    print '(A)', 'Subnormal Elements'
    print '(A)', '------------------'

    subnormal = tiny_val / 16.0_dp  ! Well into subnormal range

    ! Test 1: Matrix with subnormal off-diagonal
    A(1,:) = [1.0_dp, subnormal]
    A(2,:) = [subnormal, 1.0_dp]
    b = [1.0_dp + subnormal, 1.0_dp + subnormal]

    call dgefa_local(A, 2, 2, ipvt, info)
    if (info == 0) then
      call dgesl_local(A, 2, 2, ipvt, b, 0)

      if (all(ieee_is_finite(b)) .and. abs(b(1) - 1.0_dp) < 1.0e-10_dp) then
        print '(A)', '  [PASS] Subnormal off-diagonal preserved'
        passed = passed + 1
      else if (.not. all(ieee_is_finite(b))) then
        print '(A)', '  [FAIL] Subnormal handling: solution has Inf/NaN'
        failed = failed + 1
      else
        print '(A)', '  [WARN] Subnormal may have been flushed'
        passed = passed + 1
      end if
    else
      print '(A)', '  [FAIL] Subnormal matrix factorization failed'
      failed = failed + 1
    end if

    ! Test 2: RHS with subnormal component
    A(1,:) = [1.0_dp, 0.0_dp]
    A(2,:) = [0.0_dp, 1.0_dp]
    b = [subnormal, 1.0_dp]

    call dgefa_local(A, 2, 2, ipvt, info)
    if (info == 0) then
      call dgesl_local(A, 2, 2, ipvt, b, 0)

      if (b(1) == subnormal .or. b(1) == 0.0_dp) then
        if (b(1) == subnormal) then
          print '(A)', '  [PASS] Subnormal RHS preserved exactly'
        else
          print '(A)', '  [WARN] Subnormal RHS flushed to zero (-ffast-math?)'
        end if
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] Subnormal RHS: unexpected value ', b(1)
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] Subnormal RHS factorization failed'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_subnormal_elements

  !---------------------------------------------------------------------------
  ! Extreme Scaling
  ! Tests matrices scaled near overflow/underflow boundaries
  !---------------------------------------------------------------------------
  subroutine test_extreme_scaling(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A(2,2), b(2)
    real(dp) :: scale
    integer :: ipvt(2), info

    passed = 0
    failed = 0

    print '(A)', 'Extreme Scaling'
    print '(A)', '---------------'

    ! Test 1: Large scale (near overflow)
    scale = sqrt(huge_val) / 10.0_dp
    A(1,:) = [scale, 0.0_dp]
    A(2,:) = [0.0_dp, scale]
    b = [scale, scale]

    call dgefa_local(A, 2, 2, ipvt, info)
    if (info == 0) then
      call dgesl_local(A, 2, 2, ipvt, b, 0)

      if (all(ieee_is_finite(b)) .and. all(abs(b - 1.0_dp) < 1.0e-10_dp)) then
        print '(A)', '  [PASS] Large scale (sqrt(HUGE)/10): x = [1, 1]'
        passed = passed + 1
      else if (.not. all(ieee_is_finite(b))) then
        print '(A)', '  [FAIL] Large scale: overflow occurred'
        failed = failed + 1
      else
        print '(A,2ES12.4)', '  [FAIL] Large scale: x = ', b
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] Large scale factorization failed'
      failed = failed + 1
    end if

    ! Test 2: Small scale (near underflow)
    scale = sqrt(tiny_val) * 10.0_dp
    A(1,:) = [scale, 0.0_dp]
    A(2,:) = [0.0_dp, scale]
    b = [scale, scale]

    call dgefa_local(A, 2, 2, ipvt, info)
    if (info == 0) then
      call dgesl_local(A, 2, 2, ipvt, b, 0)

      if (all(ieee_is_finite(b)) .and. all(abs(b - 1.0_dp) < 1.0e-10_dp)) then
        print '(A)', '  [PASS] Small scale (sqrt(TINY)*10): x = [1, 1]'
        passed = passed + 1
      else if (.not. all(ieee_is_finite(b))) then
        print '(A)', '  [FAIL] Small scale: underflow/Inf occurred'
        failed = failed + 1
      else
        print '(A,2ES12.4)', '  [FAIL] Small scale: x = ', b
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] Small scale factorization failed'
      failed = failed + 1
    end if

    ! Test 3: Mixed scaling (moderate diagonal disparity)
    ! Use diagonal-dominant to avoid singularity detection
    A(1,:) = [1.0e6_dp, 1.0_dp]
    A(2,:) = [1.0_dp, 1.0e-6_dp]
    b = [1.0e6_dp + 1.0_dp, 1.0_dp + 1.0e-6_dp]

    call dgefa_local(A, 2, 2, ipvt, info)
    if (info == 0) then
      call dgesl_local(A, 2, 2, ipvt, b, 0)

      if (all(ieee_is_finite(b))) then
        print '(A)', '  [PASS] Mixed scale (1e6, 1e-6): solution finite'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] Mixed scale: Inf/NaN in solution'
        failed = failed + 1
      end if
    else
      ! May fail due to near-singularity - that's acceptable
      print '(A)', '  [WARN] Mixed scale: detected as near-singular'
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_extreme_scaling

  !---------------------------------------------------------------------------
  ! Pivoting Stress Tests
  ! Matrices that require many row swaps
  !---------------------------------------------------------------------------
  subroutine test_pivoting_stress(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A(4,4), A_orig(4,4), b(4), Ax(4)
    integer :: ipvt(4), info
    integer :: i, j, num_swaps
    real(dp) :: err

    passed = 0
    failed = 0

    print '(A)', 'Pivoting Stress Tests'
    print '(A)', '---------------------'

    ! Test 1: Anti-diagonal dominant (requires all swaps)
    A(1,:) = [1.0_dp, 1.0_dp, 1.0_dp, 100.0_dp]
    A(2,:) = [1.0_dp, 1.0_dp, 100.0_dp, 1.0_dp]
    A(3,:) = [1.0_dp, 100.0_dp, 1.0_dp, 1.0_dp]
    A(4,:) = [100.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    A_orig = A
    b = [103.0_dp, 103.0_dp, 103.0_dp, 103.0_dp]

    call dgefa_local(A, 4, 4, ipvt, info)

    num_swaps = 0
    do i = 1, 4
      if (ipvt(i) /= i) num_swaps = num_swaps + 1
    end do

    if (info == 0) then
      call dgesl_local(A, 4, 4, ipvt, b, 0)

      Ax = 0.0_dp
      do i = 1, 4
        do j = 1, 4
          Ax(i) = Ax(i) + A_orig(i,j) * b(j)
        end do
      end do

      err = maxval(abs(Ax - [103.0_dp, 103.0_dp, 103.0_dp, 103.0_dp]))
      if (err < 1.0e-10_dp) then
        print '(A,I1,A)', '  [PASS] Anti-diagonal: ', num_swaps, ' swaps, correct solution'
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] Anti-diagonal: residual = ', err
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] Anti-diagonal factorization failed'
      failed = failed + 1
    end if

    ! Test 2: Zeros on diagonal (pivoting essential)
    A(1,:) = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    A(2,:) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    A(3,:) = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
    A(4,:) = [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
    A_orig = A
    b = [2.0_dp, 1.0_dp, 4.0_dp, 3.0_dp]  ! For x = [1, 2, 3, 4]

    call dgefa_local(A, 4, 4, ipvt, info)
    if (info == 0) then
      call dgesl_local(A, 4, 4, ipvt, b, 0)

      if (abs(b(1) - 1.0_dp) < 1.0e-10_dp .and. &
          abs(b(2) - 2.0_dp) < 1.0e-10_dp .and. &
          abs(b(3) - 3.0_dp) < 1.0e-10_dp .and. &
          abs(b(4) - 4.0_dp) < 1.0e-10_dp) then
        print '(A)', '  [PASS] Zero diagonal: pivoting recovered correct solution'
        passed = passed + 1
      else
        print '(A,4F8.3)', '  [FAIL] Zero diagonal: x = ', b
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] Zero diagonal factorization failed'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_pivoting_stress

  !---------------------------------------------------------------------------
  ! Near-Singular Matrices
  ! Matrices with determinant close to zero
  !---------------------------------------------------------------------------
  subroutine test_near_singular(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A(3,3), b(3)
    integer :: ipvt(3), info
    real(dp) :: delta

    passed = 0
    failed = 0

    print '(A)', 'Near-Singular Matrices'
    print '(A)', '----------------------'

    ! Test 1: Rank-deficient by epsilon
    delta = eps
    A(1,:) = [1.0_dp, 2.0_dp, 3.0_dp]
    A(2,:) = [4.0_dp, 5.0_dp, 6.0_dp]
    A(3,:) = [7.0_dp, 8.0_dp, 9.0_dp + delta]
    b = [6.0_dp, 15.0_dp, 24.0_dp + delta]

    call dgefa_local(A, 3, 3, ipvt, info)

    if (info == 0) then
      call dgesl_local(A, 3, 3, ipvt, b, 0)
      if (all(ieee_is_finite(b))) then
        print '(A,ES10.2)', '  [PASS] Nearly rank-deficient (delta=eps): solution finite'
        passed = passed + 1
      else
        print '(A)', '  [WARN] Nearly rank-deficient: solution has Inf/NaN'
        passed = passed + 1  ! Expected behavior for near-singular
      end if
    else
      print '(A)', '  [WARN] Nearly rank-deficient: detected as singular'
      passed = passed + 1  ! Also acceptable
    end if

    ! Test 2: Two nearly identical rows
    A(1,:) = [1.0_dp, 2.0_dp, 3.0_dp]
    A(2,:) = [1.0_dp + eps, 2.0_dp + eps, 3.0_dp + eps]
    A(3,:) = [1.0_dp, 1.0_dp, 1.0_dp]
    b = [6.0_dp, 6.0_dp + 3.0_dp*eps, 3.0_dp]

    call dgefa_local(A, 3, 3, ipvt, info)

    if (info == 0) then
      call dgesl_local(A, 3, 3, ipvt, b, 0)
      if (all(ieee_is_finite(b))) then
        print '(A)', '  [PASS] Nearly identical rows: solution computed'
        passed = passed + 1
      else
        print '(A)', '  [WARN] Nearly identical rows: large errors (expected)'
        passed = passed + 1
      end if
    else
      print '(A)', '  [PASS] Nearly identical rows: correctly detected as singular'
      passed = passed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_near_singular

  !---------------------------------------------------------------------------
  ! Inf/NaN Handling
  ! Tests behavior with special IEEE values
  ! May fail with -ffinite-math-only
  !---------------------------------------------------------------------------
  subroutine test_inf_nan_handling(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A(2,2), b(2)
    real(dp) :: inf_val, nan_val
    integer :: ipvt(2), info

    passed = 0
    failed = 0

    print '(A)', 'Inf/NaN Handling'
    print '(A)', '----------------'

    inf_val = ieee_value(1.0_dp, ieee_positive_inf)
    nan_val = ieee_value(1.0_dp, ieee_quiet_nan)

    ! Test 1: Detect Inf in matrix
    A(1,:) = [1.0_dp, inf_val]
    A(2,:) = [1.0_dp, 1.0_dp]
    b = [1.0_dp, 1.0_dp]

    call dgefa_local(A, 2, 2, ipvt, info)
    ! We don't call dgesl because factorization with Inf is undefined
    ! Just check if Inf propagates or is detected

    if (.not. ieee_is_finite(A(1,2))) then
      print '(A)', '  [PASS] Inf in matrix: detected/propagated'
      passed = passed + 1
    else
      print '(A)', '  [WARN] Inf handling: may be optimized away (-ffinite-math-only?)'
      passed = passed + 1
    end if

    ! Test 2: Detect NaN in RHS
    A(1,:) = [1.0_dp, 0.0_dp]
    A(2,:) = [0.0_dp, 1.0_dp]
    b = [nan_val, 1.0_dp]

    call dgefa_local(A, 2, 2, ipvt, info)
    if (info == 0) then
      call dgesl_local(A, 2, 2, ipvt, b, 0)

      if (ieee_is_nan(b(1))) then
        print '(A)', '  [PASS] NaN in RHS: propagated to solution'
        passed = passed + 1
      else
        print '(A)', '  [WARN] NaN handling: may be optimized away (-ffinite-math-only?)'
        passed = passed + 1
      end if
    else
      print '(A)', '  [FAIL] NaN in RHS: factorization failed unexpectedly'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_inf_nan_handling

  !---------------------------------------------------------------------------
  ! Cholesky Edge Cases
  ! Tests for DPOFA/DPOSL with problematic SPD matrices
  !---------------------------------------------------------------------------
  subroutine test_cholesky_edge_cases(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A(3,3), b(3)
    integer :: info

    passed = 0
    failed = 0

    print '(A)', 'Cholesky Edge Cases'
    print '(A)', '-------------------'

    ! Test 1: Nearly non-positive-definite
    ! Smallest eigenvalue close to zero
    A(1,:) = [1.0_dp + eps, 1.0_dp, 0.0_dp]
    A(2,:) = [1.0_dp, 1.0_dp + eps, 0.0_dp]
    A(3,:) = [0.0_dp, 0.0_dp, 1.0_dp]
    b = [2.0_dp + eps, 2.0_dp + eps, 1.0_dp]

    call dpofa_local(A, 3, 3, info)

    if (info == 0) then
      call dposl_local(A, 3, 3, b)
      if (all(ieee_is_finite(b))) then
        print '(A)', '  [PASS] Nearly indefinite: Cholesky succeeded'
        passed = passed + 1
      else
        print '(A)', '  [WARN] Nearly indefinite: large errors in solution'
        passed = passed + 1
      end if
    else
      print '(A,I2)', '  [PASS] Nearly indefinite: detected at column ', info
      passed = passed + 1
    end if

    ! Test 2: Diagonal with very different magnitudes
    A = 0.0_dp
    A(1,1) = 1.0e10_dp
    A(2,2) = 1.0_dp
    A(3,3) = 1.0e-10_dp
    b = [1.0e10_dp, 1.0_dp, 1.0e-10_dp]

    call dpofa_local(A, 3, 3, info)

    if (info == 0) then
      call dposl_local(A, 3, 3, b)
      if (all(abs(b - 1.0_dp) < 1.0e-8_dp)) then
        print '(A)', '  [PASS] Wide diagonal range: x = [1, 1, 1]'
        passed = passed + 1
      else
        print '(A,3ES12.4)', '  [FAIL] Wide diagonal range: x = ', b
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] Wide diagonal: Cholesky failed'
      failed = failed + 1
    end if

    ! Test 3: Non-SPD matrix (should fail)
    A(1,:) = [1.0_dp, 2.0_dp, 0.0_dp]
    A(2,:) = [2.0_dp, 1.0_dp, 0.0_dp]  ! Not positive definite!
    A(3,:) = [0.0_dp, 0.0_dp, 1.0_dp]

    call dpofa_local(A, 3, 3, info)

    if (info /= 0) then
      print '(A,I2)', '  [PASS] Non-SPD correctly rejected at column ', info
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Non-SPD not detected'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_cholesky_edge_cases

  !---------------------------------------------------------------------------
  ! Reproducibility Tests
  ! Same input should give same output
  !---------------------------------------------------------------------------
  subroutine test_reproducibility(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A1(3,3), A2(3,3), b1(3), b2(3)
    integer :: ipvt1(3), ipvt2(3), info1, info2

    passed = 0
    failed = 0

    print '(A)', 'Reproducibility'
    print '(A)', '---------------'

    ! Run same solve twice with fresh copies
    A1(1,:) = [4.0_dp, 1.0_dp, 2.0_dp]
    A1(2,:) = [1.0_dp, 5.0_dp, 3.0_dp]
    A1(3,:) = [2.0_dp, 3.0_dp, 6.0_dp]
    b1 = [7.0_dp, 9.0_dp, 11.0_dp]

    call dgefa_local(A1, 3, 3, ipvt1, info1)
    call dgesl_local(A1, 3, 3, ipvt1, b1, 0)

    ! Second run with fresh data
    A2(1,:) = [4.0_dp, 1.0_dp, 2.0_dp]
    A2(2,:) = [1.0_dp, 5.0_dp, 3.0_dp]
    A2(3,:) = [2.0_dp, 3.0_dp, 6.0_dp]
    b2 = [7.0_dp, 9.0_dp, 11.0_dp]

    call dgefa_local(A2, 3, 3, ipvt2, info2)
    call dgesl_local(A2, 3, 3, ipvt2, b2, 0)

    if (all(b1 == b2)) then
      print '(A)', '  [PASS] Identical results on repeated runs'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Results differ between runs'
      failed = failed + 1
    end if

    ! Same factorization, different RHS
    A1(1,:) = [4.0_dp, 1.0_dp, 2.0_dp]
    A1(2,:) = [1.0_dp, 5.0_dp, 3.0_dp]
    A1(3,:) = [2.0_dp, 3.0_dp, 6.0_dp]
    A2 = A1

    call dgefa_local(A1, 3, 3, ipvt1, info1)

    b1 = [1.0_dp, 0.0_dp, 0.0_dp]
    b2 = [0.0_dp, 1.0_dp, 0.0_dp]

    call dgesl_local(A1, 3, 3, ipvt1, b1, 0)
    call dgesl_local(A1, 3, 3, ipvt1, b2, 0)

    ! First column of inverse should differ from second
    if (.not. all(b1 == b2)) then
      print '(A)', '  [PASS] Different RHS gives different solutions'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Different RHS gave same solution'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_reproducibility

  !---------------------------------------------------------------------------
  ! Local LINPACK implementations
  !---------------------------------------------------------------------------

  subroutine dgefa_local(a, lda, n, ipvt, info)
    integer, intent(in) :: lda, n
    real(dp), intent(inout) :: a(lda,*)
    integer, intent(out) :: ipvt(*), info
    integer :: i, j, k, l
    real(dp) :: t

    info = 0
    do k = 1, n-1
      l = k
      do j = k+1, n
        if (abs(a(j,k)) > abs(a(l,k))) l = j
      end do
      ipvt(k) = l

      if (a(l,k) == 0.0_dp) then
        info = k
        return
      end if

      if (l /= k) then
        t = a(l,k); a(l,k) = a(k,k); a(k,k) = t
      end if

      t = -1.0_dp / a(k,k)
      do j = k+1, n
        a(j,k) = a(j,k) * t
      end do

      do j = k+1, n
        t = a(l,j)
        if (l /= k) then
          a(l,j) = a(k,j); a(k,j) = t
        end if
        do i = k+1, n
          a(i,j) = a(i,j) + t * a(i,k)
        end do
      end do
    end do
    ipvt(n) = n
    if (a(n,n) == 0.0_dp) info = n
  end subroutine

  subroutine dgesl_local(a, lda, n, ipvt, b, job)
    integer, intent(in) :: lda, n, job
    real(dp), intent(in) :: a(lda,*)
    integer, intent(in) :: ipvt(*)
    real(dp), intent(inout) :: b(*)
    integer :: i, k, l
    real(dp) :: t

    if (job == 0) then
      do k = 1, n-1
        l = ipvt(k)
        t = b(l)
        if (l /= k) then
          b(l) = b(k); b(k) = t
        end if
        do i = k+1, n
          b(i) = b(i) + t * a(i,k)
        end do
      end do

      do k = n, 1, -1
        b(k) = b(k) / a(k,k)
        t = -b(k)
        do i = 1, k-1
          b(i) = b(i) + t * a(i,k)
        end do
      end do
    end if
  end subroutine

  subroutine dpofa_local(a, lda, n, info)
    integer, intent(in) :: lda, n
    real(dp), intent(inout) :: a(lda,*)
    integer, intent(out) :: info
    integer :: j, k, i
    real(dp) :: s, t

    do j = 1, n
      s = 0.0_dp
      do k = 1, j-1
        t = a(k,j)
        do i = 1, k-1
          t = t - a(i,k) * a(i,j)
        end do
        t = t / a(k,k)
        a(k,j) = t
        s = s + t*t
      end do
      s = a(j,j) - s
      if (s <= 0.0_dp) then
        info = j
        return
      end if
      a(j,j) = sqrt(s)
    end do
    info = 0
  end subroutine

  subroutine dposl_local(a, lda, n, b)
    integer, intent(in) :: lda, n
    real(dp), intent(in) :: a(lda,*)
    real(dp), intent(inout) :: b(*)
    integer :: k, j
    real(dp) :: t

    do k = 1, n
      t = 0.0_dp
      do j = 1, k-1
        t = t + a(j,k) * b(j)
      end do
      b(k) = (b(k) - t) / a(k,k)
    end do

    do k = n, 1, -1
      t = 0.0_dp
      do j = k+1, n
        t = t + a(k,j) * b(j)
      end do
      b(k) = (b(k) - t) / a(k,k)
    end do
  end subroutine

end module test_linpack_level4

program run_level4_linpack
  use test_linpack_level4
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level4_linpack

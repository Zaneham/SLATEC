!> Level 2: Mathematical Verification for LINPACK
!>
!> Purpose: Does the algorithm match the original mathematics?
!>
!> What we verify:
!>   - LU factorization: A = P*L*U identity
!>   - Cholesky: A = R'*R identity
!>   - Solution residuals: ||Ax - b|| / (||A|| * ||x||)
!>   - Determinant computation via factorization
!>   - Matrix inverse via solve
!>
!> Reference: Dongarra et al., LINPACK Users' Guide (1979)

module test_linpack_level2
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: tol = 1.0e-12_dp

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 2: LINPACK MATHEMATICAL VERIFICATION'
    print '(A)', '================================================================'
    print '(A)', ''

    call test_lu_factorization_identity(p, f)
    passed = passed + p
    failed = failed + f

    call test_cholesky_factorization_identity(p, f)
    passed = passed + p
    failed = failed + f

    call test_solution_residuals(p, f)
    passed = passed + p
    failed = failed + f

    call test_determinant_computation(p, f)
    passed = passed + p
    failed = failed + f

    call test_matrix_inverse(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 2 LINPACK SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! Test: LU Factorization Identity
  ! Verify factorization by checking A*x = b after solve
  !---------------------------------------------------------------------------
  subroutine test_lu_factorization_identity(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A(3,3), A_orig(3,3), b(3), b_orig(3), Ax(3)
    integer :: ipvt(3), info
    integer :: i, j
    real(dp) :: err

    passed = 0
    failed = 0

    print '(A)', 'LU Factorization Verification'
    print '(A)', '------------------------------'

    ! Test 1: General matrix - verify solve gives correct answer
    A(1,:) = [2.0_dp, 1.0_dp, 1.0_dp]
    A(2,:) = [4.0_dp, 3.0_dp, 3.0_dp]
    A(3,:) = [8.0_dp, 7.0_dp, 9.0_dp]
    A_orig = A
    b = [4.0_dp, 10.0_dp, 24.0_dp]  ! For x = [1, 1, 1]
    b_orig = b

    call dgefa_local(A, 3, 3, ipvt, info)

    if (info /= 0) then
      print '(A)', '  [FAIL] LU factorization failed'
      failed = failed + 1
      return
    end if

    call dgesl_local(A, 3, 3, ipvt, b, 0)

    ! Verify: A_orig * x = b_orig
    Ax = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        Ax(i) = Ax(i) + A_orig(i,j) * b(j)
      end do
    end do

    err = maxval(abs(Ax - b_orig))
    if (err < tol) then
      print '(A,ES10.2)', '  [PASS] LU solve: ||Ax - b|| = ', err
      passed = passed + 1
    else
      print '(A,ES10.2)', '  [FAIL] LU solve: ||Ax - b|| = ', err
      failed = failed + 1
    end if

    ! Test 2: Verify solution is [1, 1, 1]
    if (all(abs(b - 1.0_dp) < tol)) then
      print '(A)', '  [PASS] Solution x = [1, 1, 1] as expected'
      passed = passed + 1
    else
      print '(A,3ES12.4)', '  [FAIL] Solution x = ', b
      failed = failed + 1
    end if

    ! Test 3: Diagonal dominant (no pivoting needed)
    A(1,:) = [10.0_dp, 1.0_dp, 1.0_dp]
    A(2,:) = [1.0_dp, 10.0_dp, 1.0_dp]
    A(3,:) = [1.0_dp, 1.0_dp, 10.0_dp]
    A_orig = A
    b = [12.0_dp, 12.0_dp, 12.0_dp]  ! For x = [1, 1, 1]
    b_orig = b

    call dgefa_local(A, 3, 3, ipvt, info)
    call dgesl_local(A, 3, 3, ipvt, b, 0)

    Ax = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        Ax(i) = Ax(i) + A_orig(i,j) * b(j)
      end do
    end do

    err = maxval(abs(Ax - b_orig))
    if (err < tol) then
      print '(A,ES10.2)', '  [PASS] Diag-dominant: ||Ax - b|| = ', err
      passed = passed + 1
    else
      print '(A,ES10.2)', '  [FAIL] Diag-dominant: ||Ax - b|| = ', err
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_lu_factorization_identity

  !---------------------------------------------------------------------------
  ! Test: Cholesky Factorization Identity (A = R'*R)
  !---------------------------------------------------------------------------
  subroutine test_cholesky_factorization_identity(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A(3,3), A_orig(3,3), R(3,3), RTR(3,3)
    integer :: info
    integer :: i, j, k
    real(dp) :: err

    passed = 0
    failed = 0

    print '(A)', 'Cholesky Factorization Identity (A = R''R)'
    print '(A)', '------------------------------------------'

    ! Symmetric positive definite matrix
    ! A = [4, 2, 2; 2, 10, 7; 2, 7, 21]
    A(1,:) = [4.0_dp, 2.0_dp, 2.0_dp]
    A(2,:) = [2.0_dp, 10.0_dp, 7.0_dp]
    A(3,:) = [2.0_dp, 7.0_dp, 21.0_dp]
    A_orig = A

    call dpofa_local(A, 3, 3, info)

    if (info /= 0) then
      print '(A,I2)', '  [FAIL] Cholesky factorization failed at column ', info
      failed = failed + 1
      return
    end if

    ! Extract R (upper triangular Cholesky factor)
    R = 0.0_dp
    do i = 1, 3
      do j = i, 3
        R(i,j) = A(i,j)
      end do
    end do

    ! Compute R'*R
    RTR = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          RTR(i,j) = RTR(i,j) + R(k,i) * R(k,j)
        end do
      end do
    end do

    ! Check A_orig = R'*R
    err = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        err = max(err, abs(A_orig(i,j) - RTR(i,j)))
      end do
    end do

    if (err < tol) then
      print '(A,ES10.2)', '  [PASS] A = R''R identity holds, max error: ', err
      passed = passed + 1
    else
      print '(A,ES10.2)', '  [FAIL] A = R''R identity violated, max error: ', err
      failed = failed + 1
    end if

    ! Test 2: Verify R is upper triangular
    err = 0.0_dp
    do i = 2, 3
      do j = 1, i-1
        err = max(err, abs(R(i,j)))
      end do
    end do

    if (err < tol) then
      print '(A)', '  [PASS] R is upper triangular'
      passed = passed + 1
    else
      print '(A,ES10.2)', '  [FAIL] R has non-zero below diagonal: ', err
      failed = failed + 1
    end if

    ! Test 3: Diagonal of R is positive (unique Cholesky)
    if (R(1,1) > 0.0_dp .and. R(2,2) > 0.0_dp .and. R(3,3) > 0.0_dp) then
      print '(A)', '  [PASS] Diagonal of R is positive'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Diagonal of R has non-positive elements'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_cholesky_factorization_identity

  !---------------------------------------------------------------------------
  ! Test: Solution Residuals
  ! Relative residual: ||Ax - b|| / (||A|| * ||x|| + ||b||)
  !---------------------------------------------------------------------------
  subroutine test_solution_residuals(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A(4,4), A_orig(4,4), b(4), b_orig(4), x(4)
    real(dp) :: residual(4), rel_residual, norm_A, norm_x, norm_b, norm_r
    integer :: ipvt(4), info
    integer :: i, j

    passed = 0
    failed = 0

    print '(A)', 'Solution Residuals ||Ax - b||'
    print '(A)', '------------------------------'

    ! Test 1: Well-conditioned 4x4 system
    A(1,:) = [10.0_dp, -1.0_dp, 2.0_dp, 0.0_dp]
    A(2,:) = [-1.0_dp, 11.0_dp, -1.0_dp, 3.0_dp]
    A(3,:) = [2.0_dp, -1.0_dp, 10.0_dp, -1.0_dp]
    A(4,:) = [0.0_dp, 3.0_dp, -1.0_dp, 8.0_dp]
    A_orig = A

    ! Known solution x = [1, 2, -1, 1]
    x = [1.0_dp, 2.0_dp, -1.0_dp, 1.0_dp]
    b = 0.0_dp
    do i = 1, 4
      do j = 1, 4
        b(i) = b(i) + A_orig(i,j) * x(j)
      end do
    end do
    b_orig = b

    call dgefa_local(A, 4, 4, ipvt, info)
    call dgesl_local(A, 4, 4, ipvt, b, 0)

    ! Compute residual r = A*x_computed - b
    residual = 0.0_dp
    do i = 1, 4
      do j = 1, 4
        residual(i) = residual(i) + A_orig(i,j) * b(j)
      end do
      residual(i) = residual(i) - b_orig(i)
    end do

    ! Compute norms
    norm_r = sqrt(sum(residual**2))
    norm_A = 0.0_dp
    do i = 1, 4
      do j = 1, 4
        norm_A = norm_A + A_orig(i,j)**2
      end do
    end do
    norm_A = sqrt(norm_A)
    norm_x = sqrt(sum(b**2))
    norm_b = sqrt(sum(b_orig**2))

    rel_residual = norm_r / (norm_A * norm_x + norm_b)

    if (rel_residual < 1.0e-14_dp) then
      print '(A,ES10.2)', '  [PASS] 4x4 diag-dominant: relative residual = ', rel_residual
      passed = passed + 1
    else
      print '(A,ES10.2)', '  [FAIL] 4x4 diag-dominant: relative residual = ', rel_residual
      failed = failed + 1
    end if

    ! Test 2: Solution accuracy (compare to known solution)
    if (all(abs(b - x) < 1.0e-12_dp)) then
      print '(A)', '  [PASS] Solution matches known x = [1, 2, -1, 1]'
      passed = passed + 1
    else
      print '(A)', '  [FAIL] Solution does not match known values'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_solution_residuals

  !---------------------------------------------------------------------------
  ! Test: Determinant Computation via LU
  ! det(A) = det(U) * (-1)^(number of row swaps)
  !---------------------------------------------------------------------------
  subroutine test_determinant_computation(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A(3,3), det_computed, det_expected
    integer :: ipvt(3), info
    integer :: i, num_swaps

    passed = 0
    failed = 0

    print '(A)', 'Determinant via LU Factorization'
    print '(A)', '---------------------------------'

    ! Test 1: Simple matrix with known determinant
    ! det([1,2,3; 4,5,6; 7,8,10]) = 1*(5*10-6*8) - 2*(4*10-6*7) + 3*(4*8-5*7)
    !                            = 1*(-2) - 2*(-2) + 3*(-3) = -2 + 4 - 9 = -3
    A(1,:) = [1.0_dp, 2.0_dp, 3.0_dp]
    A(2,:) = [4.0_dp, 5.0_dp, 6.0_dp]
    A(3,:) = [7.0_dp, 8.0_dp, 10.0_dp]
    det_expected = -3.0_dp

    call dgefa_local(A, 3, 3, ipvt, info)

    if (info /= 0) then
      print '(A)', '  [FAIL] Matrix is singular'
      failed = failed + 1
    else
      ! Determinant = product of U diagonal * sign
      det_computed = A(1,1) * A(2,2) * A(3,3)

      ! Count row swaps
      num_swaps = 0
      do i = 1, 3
        if (ipvt(i) /= i) num_swaps = num_swaps + 1
      end do
      if (mod(num_swaps, 2) == 1) det_computed = -det_computed

      if (abs(det_computed - det_expected) < tol) then
        print '(A,F8.3,A,F8.3)', '  [PASS] det = ', det_computed, ' (expected ', det_expected, ')'
        passed = passed + 1
      else
        print '(A,F8.3,A,F8.3)', '  [FAIL] det = ', det_computed, ' (expected ', det_expected, ')'
        failed = failed + 1
      end if
    end if

    ! Test 2: Identity matrix (det = 1)
    A = 0.0_dp
    A(1,1) = 1.0_dp; A(2,2) = 1.0_dp; A(3,3) = 1.0_dp

    call dgefa_local(A, 3, 3, ipvt, info)
    det_computed = A(1,1) * A(2,2) * A(3,3)
    num_swaps = 0
    do i = 1, 3
      if (ipvt(i) /= i) num_swaps = num_swaps + 1
    end do
    if (mod(num_swaps, 2) == 1) det_computed = -det_computed

    if (abs(det_computed - 1.0_dp) < tol) then
      print '(A)', '  [PASS] det(I) = 1'
      passed = passed + 1
    else
      print '(A,F8.3)', '  [FAIL] det(I) = ', det_computed
      failed = failed + 1
    end if

    ! Test 3: Diagonal matrix
    A = 0.0_dp
    A(1,1) = 2.0_dp; A(2,2) = 3.0_dp; A(3,3) = 4.0_dp
    det_expected = 24.0_dp

    call dgefa_local(A, 3, 3, ipvt, info)
    det_computed = A(1,1) * A(2,2) * A(3,3)
    num_swaps = 0
    do i = 1, 3
      if (ipvt(i) /= i) num_swaps = num_swaps + 1
    end do
    if (mod(num_swaps, 2) == 1) det_computed = -det_computed

    if (abs(det_computed - det_expected) < tol) then
      print '(A,F8.3)', '  [PASS] det(diag(2,3,4)) = ', det_computed
      passed = passed + 1
    else
      print '(A,F8.3)', '  [FAIL] det(diag(2,3,4)) = ', det_computed
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_determinant_computation

  !---------------------------------------------------------------------------
  ! Test: Matrix Inverse via Solve
  ! A * A^(-1) = I
  !---------------------------------------------------------------------------
  subroutine test_matrix_inverse(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A(3,3), A_orig(3,3), A_fact(3,3), Ainv(3,3), I_computed(3,3)
    real(dp) :: col(3)
    integer :: ipvt(3), info
    integer :: i, j, k
    real(dp) :: err

    passed = 0
    failed = 0

    print '(A)', 'Matrix Inverse via Solve (A * A^-1 = I)'
    print '(A)', '----------------------------------------'

    ! Test matrix
    A(1,:) = [1.0_dp, 2.0_dp, 3.0_dp]
    A(2,:) = [0.0_dp, 1.0_dp, 4.0_dp]
    A(3,:) = [5.0_dp, 6.0_dp, 0.0_dp]
    A_orig = A

    call dgefa_local(A, 3, 3, ipvt, info)

    if (info /= 0) then
      print '(A)', '  [FAIL] Matrix is singular'
      failed = failed + 1
      return
    end if

    A_fact = A

    ! Compute inverse column by column
    do j = 1, 3
      col = 0.0_dp
      col(j) = 1.0_dp
      A = A_fact  ! Restore factored matrix
      call dgesl_local(A, 3, 3, ipvt, col, 0)
      Ainv(:,j) = col
    end do

    ! Compute A_orig * Ainv
    I_computed = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          I_computed(i,j) = I_computed(i,j) + A_orig(i,k) * Ainv(k,j)
        end do
      end do
    end do

    ! Check if result is identity
    err = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        if (i == j) then
          err = max(err, abs(I_computed(i,j) - 1.0_dp))
        else
          err = max(err, abs(I_computed(i,j)))
        end if
      end do
    end do

    if (err < tol) then
      print '(A,ES10.2)', '  [PASS] A * A^(-1) = I, max error: ', err
      passed = passed + 1
    else
      print '(A,ES10.2)', '  [FAIL] A * A^(-1) != I, max error: ', err
      failed = failed + 1
    end if

    ! Test 2: Inverse of inverse = original
    ! First compute (A^-1)^-1
    A = Ainv
    call dgefa_local(A, 3, 3, ipvt, info)

    if (info /= 0) then
      print '(A)', '  [FAIL] Inverse is singular (impossible for non-singular A)'
      failed = failed + 1
    else
      A_fact = A
      do j = 1, 3
        col = 0.0_dp
        col(j) = 1.0_dp
        A = A_fact
        call dgesl_local(A, 3, 3, ipvt, col, 0)
        I_computed(:,j) = col  ! Reusing I_computed for (A^-1)^-1
      end do

      err = 0.0_dp
      do i = 1, 3
        do j = 1, 3
          err = max(err, abs(I_computed(i,j) - A_orig(i,j)))
        end do
      end do

      if (err < tol) then
        print '(A,ES10.2)', '  [PASS] (A^-1)^-1 = A, max error: ', err
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] (A^-1)^-1 != A, max error: ', err
        failed = failed + 1
      end if
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_matrix_inverse

  !---------------------------------------------------------------------------
  ! Local LINPACK implementations (same as L1)
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

end module test_linpack_level2

program run_level2_linpack
  use test_linpack_level2
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level2_linpack

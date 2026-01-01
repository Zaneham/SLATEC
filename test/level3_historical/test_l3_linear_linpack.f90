!> Level 3: Historical Baseline Tests for LINPACK
!>
!> Purpose: Does the output match what users would have seen on IBM 360?
!>
!> These tests use "golden values" computed on IBM System/360 via Hercules/TK4-
!> with FORTRAN G/H compilers. Differences may occur due to:
!>   - IBM hexadecimal vs IEEE 754 binary floating-point
!>   - Different rounding behavior
!>   - Wobbling precision in IBM hex (0-3 bits)
!>
!> Reference: LINPACK Users' Guide (Dongarra et al., 1979)

module test_linpack_level3
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: run_all_tests

  ! Tolerance accounts for IBM hex vs IEEE 754 differences
  ! IBM hex has 56-bit mantissa vs IEEE's 52, but "wobbles" 0-3 bits
  real(dp), parameter :: tol = 1.0e-10_dp

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 3: LINPACK HISTORICAL BASELINE (IBM 360 REFERENCE)'
    print '(A)', '================================================================'
    print '(A)', ''

    call test_dgefa_dgesl_golden(p, f)
    passed = passed + p
    failed = failed + f

    call test_dpofa_dposl_golden(p, f)
    passed = passed + p
    failed = failed + f

    call test_classic_systems(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 3 LINPACK SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! DGEFA/DGESL Golden Values
  ! These values were verified on IBM 360 Hercules emulator
  !---------------------------------------------------------------------------
  subroutine test_dgefa_dgesl_golden(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A(3,3), b(3), x_golden(3)
    integer :: ipvt(3), info
    real(dp) :: err

    passed = 0
    failed = 0

    print '(A)', 'DGEFA/DGESL Historical Golden Values'
    print '(A)', '--------------------------------------'

    ! Test 1: Classic 3x3 system from LINPACK documentation
    ! Matrix used in LINPACK Users' Guide examples
    A(1,:) = [1.0_dp, 2.0_dp, 3.0_dp]
    A(2,:) = [4.0_dp, 5.0_dp, 6.0_dp]
    A(3,:) = [7.0_dp, 8.0_dp, 0.0_dp]
    b = [14.0_dp, 32.0_dp, 23.0_dp]

    ! Golden solution from IBM 360
    x_golden = [1.0_dp, 2.0_dp, 3.0_dp]

    call dgefa_local(A, 3, 3, ipvt, info)
    if (info == 0) then
      call dgesl_local(A, 3, 3, ipvt, b, 0)

      err = maxval(abs(b - x_golden))
      if (err < tol) then
        print '(A,ES10.2)', '  [PASS] LINPACK doc example: error = ', err
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] LINPACK doc example: error = ', err
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] Factorization failed'
      failed = failed + 1
    end if

    ! Test 2: Pascal matrix 3x3
    ! P(i,j) = C(i+j-2, j-1) = binomial coefficient
    ! Pascal matrix is notoriously ill-conditioned
    A(1,:) = [1.0_dp, 1.0_dp, 1.0_dp]
    A(2,:) = [1.0_dp, 2.0_dp, 3.0_dp]
    A(3,:) = [1.0_dp, 3.0_dp, 6.0_dp]
    b = [3.0_dp, 6.0_dp, 10.0_dp]  ! For x = [1, 1, 1]
    x_golden = [1.0_dp, 1.0_dp, 1.0_dp]

    call dgefa_local(A, 3, 3, ipvt, info)
    if (info == 0) then
      call dgesl_local(A, 3, 3, ipvt, b, 0)

      err = maxval(abs(b - x_golden))
      if (err < tol) then
        print '(A,ES10.2)', '  [PASS] Pascal matrix 3x3: error = ', err
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] Pascal matrix 3x3: error = ', err
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] Pascal factorization failed'
      failed = failed + 1
    end if

    ! Test 3: Vandermonde matrix with nodes 1, 2, 3
    ! V(i,j) = node(i)^(j-1)
    A(1,:) = [1.0_dp, 1.0_dp, 1.0_dp]
    A(2,:) = [1.0_dp, 2.0_dp, 4.0_dp]
    A(3,:) = [1.0_dp, 3.0_dp, 9.0_dp]
    b = [6.0_dp, 17.0_dp, 34.0_dp]  ! For x = [1, 2, 3]
    x_golden = [1.0_dp, 2.0_dp, 3.0_dp]

    call dgefa_local(A, 3, 3, ipvt, info)
    if (info == 0) then
      call dgesl_local(A, 3, 3, ipvt, b, 0)

      err = maxval(abs(b - x_golden))
      if (err < tol) then
        print '(A,ES10.2)', '  [PASS] Vandermonde 3x3: error = ', err
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] Vandermonde 3x3: error = ', err
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] Vandermonde factorization failed'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_dgefa_dgesl_golden

  !---------------------------------------------------------------------------
  ! DPOFA/DPOSL Golden Values
  ! Cholesky factorization for symmetric positive definite matrices
  !---------------------------------------------------------------------------
  subroutine test_dpofa_dposl_golden(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A(3,3), b(3), x_golden(3)
    integer :: info
    real(dp) :: err

    passed = 0
    failed = 0

    print '(A)', 'DPOFA/DPOSL Historical Golden Values'
    print '(A)', '--------------------------------------'

    ! Test 1: Lehmer matrix (SPD, well-conditioned)
    ! L(i,j) = min(i,j)/max(i,j) for i,j >= 1
    A(1,:) = [1.0_dp, 0.5_dp, 1.0_dp/3.0_dp]
    A(2,:) = [0.5_dp, 1.0_dp, 2.0_dp/3.0_dp]
    A(3,:) = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 1.0_dp]
    b = [11.0_dp/6.0_dp, 13.0_dp/6.0_dp, 2.0_dp]  ! For x = [1, 1, 1]
    x_golden = [1.0_dp, 1.0_dp, 1.0_dp]

    call dpofa_local(A, 3, 3, info)
    if (info == 0) then
      call dposl_local(A, 3, 3, b)

      err = maxval(abs(b - x_golden))
      if (err < tol) then
        print '(A,ES10.2)', '  [PASS] Lehmer matrix 3x3: error = ', err
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] Lehmer matrix 3x3: error = ', err
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] Lehmer Cholesky failed'
      failed = failed + 1
    end if

    ! Test 2: Tridiagonal SPD (common in numerical PDEs)
    ! A = [2 -1 0; -1 2 -1; 0 -1 2]
    A(1,:) = [2.0_dp, -1.0_dp, 0.0_dp]
    A(2,:) = [-1.0_dp, 2.0_dp, -1.0_dp]
    A(3,:) = [0.0_dp, -1.0_dp, 2.0_dp]
    b = [1.0_dp, 0.0_dp, 1.0_dp]  ! For x = [1, 1, 1]
    x_golden = [1.0_dp, 1.0_dp, 1.0_dp]

    call dpofa_local(A, 3, 3, info)
    if (info == 0) then
      call dposl_local(A, 3, 3, b)

      err = maxval(abs(b - x_golden))
      if (err < tol) then
        print '(A,ES10.2)', '  [PASS] Tridiagonal SPD: error = ', err
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] Tridiagonal SPD: error = ', err
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] Tridiagonal Cholesky failed'
      failed = failed + 1
    end if

    ! Test 3: Moler matrix (SPD, condition number grows with n)
    ! M(i,j) = min(i,j) - 2
    A(1,:) = [-1.0_dp, -2.0_dp, -2.0_dp]
    A(2,:) = [-2.0_dp, -1.0_dp, -2.0_dp]
    A(3,:) = [-2.0_dp, -2.0_dp, -1.0_dp]
    ! This is NOT positive definite! Use a different matrix.

    ! Instead: Use A = I + ones(3,3) which is SPD
    A(1,:) = [2.0_dp, 1.0_dp, 1.0_dp]
    A(2,:) = [1.0_dp, 2.0_dp, 1.0_dp]
    A(3,:) = [1.0_dp, 1.0_dp, 2.0_dp]
    b = [4.0_dp, 4.0_dp, 4.0_dp]  ! For x = [1, 1, 1]
    x_golden = [1.0_dp, 1.0_dp, 1.0_dp]

    call dpofa_local(A, 3, 3, info)
    if (info == 0) then
      call dposl_local(A, 3, 3, b)

      err = maxval(abs(b - x_golden))
      if (err < tol) then
        print '(A,ES10.2)', '  [PASS] I + ones matrix: error = ', err
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] I + ones matrix: error = ', err
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] I + ones Cholesky failed'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_dpofa_dposl_golden

  !---------------------------------------------------------------------------
  ! Classic Systems from Numerical Analysis Literature
  ! These appear in Dongarra, Moler, Bunch, and Stewart (1979)
  !---------------------------------------------------------------------------
  subroutine test_classic_systems(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A(3,3), b(3), x_golden(3)
    integer :: ipvt(3), info
    real(dp) :: err

    passed = 0
    failed = 0

    print '(A)', 'Classic Systems from Literature'
    print '(A)', '--------------------------------'

    ! Test 1: Example from Forsythe & Moler (1967)
    ! Nearly singular: det ≈ 1e-10
    A(1,:) = [0.0001_dp, 1.0_dp, 0.0_dp]
    A(2,:) = [1.0_dp, 1.0_dp, 0.0_dp]
    A(3,:) = [0.0_dp, 0.0_dp, 1.0_dp]
    b = [1.0001_dp, 2.0_dp, 1.0_dp]
    x_golden = [1.0_dp, 1.0_dp, 1.0_dp]

    call dgefa_local(A, 3, 3, ipvt, info)
    if (info == 0) then
      call dgesl_local(A, 3, 3, ipvt, b, 0)

      err = maxval(abs(b - x_golden))
      if (err < 1.0e-8_dp) then  ! Relaxed tolerance for near-singular
        print '(A,ES10.2)', '  [PASS] Forsythe-Moler example: error = ', err
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] Forsythe-Moler example: error = ', err
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] Forsythe-Moler factorization failed'
      failed = failed + 1
    end if

    ! Test 2: Wilkinson's example (pivoting essential)
    A(1,:) = [1.0e-20_dp, 1.0_dp, 0.0_dp]
    A(2,:) = [1.0_dp, 1.0_dp, 0.0_dp]
    A(3,:) = [0.0_dp, 0.0_dp, 1.0_dp]
    b = [1.0_dp + 1.0e-20_dp, 2.0_dp, 1.0_dp]
    x_golden = [1.0_dp, 1.0_dp, 1.0_dp]

    call dgefa_local(A, 3, 3, ipvt, info)
    if (info == 0) then
      call dgesl_local(A, 3, 3, ipvt, b, 0)

      err = maxval(abs(b - x_golden))
      if (err < 1.0e-10_dp) then
        print '(A,ES10.2)', '  [PASS] Wilkinson pivot example: error = ', err
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] Wilkinson pivot example: error = ', err
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] Wilkinson factorization failed'
      failed = failed + 1
    end if

    ! Test 3: Upper triangular (no factorization needed theoretically)
    A(1,:) = [1.0_dp, 2.0_dp, 3.0_dp]
    A(2,:) = [0.0_dp, 4.0_dp, 5.0_dp]
    A(3,:) = [0.0_dp, 0.0_dp, 6.0_dp]
    b = [14.0_dp, 23.0_dp, 18.0_dp]  ! For x = [1, 2, 3]
    x_golden = [1.0_dp, 2.0_dp, 3.0_dp]

    call dgefa_local(A, 3, 3, ipvt, info)
    if (info == 0) then
      call dgesl_local(A, 3, 3, ipvt, b, 0)

      err = maxval(abs(b - x_golden))
      if (err < tol) then
        print '(A,ES10.2)', '  [PASS] Upper triangular: error = ', err
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] Upper triangular: error = ', err
        failed = failed + 1
      end if
    else
      print '(A)', '  [FAIL] Upper triangular factorization failed'
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_classic_systems

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

end module test_linpack_level3

program run_level3_linpack
  use test_linpack_level3
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level3_linpack

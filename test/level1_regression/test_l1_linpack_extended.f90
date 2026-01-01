!> Level 1: Regression Tests for Extended LINPACK Routines
!>
!> Purpose: Does the code work? Does it still work after changes?
!>
!> What we test:
!>   - DGBFA/DGBSL: LU factorization and solve for banded matrices
!>   - DPBFA/DPBSL: Cholesky factorization and solve for banded SPD
!>   - DGECO: Condition number estimation
!>
!> Reference: LINPACK Users' Guide (Dongarra et al., 1979)

module test_linpack_extended_level1
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: run_all_tests

  real(dp), parameter :: tol = 1.0e-10_dp

contains

  subroutine run_all_tests(passed, failed)
    integer, intent(out) :: passed, failed
    integer :: p, f

    passed = 0
    failed = 0

    print '(A)', '================================================================'
    print '(A)', 'LEVEL 1: LINPACK EXTENDED REGRESSION TESTS'
    print '(A)', '================================================================'
    print '(A)', ''

    call test_dgbfa_dgbsl_suite(p, f)
    passed = passed + p
    failed = failed + f

    call test_dpbfa_dpbsl_suite(p, f)
    passed = passed + p
    failed = failed + f

    call test_dgeco_suite(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 1 LINPACK EXT SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! DGBFA/DGBSL Test Suite: Banded LU factorization and solve
  ! For matrices with limited bandwidth (ML subdiagonals, MU superdiagonals)
  !---------------------------------------------------------------------------
  subroutine test_dgbfa_dgbsl_suite(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'DGBFA/DGBSL (Banded LU factorization and solve)'
    print '(A)', '------------------------------------------------'

    ! Test 1: Tridiagonal system (ML=1, MU=1)
    ! [ 2 -1  0  0] [x1]   [1]
    ! [-1  2 -1  0] [x2] = [0]   x = [1, 1, 1, 1]
    ! [ 0 -1  2 -1] [x3]   [0]
    ! [ 0  0 -1  2] [x4]   [1]
    block
      integer, parameter :: n = 4, ml = 1, mu = 1
      integer, parameter :: lda = 2*ml + mu + 1  ! = 4
      real(dp) :: abd(lda, n), b(n)
      integer :: ipvt(n), info
      integer :: i

      ! Initialize banded storage to zero
      abd = 0.0_dp

      ! Store tridiagonal matrix in banded format
      ! Diagonal at row ML+MU+1 = 3
      ! Superdiagonal at row ML+MU = 2
      ! Subdiagonal at row ML+MU+2 = 4
      do i = 1, n
        abd(ml+mu+1, i) = 2.0_dp  ! diagonal
      end do
      do i = 1, n-1
        abd(ml+mu, i+1) = -1.0_dp    ! superdiagonal
        abd(ml+mu+2, i) = -1.0_dp    ! subdiagonal
      end do

      b = [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]

      call dgbfa_local(abd, lda, n, ml, mu, ipvt, info)

      if (info == 0) then
        call dgbsl_local(abd, lda, n, ml, mu, ipvt, b, 0)

        if (all(abs(b - 1.0_dp) < tol)) then
          print '(A)', '  [PASS] Tridiagonal 4x4: x = [1, 1, 1, 1]'
          passed = passed + 1
        else
          print '(A,4F8.4)', '  [FAIL] Tridiagonal 4x4: x = ', b
          failed = failed + 1
        end if
      else
        print '(A,I2)', '  [FAIL] Tridiagonal factorization failed at ', info
        failed = failed + 1
      end if
    end block

    ! Test 2: Pentadiagonal system (ML=2, MU=2)
    block
      integer, parameter :: n = 5, ml = 2, mu = 2
      integer, parameter :: lda = 2*ml + mu + 1  ! = 7
      real(dp) :: abd(lda, n), b(n), x_expected(n)
      integer :: ipvt(n), info

      abd = 0.0_dp

      ! Build a diagonally dominant pentadiagonal matrix
      ! Main diagonal = 6, first off-diagonals = -2, second off-diagonals = -1
      abd(ml+mu+1, :) = 6.0_dp        ! diagonal (row 5)
      abd(ml+mu, 2:n) = -2.0_dp       ! superdiagonal 1 (row 4)
      abd(ml+mu-1, 3:n) = -1.0_dp     ! superdiagonal 2 (row 3)
      abd(ml+mu+2, 1:n-1) = -2.0_dp   ! subdiagonal 1 (row 6)
      abd(ml+mu+3, 1:n-2) = -1.0_dp   ! subdiagonal 2 (row 7)

      ! b for x = [1, 1, 1, 1, 1]
      b = [3.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 3.0_dp]
      x_expected = 1.0_dp

      call dgbfa_local(abd, lda, n, ml, mu, ipvt, info)

      if (info == 0) then
        call dgbsl_local(abd, lda, n, ml, mu, ipvt, b, 0)

        if (all(abs(b - x_expected) < tol)) then
          print '(A)', '  [PASS] Pentadiagonal 5x5: x = [1, 1, 1, 1, 1]'
          passed = passed + 1
        else
          print '(A,5F8.4)', '  [FAIL] Pentadiagonal 5x5: x = ', b
          failed = failed + 1
        end if
      else
        print '(A,I2)', '  [FAIL] Pentadiagonal factorization failed at ', info
        failed = failed + 1
      end if
    end block

    ! Test 3: Diagonal matrix (ML=0, MU=0) - edge case
    block
      integer, parameter :: n = 3, ml = 0, mu = 0
      integer, parameter :: lda = 2*ml + mu + 1  ! = 1
      real(dp) :: abd(lda, n), b(n)
      integer :: ipvt(n), info

      abd(1, :) = [2.0_dp, 4.0_dp, 8.0_dp]
      b = [4.0_dp, 12.0_dp, 32.0_dp]  ! x = [2, 3, 4]

      call dgbfa_local(abd, lda, n, ml, mu, ipvt, info)

      if (info == 0) then
        call dgbsl_local(abd, lda, n, ml, mu, ipvt, b, 0)

        if (abs(b(1) - 2.0_dp) < tol .and. &
            abs(b(2) - 3.0_dp) < tol .and. &
            abs(b(3) - 4.0_dp) < tol) then
          print '(A)', '  [PASS] Diagonal 3x3: x = [2, 3, 4]'
          passed = passed + 1
        else
          print '(A,3F8.4)', '  [FAIL] Diagonal 3x3: x = ', b
          failed = failed + 1
        end if
      else
        print '(A)', '  [FAIL] Diagonal factorization failed'
        failed = failed + 1
      end if
    end block

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_dgbfa_dgbsl_suite

  !---------------------------------------------------------------------------
  ! DPBFA/DPBSL Test Suite: Banded Cholesky factorization and solve
  ! For symmetric positive definite banded matrices
  !---------------------------------------------------------------------------
  subroutine test_dpbfa_dpbsl_suite(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'DPBFA/DPBSL (Banded Cholesky factorization and solve)'
    print '(A)', '------------------------------------------------------'

    ! Test 1: Tridiagonal SPD (ML=MU=1, stored as upper bandwidth M=1)
    ! [ 2 -1  0  0]
    ! [-1  2 -1  0]  This is SPD (eigenvalues > 0)
    ! [ 0 -1  2 -1]
    ! [ 0  0 -1  2]
    block
      integer, parameter :: n = 4, m = 1  ! m = bandwidth (upper)
      integer, parameter :: lda = m + 1   ! = 2
      real(dp) :: abd(lda, n), b(n)
      integer :: info

      ! Store upper triangle in banded format
      ! Diagonal at row M+1 = 2
      ! Superdiagonal at row M = 1
      abd(m+1, :) = 2.0_dp           ! diagonal
      abd(m, 2:n) = -1.0_dp          ! superdiagonal

      b = [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]  ! For x = [1, 1, 1, 1]

      call dpbfa_local(abd, lda, n, m, info)

      if (info == 0) then
        call dpbsl_local(abd, lda, n, m, b)

        if (all(abs(b - 1.0_dp) < tol)) then
          print '(A)', '  [PASS] Tridiagonal SPD 4x4: x = [1, 1, 1, 1]'
          passed = passed + 1
        else
          print '(A,4F8.4)', '  [FAIL] Tridiagonal SPD 4x4: x = ', b
          failed = failed + 1
        end if
      else
        print '(A,I2)', '  [FAIL] Tridiagonal SPD Cholesky failed at ', info
        failed = failed + 1
      end if
    end block

    ! Test 2: Diagonal SPD (M=0) - edge case
    block
      integer, parameter :: n = 3, m = 0
      integer, parameter :: lda = m + 1  ! = 1
      real(dp) :: abd(lda, n), b(n)
      integer :: info

      abd(1, :) = [4.0_dp, 9.0_dp, 16.0_dp]
      b = [8.0_dp, 27.0_dp, 64.0_dp]  ! x = [2, 3, 4]

      call dpbfa_local(abd, lda, n, m, info)

      if (info == 0) then
        call dpbsl_local(abd, lda, n, m, b)

        if (abs(b(1) - 2.0_dp) < tol .and. &
            abs(b(2) - 3.0_dp) < tol .and. &
            abs(b(3) - 4.0_dp) < tol) then
          print '(A)', '  [PASS] Diagonal SPD 3x3: x = [2, 3, 4]'
          passed = passed + 1
        else
          print '(A,3F8.4)', '  [FAIL] Diagonal SPD 3x3: x = ', b
          failed = failed + 1
        end if
      else
        print '(A)', '  [FAIL] Diagonal SPD Cholesky failed'
        failed = failed + 1
      end if
    end block

    ! Test 3: Bandwidth 2 SPD
    block
      integer, parameter :: n = 4, m = 2
      integer, parameter :: lda = m + 1  ! = 3
      real(dp) :: abd(lda, n), b(n)
      integer :: info

      ! Build SPD banded matrix with bandwidth 2
      ! Main diag = 10, first super = -3, second super = 1
      abd = 0.0_dp
      abd(m+1, :) = 10.0_dp          ! diagonal (row 3)
      abd(m, 2:n) = -3.0_dp          ! superdiagonal 1 (row 2)
      abd(m-1, 3:n) = 1.0_dp         ! superdiagonal 2 (row 1)

      ! b for x = [1, 1, 1, 1]
      b = [8.0_dp, 5.0_dp, 5.0_dp, 8.0_dp]

      call dpbfa_local(abd, lda, n, m, info)

      if (info == 0) then
        call dpbsl_local(abd, lda, n, m, b)

        if (all(abs(b - 1.0_dp) < tol)) then
          print '(A)', '  [PASS] Bandwidth-2 SPD 4x4: x = [1, 1, 1, 1]'
          passed = passed + 1
        else
          print '(A,4F8.4)', '  [FAIL] Bandwidth-2 SPD 4x4: x = ', b
          failed = failed + 1
        end if
      else
        print '(A,I2)', '  [FAIL] Bandwidth-2 SPD Cholesky failed at ', info
        failed = failed + 1
      end if
    end block

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_dpbfa_dpbsl_suite

  !---------------------------------------------------------------------------
  ! DGECO Test Suite: Condition number estimation
  ! Returns RCOND = 1/(||A|| * ||A^-1||) estimate
  !---------------------------------------------------------------------------
  subroutine test_dgeco_suite(passed, failed)
    integer, intent(out) :: passed, failed
    real(dp) :: A(3,3), z(3), rcond
    integer :: ipvt(3), info

    passed = 0
    failed = 0

    print '(A)', 'DGECO (Condition number estimation)'
    print '(A)', '------------------------------------'

    ! Test 1: Identity matrix (condition number = 1, rcond = 1)
    A = 0.0_dp
    A(1,1) = 1.0_dp; A(2,2) = 1.0_dp; A(3,3) = 1.0_dp

    call dgeco_local(A, 3, 3, ipvt, rcond, z)

    if (abs(rcond - 1.0_dp) < 0.1_dp) then
      print '(A,F8.4)', '  [PASS] Identity: rcond = ', rcond
      passed = passed + 1
    else
      print '(A,F8.4)', '  [FAIL] Identity: rcond = ', rcond
      failed = failed + 1
    end if

    ! Test 2: Diagonal matrix (condition number = max/min diagonal)
    A = 0.0_dp
    A(1,1) = 1.0_dp; A(2,2) = 2.0_dp; A(3,3) = 4.0_dp
    ! Condition number = 4/1 = 4, so rcond ≈ 0.25

    call dgeco_local(A, 3, 3, ipvt, rcond, z)

    if (rcond > 0.1_dp .and. rcond < 0.5_dp) then
      print '(A,F8.4)', '  [PASS] Diagonal (cond~4): rcond = ', rcond
      passed = passed + 1
    else
      print '(A,F8.4)', '  [FAIL] Diagonal: rcond = ', rcond
      failed = failed + 1
    end if

    ! Test 3: Ill-conditioned Hilbert 3x3 (rcond should be small)
    A(1,:) = [1.0_dp, 0.5_dp, 1.0_dp/3.0_dp]
    A(2,:) = [0.5_dp, 1.0_dp/3.0_dp, 0.25_dp]
    A(3,:) = [1.0_dp/3.0_dp, 0.25_dp, 0.2_dp]

    call dgeco_local(A, 3, 3, ipvt, rcond, z)

    if (rcond < 0.01_dp .and. rcond > 0.0_dp) then
      print '(A,ES10.2)', '  [PASS] Hilbert 3x3 (ill-cond): rcond = ', rcond
      passed = passed + 1
    else
      print '(A,ES10.2)', '  [FAIL] Hilbert 3x3: rcond = ', rcond
      failed = failed + 1
    end if

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_dgeco_suite

  !---------------------------------------------------------------------------
  ! Local implementations
  !---------------------------------------------------------------------------

  ! DGBFA: LU factorization for banded matrices
  subroutine dgbfa_local(abd, lda, n, ml, mu, ipvt, info)
    integer, intent(in) :: lda, n, ml, mu
    real(dp), intent(inout) :: abd(lda, *)
    integer, intent(out) :: ipvt(*), info
    integer :: i, j, k, l, m, lm, mm, ju, jz
    real(dp) :: t

    m = ml + mu + 1
    info = 0

    ! Zero initial fill-in area
    jz = min(mu+1, n)
    if (jz > 1) then
      do j = 1, jz - 1
        do i = ml + 2 - j, ml
          abd(i, j) = 0.0_dp
        end do
      end do
    end if

    jz = jz + 1
    ju = 0

    do k = 1, n - 1
      ! Zero next fill-in column
      if (jz <= n) then
        do i = 1, ml
          abd(i, jz) = 0.0_dp
        end do
        jz = jz + 1
      end if

      ! Find pivot
      lm = min(ml, n - k)
      l = m
      do j = m + 1, m + lm
        if (abs(abd(j, k)) > abs(abd(l, k))) l = j
      end do
      ipvt(k) = l + k - m

      if (abd(l, k) == 0.0_dp) then
        info = k
        cycle
      end if

      ! Swap
      if (l /= m) then
        t = abd(l, k)
        abd(l, k) = abd(m, k)
        abd(m, k) = t
      end if

      ! Compute multipliers
      t = -1.0_dp / abd(m, k)
      do i = m + 1, m + lm
        abd(i, k) = abd(i, k) * t
      end do

      ! Row elimination
      ju = max(ju, min(mu + ipvt(k), n))
      mm = m
      do j = k + 1, ju
        l = l - 1
        mm = mm - 1
        t = abd(l, j)
        if (l /= mm) then
          abd(l, j) = abd(mm, j)
          abd(mm, j) = t
        end if
        do i = 1, lm
          abd(mm + i, j) = abd(mm + i, j) + t * abd(m + i, k)
        end do
      end do
    end do

    ipvt(n) = n
    if (abd(m, n) == 0.0_dp) info = n

  end subroutine dgbfa_local

  ! DGBSL: Solve using banded LU factorization
  subroutine dgbsl_local(abd, lda, n, ml, mu, ipvt, b, job)
    integer, intent(in) :: lda, n, ml, mu, job
    real(dp), intent(in) :: abd(lda, *)
    integer, intent(in) :: ipvt(*)
    real(dp), intent(inout) :: b(*)
    integer :: k, l, m, lm
    real(dp) :: t

    m = mu + ml + 1

    if (job == 0) then
      ! Solve A*x = b

      ! Forward substitution
      if (ml > 0) then
        do k = 1, n - 1
          lm = min(ml, n - k)
          l = ipvt(k)
          t = b(l)
          if (l /= k) then
            b(l) = b(k)
            b(k) = t
          end if
          do l = 1, lm
            b(k + l) = b(k + l) + t * abd(m + l, k)
          end do
        end do
      end if

      ! Back substitution
      do k = n, 1, -1
        b(k) = b(k) / abd(m, k)
        t = -b(k)
        lm = min(k - 1, m - 1)
        do l = 1, lm
          b(k - l) = b(k - l) + t * abd(m - l, k)
        end do
      end do
    end if

  end subroutine dgbsl_local

  ! DPBFA: Cholesky factorization for banded SPD matrices
  subroutine dpbfa_local(abd, lda, n, m, info)
    integer, intent(in) :: lda, n, m
    real(dp), intent(inout) :: abd(lda, *)
    integer, intent(out) :: info
    integer :: j, k, ik, jk, mu
    real(dp) :: s, t

    do j = 1, n
      s = 0.0_dp
      ik = m + 1
      jk = max(j - m, 1)
      mu = max(m + 2 - j, 1)

      if (m >= mu) then
        do k = mu, m
          t = abd(k, j) - sum(abd(ik:m, jk) * abd(mu:k-1, j))
          t = t / abd(m + 1, jk)
          abd(k, j) = t
          s = s + t * t
          ik = ik - 1
          jk = jk + 1
        end do
      end if

      s = abd(m + 1, j) - s
      if (s <= 0.0_dp) then
        info = j
        return
      end if
      abd(m + 1, j) = sqrt(s)
    end do
    info = 0

  end subroutine dpbfa_local

  ! DPBSL: Solve using banded Cholesky factorization
  subroutine dpbsl_local(abd, lda, n, m, b)
    integer, intent(in) :: lda, n, m
    real(dp), intent(in) :: abd(lda, *)
    real(dp), intent(inout) :: b(*)
    integer :: k, l, lm, la, lb
    real(dp) :: t

    ! Solve R'*y = b (forward substitution)
    do k = 1, n
      lm = min(k - 1, m)
      la = m + 1 - lm
      lb = k - lm
      t = sum(abd(la:m, k) * b(lb:k-1))
      b(k) = (b(k) - t) / abd(m + 1, k)
    end do

    ! Solve R*x = y (back substitution)
    ! R(k, k+l) is stored at ABD(M+1-l, k+l)
    do k = n, 1, -1
      lm = min(m, n - k)
      t = 0.0_dp
      do l = 1, lm
        t = t + abd(m + 1 - l, k + l) * b(k + l)
      end do
      b(k) = (b(k) - t) / abd(m + 1, k)
    end do

  end subroutine dpbsl_local

  ! DGECO: Factor and estimate condition number
  subroutine dgeco_local(a, lda, n, ipvt, rcond, z)
    integer, intent(in) :: lda, n
    real(dp), intent(inout) :: a(lda, *)
    integer, intent(out) :: ipvt(*)
    real(dp), intent(out) :: rcond, z(*)
    real(dp) :: anorm, ynorm, s, sm, t, wk, wkm, ek
    integer :: i, j, k, l, info

    ! Compute 1-norm of A
    anorm = 0.0_dp
    do j = 1, n
      s = 0.0_dp
      do i = 1, n
        s = s + abs(a(i, j))
      end do
      anorm = max(anorm, s)
    end do

    ! Factor
    call dgefa_internal(a, lda, n, ipvt, info)

    ! Estimate norm of inverse
    ek = 1.0_dp
    z(1:n) = 0.0_dp

    do k = 1, n
      if (z(k) /= 0.0_dp) ek = sign(ek, -z(k))
      if (abs(ek - z(k)) > abs(a(k, k))) then
        s = abs(a(k, k)) / abs(ek - z(k))
        z(1:n) = s * z(1:n)
        ek = s * ek
      end if
      wk = ek - z(k)
      wkm = -ek - z(k)
      s = abs(wk)
      sm = abs(wkm)
      if (a(k, k) /= 0.0_dp) then
        wk = wk / a(k, k)
        wkm = wkm / a(k, k)
      else
        wk = 1.0_dp
        wkm = 1.0_dp
      end if
      if (k + 1 <= n) then
        do j = k + 1, n
          sm = sm + abs(z(j) + wkm * a(k, j))
          z(j) = z(j) + wk * a(k, j)
          s = s + abs(z(j))
        end do
        if (s < sm) then
          t = wkm - wk
          wk = wkm
          do j = k + 1, n
            z(j) = z(j) + t * a(k, j)
          end do
        end if
      end if
      z(k) = wk
    end do

    s = sum(abs(z(1:n)))
    z(1:n) = z(1:n) / s

    ! Solve L'*y = z
    do k = n, 1, -1
      if (k < n) then
        z(k) = z(k) + sum(a(k+1:n, k) * z(k+1:n))
      end if
      if (abs(z(k)) > 1.0_dp) then
        s = 1.0_dp / abs(z(k))
        z(1:n) = s * z(1:n)
      end if
      l = ipvt(k)
      t = z(l)
      z(l) = z(k)
      z(k) = t
    end do

    s = sum(abs(z(1:n)))
    z(1:n) = z(1:n) / s
    ynorm = 1.0_dp

    ! Solve L*v = y
    do k = 1, n
      l = ipvt(k)
      t = z(l)
      z(l) = z(k)
      z(k) = t
      if (k < n) then
        do i = k + 1, n
          z(i) = z(i) + t * a(i, k)
        end do
      end if
      if (abs(z(k)) > 1.0_dp) then
        s = 1.0_dp / abs(z(k))
        z(1:n) = s * z(1:n)
        ynorm = s * ynorm
      end if
    end do

    s = sum(abs(z(1:n)))
    z(1:n) = z(1:n) / s
    ynorm = ynorm / s

    ! Solve U*z = v
    do k = n, 1, -1
      if (abs(z(k)) > abs(a(k, k))) then
        s = abs(a(k, k)) / abs(z(k))
        z(1:n) = s * z(1:n)
        ynorm = s * ynorm
      end if
      if (a(k, k) /= 0.0_dp) then
        z(k) = z(k) / a(k, k)
      else
        z(k) = 1.0_dp
      end if
      t = -z(k)
      do i = 1, k - 1
        z(i) = z(i) + t * a(i, k)
      end do
    end do

    s = sum(abs(z(1:n)))
    z(1:n) = z(1:n) / s
    ynorm = ynorm / s

    if (anorm /= 0.0_dp) then
      rcond = ynorm / anorm
    else
      rcond = 0.0_dp
    end if

  contains
    subroutine dgefa_internal(a, lda, n, ipvt, info)
      integer, intent(in) :: lda, n
      real(dp), intent(inout) :: a(lda, *)
      integer, intent(out) :: ipvt(*), info
      integer :: ii, jj, kk, ll
      real(dp) :: tt

      info = 0
      do kk = 1, n - 1
        ll = kk
        do jj = kk + 1, n
          if (abs(a(jj, kk)) > abs(a(ll, kk))) ll = jj
        end do
        ipvt(kk) = ll
        if (a(ll, kk) == 0.0_dp) then
          info = kk
          cycle
        end if
        if (ll /= kk) then
          tt = a(ll, kk); a(ll, kk) = a(kk, kk); a(kk, kk) = tt
        end if
        tt = -1.0_dp / a(kk, kk)
        do jj = kk + 1, n
          a(jj, kk) = a(jj, kk) * tt
        end do
        do jj = kk + 1, n
          tt = a(ll, jj)
          if (ll /= kk) then
            a(ll, jj) = a(kk, jj); a(kk, jj) = tt
          end if
          do ii = kk + 1, n
            a(ii, jj) = a(ii, jj) + tt * a(ii, kk)
          end do
        end do
      end do
      ipvt(n) = n
      if (a(n, n) == 0.0_dp) info = n
    end subroutine
  end subroutine dgeco_local

end module test_linpack_extended_level1

program run_level1_linpack_extended
  use test_linpack_extended_level1
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level1_linpack_extended

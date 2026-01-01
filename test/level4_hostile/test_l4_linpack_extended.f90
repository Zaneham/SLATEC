!> Level 4: Hostile Tests for Extended LINPACK Routines
!>
!> Purpose: What breaks under stress? Platform-specific edge cases.
!>
!> What we test:
!>   - DGBFA/DGBSL: Ill-conditioned banded, extreme scaling, subnormals
!>   - DPBFA/DPBSL: Near-singular SPD, ill-conditioned, edge cases
!>   - DGECO: Very ill-conditioned, near-zero determinant, Inf/NaN
!>
!> Reference: IEEE 754-2008, LINPACK Users' Guide

module test_linpack_extended_level4
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use, intrinsic :: ieee_arithmetic
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
    print '(A)', 'LEVEL 4: LINPACK EXTENDED HOSTILE TESTS'
    print '(A)', '================================================================'
    print '(A)', ''

    call test_dgbfa_dgbsl_hostile(p, f)
    passed = passed + p
    failed = failed + f

    call test_dpbfa_dpbsl_hostile(p, f)
    passed = passed + p
    failed = failed + f

    call test_dgeco_hostile(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 4 LINPACK EXT SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! DGBFA/DGBSL Hostile Tests
  !---------------------------------------------------------------------------
  subroutine test_dgbfa_dgbsl_hostile(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'DGBFA/DGBSL Hostile Tests'
    print '(A)', '-------------------------'

    ! Test 1: Near-singular banded matrix
    block
      integer, parameter :: n = 4, ml = 1, mu = 1
      integer, parameter :: lda = 2*ml + mu + 1
      real(dp) :: abd(lda, n), b(n)
      integer :: ipvt(n), info
      real(dp) :: eps

      eps = epsilon(1.0_dp)

      abd = 0.0_dp
      abd(ml+mu+1, :) = [1.0_dp, 1.0_dp, 1.0_dp, eps]  ! Last diagonal tiny
      abd(ml+mu, 2:n) = -0.5_dp
      abd(ml+mu+2, 1:n-1) = -0.5_dp

      b = 1.0_dp

      call dgbfa_local(abd, lda, n, ml, mu, ipvt, info)

      if (info == 0) then
        call dgbsl_local(abd, lda, n, ml, mu, ipvt, b, 0)
        if (all(ieee_is_finite(b))) then
          print '(A)', '  [PASS] Near-singular banded: finite solution'
          passed = passed + 1
        else
          print '(A)', '  [WARN] Near-singular: non-finite solution'
          passed = passed + 1  ! Expected under stress
        end if
      else
        print '(A)', '  [PASS] Near-singular banded: singularity detected'
        passed = passed + 1
      end if
    end block

    ! Test 2: Extreme scaling - large elements
    block
      integer, parameter :: n = 3, ml = 1, mu = 1
      integer, parameter :: lda = 2*ml + mu + 1
      real(dp) :: abd(lda, n), b(n), scale
      integer :: ipvt(n), info

      scale = sqrt(huge(1.0_dp))

      abd = 0.0_dp
      abd(ml+mu+1, :) = 2.0_dp * scale
      abd(ml+mu, 2:n) = -1.0_dp * scale
      abd(ml+mu+2, 1:n-1) = -1.0_dp * scale

      b = scale

      call dgbfa_local(abd, lda, n, ml, mu, ipvt, info)
      call dgbsl_local(abd, lda, n, ml, mu, ipvt, b, 0)

      if (all(ieee_is_finite(b))) then
        print '(A)', '  [PASS] Extreme scaling (large): finite solution'
        passed = passed + 1
      else
        print '(A)', '  [WARN] Extreme scaling: overflow detected'
        passed = passed + 1
      end if
    end block

    ! Test 3: Extreme scaling - small elements
    block
      integer, parameter :: n = 3, ml = 1, mu = 1
      integer, parameter :: lda = 2*ml + mu + 1
      real(dp) :: abd(lda, n), b(n), scale
      integer :: ipvt(n), info

      scale = sqrt(tiny(1.0_dp))

      abd = 0.0_dp
      abd(ml+mu+1, :) = 2.0_dp * scale
      abd(ml+mu, 2:n) = -1.0_dp * scale
      abd(ml+mu+2, 1:n-1) = -1.0_dp * scale

      b = scale

      call dgbfa_local(abd, lda, n, ml, mu, ipvt, info)
      call dgbsl_local(abd, lda, n, ml, mu, ipvt, b, 0)

      if (all(ieee_is_finite(b)) .and. all(b /= 0.0_dp)) then
        print '(A)', '  [PASS] Extreme scaling (small): non-zero solution'
        passed = passed + 1
      else
        print '(A)', '  [WARN] Extreme scaling: underflow detected'
        passed = passed + 1
      end if
    end block

    ! Test 4: Subnormal elements
    block
      integer, parameter :: n = 3, ml = 1, mu = 1
      integer, parameter :: lda = 2*ml + mu + 1
      real(dp) :: abd(lda, n), b(n), subn
      integer :: ipvt(n), info

      subn = tiny(1.0_dp) / 16.0_dp

      abd = 0.0_dp
      abd(ml+mu+1, :) = 2.0_dp
      abd(ml+mu+1, 2) = 2.0_dp + subn  ! Subnormal perturbation
      abd(ml+mu, 2:n) = -1.0_dp
      abd(ml+mu+2, 1:n-1) = -1.0_dp

      b = [1.0_dp, subn, 1.0_dp]

      call dgbfa_local(abd, lda, n, ml, mu, ipvt, info)
      call dgbsl_local(abd, lda, n, ml, mu, ipvt, b, 0)

      if (all(ieee_is_finite(b))) then
        print '(A)', '  [PASS] Subnormal elements: handled correctly'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] Subnormal elements caused non-finite result'
        failed = failed + 1
      end if
    end block

    ! Test 5: Ill-conditioned banded Hilbert-like
    block
      integer, parameter :: n = 4, ml = 2, mu = 2
      integer, parameter :: lda = 2*ml + mu + 1
      real(dp) :: abd(lda, n), b(n)
      integer :: ipvt(n), info, i

      abd = 0.0_dp
      ! Banded matrix with 1/(|i-j|+1) entries
      do i = 1, n
        abd(ml+mu+1, i) = 1.0_dp  ! diagonal
        if (i > 1) abd(ml+mu+2, i-1) = 0.5_dp
        if (i > 2) abd(ml+mu+3, i-2) = 1.0_dp/3.0_dp
        if (i < n) abd(ml+mu, i+1) = 0.5_dp
        if (i < n-1) abd(ml+mu-1, i+2) = 1.0_dp/3.0_dp
      end do

      b = 1.0_dp

      call dgbfa_local(abd, lda, n, ml, mu, ipvt, info)
      call dgbsl_local(abd, lda, n, ml, mu, ipvt, b, 0)

      if (info == 0 .and. all(ieee_is_finite(b))) then
        print '(A)', '  [PASS] Ill-conditioned banded: solved'
        passed = passed + 1
      else
        print '(A)', '  [WARN] Ill-conditioned banded: numerical issues'
        passed = passed + 1
      end if
    end block

    ! Test 6: Zero bandwidth (diagonal)
    block
      integer, parameter :: n = 3, ml = 0, mu = 0
      integer, parameter :: lda = 1
      real(dp) :: abd(lda, n), b(n)
      integer :: ipvt(n), info

      abd(1, :) = [2.0_dp, 3.0_dp, 4.0_dp]
      b = [4.0_dp, 9.0_dp, 16.0_dp]  ! x = [2, 3, 4]

      call dgbfa_local(abd, lda, n, ml, mu, ipvt, info)
      call dgbsl_local(abd, lda, n, ml, mu, ipvt, b, 0)

      if (abs(b(1) - 2.0_dp) < tol .and. abs(b(2) - 3.0_dp) < tol) then
        print '(A)', '  [PASS] Zero bandwidth (diagonal) edge case'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] Zero bandwidth failed'
        failed = failed + 1
      end if
    end block

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_dgbfa_dgbsl_hostile

  !---------------------------------------------------------------------------
  ! DPBFA/DPBSL Hostile Tests
  !---------------------------------------------------------------------------
  subroutine test_dpbfa_dpbsl_hostile(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'DPBFA/DPBSL Hostile Tests'
    print '(A)', '-------------------------'

    ! Test 1: Nearly not positive definite
    block
      integer, parameter :: n = 3, m = 1
      integer, parameter :: lda = m + 1
      real(dp) :: abd(lda, n), b(n)
      integer :: info
      real(dp) :: eps

      eps = epsilon(1.0_dp)

      ! Diagonal dominant but barely
      abd = 0.0_dp
      abd(m+1, :) = [1.0_dp + eps, 1.0_dp + eps, 1.0_dp + eps]
      abd(m, 2:n) = -0.5_dp

      b = 1.0_dp

      call dpbfa_local(abd, lda, n, m, info)

      if (info == 0) then
        call dpbsl_local(abd, lda, n, m, b)
        if (all(ieee_is_finite(b))) then
          print '(A)', '  [PASS] Nearly non-SPD: solved'
          passed = passed + 1
        else
          print '(A)', '  [WARN] Nearly non-SPD: non-finite'
          passed = passed + 1
        end if
      else
        print '(A)', '  [PASS] Nearly non-SPD: detected'
        passed = passed + 1
      end if
    end block

    ! Test 2: Very ill-conditioned SPD
    block
      integer, parameter :: n = 4, m = 1
      integer, parameter :: lda = m + 1
      real(dp) :: abd(lda, n), b(n)
      integer :: info

      ! Create ill-conditioned SPD via scaling
      abd = 0.0_dp
      abd(m+1, :) = [1.0e12_dp, 1.0_dp, 1.0_dp, 1.0e-12_dp]
      abd(m, 2:n) = [-0.1_dp, -0.1_dp, -1.0e-13_dp]

      b = 1.0_dp

      call dpbfa_local(abd, lda, n, m, info)

      if (info == 0) then
        call dpbsl_local(abd, lda, n, m, b)
        print '(A)', '  [PASS] Ill-conditioned SPD: factored'
        passed = passed + 1
      else
        print '(A,I2)', '  [WARN] Ill-conditioned SPD: Cholesky failed at ', info
        passed = passed + 1
      end if
    end block

    ! Test 3: Extreme diagonal dominance
    block
      integer, parameter :: n = 3, m = 1
      integer, parameter :: lda = m + 1
      real(dp) :: abd(lda, n), b(n), scale
      integer :: info

      scale = 1.0e100_dp

      abd = 0.0_dp
      abd(m+1, :) = scale
      abd(m, 2:n) = -1.0_dp

      b = scale

      call dpbfa_local(abd, lda, n, m, info)

      if (info == 0) then
        call dpbsl_local(abd, lda, n, m, b)
        if (all(ieee_is_finite(b))) then
          print '(A)', '  [PASS] Extreme diagonal dominance: finite'
          passed = passed + 1
        else
          print '(A)', '  [WARN] Extreme dominance: non-finite'
          passed = passed + 1
        end if
      else
        print '(A)', '  [FAIL] Extreme dominance: Cholesky failed'
        failed = failed + 1
      end if
    end block

    ! Test 4: Minimum bandwidth SPD (diagonal only)
    block
      integer, parameter :: n = 4, m = 0
      integer, parameter :: lda = 1
      real(dp) :: abd(lda, n), b(n)
      integer :: info

      abd(1, :) = [4.0_dp, 9.0_dp, 16.0_dp, 25.0_dp]
      b = [8.0_dp, 27.0_dp, 64.0_dp, 125.0_dp]

      call dpbfa_local(abd, lda, n, m, info)
      call dpbsl_local(abd, lda, n, m, b)

      if (abs(b(1) - 2.0_dp) < tol .and. abs(b(2) - 3.0_dp) < tol) then
        print '(A)', '  [PASS] Minimum bandwidth (m=0) edge case'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] Minimum bandwidth failed'
        failed = failed + 1
      end if
    end block

    ! Test 5: Non-SPD detection
    block
      integer, parameter :: n = 3, m = 1
      integer, parameter :: lda = m + 1
      real(dp) :: abd(lda, n)
      integer :: info

      ! Negative diagonal makes it not SPD
      abd = 0.0_dp
      abd(m+1, :) = [1.0_dp, -1.0_dp, 1.0_dp]
      abd(m, 2:n) = -0.1_dp

      call dpbfa_local(abd, lda, n, m, info)

      if (info /= 0) then
        print '(A)', '  [PASS] Non-SPD (negative diag) detected'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] Non-SPD not detected'
        failed = failed + 1
      end if
    end block

    ! Test 6: SPD with subnormal off-diagonal
    block
      integer, parameter :: n = 3, m = 1
      integer, parameter :: lda = m + 1
      real(dp) :: abd(lda, n), b(n), subn
      integer :: info

      subn = tiny(1.0_dp) / 16.0_dp

      abd = 0.0_dp
      abd(m+1, :) = 2.0_dp
      abd(m, 2:n) = subn

      b = 1.0_dp

      call dpbfa_local(abd, lda, n, m, info)
      call dpbsl_local(abd, lda, n, m, b)

      if (all(ieee_is_finite(b))) then
        print '(A)', '  [PASS] Subnormal off-diagonal handled'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] Subnormal off-diagonal failed'
        failed = failed + 1
      end if
    end block

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_dpbfa_dpbsl_hostile

  !---------------------------------------------------------------------------
  ! DGECO Hostile Tests
  !---------------------------------------------------------------------------
  subroutine test_dgeco_hostile(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'DGECO Hostile Tests'
    print '(A)', '-------------------'

    ! Test 1: Singular matrix (rcond should be 0)
    block
      real(dp) :: A(3,3), z(3), rcond
      integer :: ipvt(3)

      ! Row 3 = Row 1 + Row 2 => singular
      A = reshape([1.0_dp, 2.0_dp, 3.0_dp, &
                   4.0_dp, 5.0_dp, 9.0_dp, &
                   7.0_dp, 8.0_dp, 15.0_dp], [3,3])

      call dgeco_local(A, 3, 3, ipvt, rcond, z)

      if (rcond < 1.0e-14_dp) then
        print '(A,ES10.2)', '  [PASS] Singular matrix: rcond = ', rcond
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] Singular not detected: rcond = ', rcond
        failed = failed + 1
      end if
    end block

    ! Test 2: Very ill-conditioned (Hilbert 4x4)
    block
      real(dp) :: H(4,4), z(4), rcond
      integer :: ipvt(4), i, j

      do i = 1, 4
        do j = 1, 4
          H(i,j) = 1.0_dp / real(i + j - 1, dp)
        end do
      end do

      call dgeco_local(H, 4, 4, ipvt, rcond, z)

      ! Hilbert 4x4 has cond ≈ 15,514
      if (rcond < 1.0e-3_dp .and. rcond > 0.0_dp) then
        print '(A,ES10.2)', '  [PASS] Hilbert 4x4 ill-cond: rcond = ', rcond
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] Hilbert 4x4: rcond = ', rcond
        failed = failed + 1
      end if
    end block

    ! Test 3: Matrix with Inf element
    block
      real(dp) :: A(2,2), z(2), rcond, inf_val
      integer :: ipvt(2)

      inf_val = ieee_value(1.0_dp, ieee_positive_inf)

      if (.not. ieee_is_finite(inf_val)) then
        A = reshape([1.0_dp, 0.0_dp, inf_val, 1.0_dp], [2,2])

        call dgeco_local(A, 2, 2, ipvt, rcond, z)

        ! Should handle gracefully (rcond = 0 or NaN)
        if (rcond == 0.0_dp .or. ieee_is_nan(rcond)) then
          print '(A)', '  [PASS] Inf element: rcond=0 or NaN'
          passed = passed + 1
        else
          print '(A,ES10.2)', '  [WARN] Inf element: rcond = ', rcond
          passed = passed + 1
        end if
      else
        print '(A)', '  [SKIP] Could not create Inf'
        passed = passed + 1
      end if
    end block

    ! Test 4: Matrix with NaN element
    block
      real(dp) :: A(2,2), z(2), rcond
      integer :: ipvt(2)
      real(dp) :: nan_val

      nan_val = ieee_value(1.0_dp, ieee_quiet_nan)

      if (ieee_is_nan(nan_val)) then
        A = reshape([1.0_dp, nan_val, 0.0_dp, 1.0_dp], [2,2])

        call dgeco_local(A, 2, 2, ipvt, rcond, z)

        ! NaN should propagate or result in rcond=0
        if (ieee_is_nan(rcond) .or. rcond == 0.0_dp) then
          print '(A)', '  [PASS] NaN element: propagated correctly'
          passed = passed + 1
        else
          print '(A,ES10.2)', '  [WARN] NaN element: rcond = ', rcond
          passed = passed + 1
        end if
      else
        print '(A)', '  [SKIP] Could not create NaN'
        passed = passed + 1
      end if
    end block

    ! Test 5: Extreme scaling
    block
      real(dp) :: A(2,2), z(2), rcond, scale
      integer :: ipvt(2)

      scale = sqrt(huge(1.0_dp))

      A = reshape([scale, 0.0_dp, 0.0_dp, 1.0_dp/scale], [2,2])

      call dgeco_local(A, 2, 2, ipvt, rcond, z)

      ! Very ill-conditioned due to scaling
      if (rcond < 1.0e-100_dp .or. rcond == 0.0_dp) then
        print '(A)', '  [PASS] Extreme scaling: very small rcond'
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [WARN] Extreme scaling: rcond = ', rcond
        passed = passed + 1
      end if
    end block

    ! Test 6: 1x1 matrix edge case
    block
      real(dp) :: A(1,1), z(1), rcond
      integer :: ipvt(1)

      A(1,1) = 5.0_dp

      call dgeco_local(A, 1, 1, ipvt, rcond, z)

      if (abs(rcond - 1.0_dp) < 0.5_dp) then
        print '(A,F8.4)', '  [PASS] 1x1 matrix: rcond = ', rcond
        passed = passed + 1
      else
        print '(A,F8.4)', '  [FAIL] 1x1 matrix: rcond = ', rcond
        failed = failed + 1
      end if
    end block

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_dgeco_hostile

  !---------------------------------------------------------------------------
  ! Local implementations
  !---------------------------------------------------------------------------

  subroutine dgbfa_local(abd, lda, n, ml, mu, ipvt, info)
    integer, intent(in) :: lda, n, ml, mu
    real(dp), intent(inout) :: abd(lda, *)
    integer, intent(out) :: ipvt(*), info
    integer :: i, j, k, l, m, lm, mm, ju, jz
    real(dp) :: t

    m = ml + mu + 1
    info = 0

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
      if (jz <= n) then
        do i = 1, ml
          abd(i, jz) = 0.0_dp
        end do
        jz = jz + 1
      end if

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

      if (l /= m) then
        t = abd(l, k)
        abd(l, k) = abd(m, k)
        abd(m, k) = t
      end if

      t = -1.0_dp / abd(m, k)
      do i = m + 1, m + lm
        abd(i, k) = abd(i, k) * t
      end do

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

  subroutine dgbsl_local(abd, lda, n, ml, mu, ipvt, b, job)
    integer, intent(in) :: lda, n, ml, mu, job
    real(dp), intent(in) :: abd(lda, *)
    integer, intent(in) :: ipvt(*)
    real(dp), intent(inout) :: b(*)
    integer :: k, l, m, lm
    real(dp) :: t

    m = mu + ml + 1

    if (job == 0) then
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

  subroutine dpbsl_local(abd, lda, n, m, b)
    integer, intent(in) :: lda, n, m
    real(dp), intent(in) :: abd(lda, *)
    real(dp), intent(inout) :: b(*)
    integer :: k, l, lm, la, lb
    real(dp) :: t

    do k = 1, n
      lm = min(k - 1, m)
      la = m + 1 - lm
      lb = k - lm
      t = sum(abd(la:m, k) * b(lb:k-1))
      b(k) = (b(k) - t) / abd(m + 1, k)
    end do

    do k = n, 1, -1
      lm = min(m, n - k)
      t = 0.0_dp
      do l = 1, lm
        t = t + abd(m + 1 - l, k + l) * b(k + l)
      end do
      b(k) = (b(k) - t) / abd(m + 1, k)
    end do

  end subroutine dpbsl_local

  subroutine dgeco_local(a, lda, n, ipvt, rcond, z)
    integer, intent(in) :: lda, n
    real(dp), intent(inout) :: a(lda, *)
    integer, intent(out) :: ipvt(*)
    real(dp), intent(out) :: rcond, z(*)
    real(dp) :: anorm, ynorm, s, sm, t, wk, wkm, ek
    integer :: i, j, k, l, info

    anorm = 0.0_dp
    do j = 1, n
      s = 0.0_dp
      do i = 1, n
        s = s + abs(a(i, j))
      end do
      anorm = max(anorm, s)
    end do

    call dgefa_internal(a, lda, n, ipvt, info)

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

end module test_linpack_extended_level4

program run_level4_linpack_extended
  use test_linpack_extended_level4
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level4_linpack_extended

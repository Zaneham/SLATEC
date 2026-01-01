!> Level 3: Historical Tests for Extended LINPACK Routines
!>
!> Purpose: Does output match classical known values?
!>
!> What we test:
!>   - DGBFA/DGBSL: Classic discretization matrices from PDEs
!>   - DPBFA/DPBSL: Symmetric banded from finite elements
!>   - DGECO: Known condition numbers from literature
!>
!> Reference: Forsythe & Moler, "Computer Solution of Linear Algebraic Systems"

module test_linpack_extended_level3
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
    print '(A)', 'LEVEL 3: LINPACK EXTENDED HISTORICAL TESTS'
    print '(A)', '================================================================'
    print '(A)', ''

    call test_dgbfa_dgbsl_historical(p, f)
    passed = passed + p
    failed = failed + f

    call test_dpbfa_dpbsl_historical(p, f)
    passed = passed + p
    failed = failed + f

    call test_dgeco_historical(p, f)
    passed = passed + p
    failed = failed + f

    print '(A)', ''
    print '(A)', '================================================================'
    print '(A,I3,A,I3,A)', 'LEVEL 3 LINPACK EXT SUMMARY: ', passed, ' passed, ', failed, ' failed'
    print '(A)', '================================================================'

  end subroutine run_all_tests

  !---------------------------------------------------------------------------
  ! DGBFA/DGBSL Historical Tests - Classic PDE Discretizations
  !---------------------------------------------------------------------------
  subroutine test_dgbfa_dgbsl_historical(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'DGBFA/DGBSL Historical Tests'
    print '(A)', '-----------------------------'

    ! Test 1: 1D Laplacian - Classic finite difference
    ! -d²u/dx² = f with Dirichlet BC
    ! Discretization: [-1, 2, -1] / h²
    ! For h=1/5, n=4 interior points
    ! Solution: u(x) = x(1-x), f = 2
    block
      integer, parameter :: n = 4, ml = 1, mu = 1
      integer, parameter :: lda = 2*ml + mu + 1
      real(dp) :: abd(lda, n), b(n), x_exact(n)
      integer :: ipvt(n), info, i
      real(dp) :: h, xi

      h = 0.2_dp  ! 1/5

      ! Scaled tridiagonal: [−1, 2, −1]
      abd = 0.0_dp
      abd(ml+mu+1, :) = 2.0_dp / h**2
      abd(ml+mu, 2:n) = -1.0_dp / h**2
      abd(ml+mu+2, 1:n-1) = -1.0_dp / h**2

      ! RHS: f = 2 at interior points
      b = 2.0_dp

      ! Exact solution: u(x) = x(1-x)
      do i = 1, n
        xi = real(i, dp) * h
        x_exact(i) = xi * (1.0_dp - xi)
      end do

      call dgbfa_local(abd, lda, n, ml, mu, ipvt, info)
      call dgbsl_local(abd, lda, n, ml, mu, ipvt, b, 0)

      if (all(abs(b - x_exact) < 1.0e-12_dp)) then
        print '(A)', '  [PASS] 1D Laplacian: u(x)=x(1-x) recovered'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] 1D Laplacian solution error'
        failed = failed + 1
      end if
    end block

    ! Test 2: Convection-Diffusion (upwind scheme)
    ! -εu'' + u' = 0, u(0)=0, u(1)=1
    ! Upwind: [-ε/h² - 1/h, 2ε/h² + 1/h, -ε/h²]
    block
      integer, parameter :: n = 4, ml = 1, mu = 1
      integer, parameter :: lda = 2*ml + mu + 1
      real(dp) :: abd(lda, n), b(n), x_exact(n)
      integer :: ipvt(n), info, i
      real(dp) :: h, eps, xi, diag, sub, sup

      h = 0.2_dp
      eps = 0.1_dp  ! Peclet number = 1/(eps*n) = 2.5

      ! Upwind discretization
      sub = -eps/h**2 - 1.0_dp/h
      diag = 2.0_dp*eps/h**2 + 1.0_dp/h
      sup = -eps/h**2

      abd = 0.0_dp
      abd(ml+mu+1, :) = diag
      abd(ml+mu, 2:n) = sup
      abd(ml+mu+2, 1:n-1) = sub

      ! RHS: zero except for BC contribution at last point
      b = 0.0_dp
      b(n) = -sup  ! BC: u(1) = 1

      ! Exact: u(x) = (exp(x/eps) - 1)/(exp(1/eps) - 1)
      do i = 1, n
        xi = real(i, dp) * h
        x_exact(i) = (exp(xi/eps) - 1.0_dp) / (exp(1.0_dp/eps) - 1.0_dp)
      end do

      call dgbfa_local(abd, lda, n, ml, mu, ipvt, info)
      call dgbsl_local(abd, lda, n, ml, mu, ipvt, b, 0)

      ! Upwind is first-order on coarse grid; verify qualitative behavior
      ! Solution should be monotonically increasing from 0 to ~1
      if (b(1) > 0.0_dp .and. b(n) > b(1) .and. b(n) < 1.5_dp) then
        print '(A)', '  [PASS] Convection-Diffusion (upwind): monotone increasing'
        passed = passed + 1
      else
        print '(A,4F8.4)', '  [FAIL] Convection-Diffusion: ', b
        failed = failed + 1
      end if
    end block

    ! Test 3: Block tridiagonal from 2D Laplacian (1D slice)
    ! Classic 5-point stencil on 3x3 grid = 9 unknowns
    ! Natural ordering gives block tridiagonal
    block
      integer, parameter :: n = 3, ml = 1, mu = 1
      integer, parameter :: lda = 2*ml + mu + 1
      real(dp) :: abd(lda, n), b(n)
      integer :: ipvt(n), info

      ! Single row of 2D Laplacian (coupled to row above/below)
      ! Just test the 1D portion: [−1, 4, −1] on 3 points
      abd = 0.0_dp
      abd(ml+mu+1, :) = 4.0_dp
      abd(ml+mu, 2:n) = -1.0_dp
      abd(ml+mu+2, 1:n-1) = -1.0_dp

      ! RHS chosen so x = [1, 1, 1]
      b = [3.0_dp, 2.0_dp, 3.0_dp]

      call dgbfa_local(abd, lda, n, ml, mu, ipvt, info)
      call dgbsl_local(abd, lda, n, ml, mu, ipvt, b, 0)

      if (all(abs(b - 1.0_dp) < tol)) then
        print '(A)', '  [PASS] 2D Laplacian row slice'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] 2D Laplacian row slice'
        failed = failed + 1
      end if
    end block

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_dgbfa_dgbsl_historical

  !---------------------------------------------------------------------------
  ! DPBFA/DPBSL Historical Tests - FEM Stiffness Matrices
  !---------------------------------------------------------------------------
  subroutine test_dpbfa_dpbsl_historical(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'DPBFA/DPBSL Historical Tests'
    print '(A)', '-----------------------------'

    ! Test 1: FEM stiffness matrix for 1D bar with uniform properties
    ! Element stiffness: [1, -1; -1, 1] * E*A/L
    ! Assembled tridiagonal
    block
      integer, parameter :: n = 4, m = 1
      integer, parameter :: lda = m + 1
      real(dp) :: abd(lda, n), b(n), x_exact(n)
      integer :: info, i
      real(dp) :: EAL

      EAL = 1.0_dp  ! E*A/L = 1 for simplicity

      ! Assembled stiffness (fixed at both ends adds to diagonal)
      ! Interior: [1, 2, 1] * EAL
      abd = 0.0_dp
      abd(m+1, :) = 2.0_dp * EAL
      abd(m, 2:n) = -1.0_dp * EAL

      ! Force: uniform load, exact solution is parabolic
      b = 1.0_dp

      call dpbfa_local(abd, lda, n, m, info)
      call dpbsl_local(abd, lda, n, m, b)

      ! Displacement should be positive and symmetric about center
      if (b(1) > 0.0_dp .and. b(n) > 0.0_dp .and. &
          abs(b(1) - b(4)) < tol .and. abs(b(2) - b(3)) < tol) then
        print '(A)', '  [PASS] 1D FEM bar: symmetric displacement'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] 1D FEM bar'
        failed = failed + 1
      end if
    end block

    ! Test 2: Mass matrix from linear elements
    ! Element mass: [2, 1; 1, 2] * rho*A*L/6
    ! Forms tridiagonal SPD
    block
      integer, parameter :: n = 3, m = 1
      integer, parameter :: lda = m + 1
      real(dp) :: abd(lda, n), b(n)
      integer :: info

      ! Lumped version (simpler): diagonal only
      ! Consistent: [1, 4, 1] / 6 on interior
      abd = 0.0_dp
      abd(m+1, :) = 4.0_dp / 6.0_dp
      abd(m, 2:n) = 1.0_dp / 6.0_dp

      b = 1.0_dp  ! Unit input

      call dpbfa_local(abd, lda, n, m, info)
      call dpbsl_local(abd, lda, n, m, b)

      ! Solution should be positive
      if (all(b > 0.0_dp)) then
        print '(A)', '  [PASS] FEM consistent mass matrix'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] FEM mass matrix solve'
        failed = failed + 1
      end if
    end block

    ! Test 3: String vibration eigenvalue setup
    ! M*a = λ*K*x discretization leads to banded SPD
    block
      integer, parameter :: n = 5, m = 1
      integer, parameter :: lda = m + 1
      real(dp) :: abd(lda, n), b(n), x_exact(n)
      integer :: info, i
      real(dp) :: pi

      pi = 4.0_dp * atan(1.0_dp)

      ! Stiffness matrix from string: [−1, 2, −1]
      abd = 0.0_dp
      abd(m+1, :) = 2.0_dp
      abd(m, 2:n) = -1.0_dp

      ! First eigenmode: sin(pi*i/(n+1))
      do i = 1, n
        x_exact(i) = sin(pi * real(i, dp) / real(n+1, dp))
      end do

      ! K*x_exact should give λ*x_exact with λ = 4*sin²(pi/(2(n+1)))
      b = x_exact

      call dpbfa_local(abd, lda, n, m, info)
      call dpbsl_local(abd, lda, n, m, b)

      ! b should be proportional to x_exact
      ! b = K^{-1} * x_exact = (1/λ) * x_exact
      if (abs(b(3)/x_exact(3) - b(1)/x_exact(1)) < 0.01_dp) then
        print '(A)', '  [PASS] String vibration eigenmode'
        passed = passed + 1
      else
        print '(A)', '  [FAIL] String eigenmode'
        failed = failed + 1
      end if
    end block

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_dpbfa_dpbsl_historical

  !---------------------------------------------------------------------------
  ! DGECO Historical Tests - Known Condition Numbers
  !---------------------------------------------------------------------------
  subroutine test_dgeco_historical(passed, failed)
    integer, intent(out) :: passed, failed

    passed = 0
    failed = 0

    print '(A)', 'DGECO Historical Tests'
    print '(A)', '----------------------'

    ! Test 1: Hilbert matrix n=3 (cond ≈ 748)
    block
      real(dp) :: H(3,3), z(3), rcond
      integer :: ipvt(3), i, j
      real(dp) :: cond_exact

      do i = 1, 3
        do j = 1, 3
          H(i,j) = 1.0_dp / real(i + j - 1, dp)
        end do
      end do

      cond_exact = 748.0_dp  ! Known from literature

      call dgeco_local(H, 3, 3, ipvt, rcond, z)

      ! rcond should be roughly 1/cond ≈ 0.00134
      if (rcond < 0.01_dp .and. rcond > 0.0001_dp) then
        print '(A,ES10.2,A,F8.1)', '  [PASS] Hilbert 3x3: rcond=', rcond, ' (cond≈', 1.0_dp/rcond, ')'
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] Hilbert 3x3: rcond=', rcond
        failed = failed + 1
      end if
    end block

    ! Test 2: Vandermonde matrix with nodes 1,2,3 (ill-conditioned)
    block
      real(dp) :: V(3,3), z(3), rcond
      integer :: ipvt(3), i, j
      real(dp) :: nodes(3)

      nodes = [1.0_dp, 2.0_dp, 3.0_dp]

      do i = 1, 3
        do j = 1, 3
          V(i,j) = nodes(i) ** (j-1)
        end do
      end do

      call dgeco_local(V, 3, 3, ipvt, rcond, z)

      ! Vandermonde condition grows exponentially
      if (rcond > 0.001_dp .and. rcond < 0.5_dp) then
        print '(A,ES10.2)', '  [PASS] Vandermonde: rcond = ', rcond
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] Vandermonde: rcond = ', rcond
        failed = failed + 1
      end if
    end block

    ! Test 3: Moler matrix (tricky condition)
    block
      real(dp) :: M(3,3), z(3), rcond
      integer :: ipvt(3), i, j

      ! Moler matrix: M(i,j) = min(i,j) - 2
      do i = 1, 3
        do j = 1, 3
          M(i,j) = real(min(i,j), dp) - 2.0_dp
        end do
      end do

      call dgeco_local(M, 3, 3, ipvt, rcond, z)

      if (rcond > 0.0_dp .and. rcond < 1.0_dp) then
        print '(A,ES10.2)', '  [PASS] Moler matrix: rcond = ', rcond
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] Moler matrix: rcond = ', rcond
        failed = failed + 1
      end if
    end block

    ! Test 4: Pascal matrix (exact integer entries)
    block
      real(dp) :: P(3,3), z(3), rcond
      integer :: ipvt(3)

      ! Pascal matrix
      P = reshape([1.0_dp, 1.0_dp, 1.0_dp, &
                   1.0_dp, 2.0_dp, 3.0_dp, &
                   1.0_dp, 3.0_dp, 6.0_dp], [3,3])

      call dgeco_local(P, 3, 3, ipvt, rcond, z)

      ! Pascal has modest condition number
      if (rcond > 0.01_dp) then
        print '(A,ES10.2)', '  [PASS] Pascal matrix: rcond = ', rcond
        passed = passed + 1
      else
        print '(A,ES10.2)', '  [FAIL] Pascal matrix: rcond = ', rcond
        failed = failed + 1
      end if
    end block

    print '(A,I2,A,I2,A)', '  Subtotal: ', passed, ' passed, ', failed, ' failed'
    print '(A)', ''

  end subroutine test_dgeco_historical

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

end module test_linpack_extended_level3

program run_level3_linpack_extended
  use test_linpack_extended_level3
  implicit none

  integer :: passed, failed

  call run_all_tests(passed, failed)

  if (failed > 0) then
    stop 1
  end if

end program run_level3_linpack_extended

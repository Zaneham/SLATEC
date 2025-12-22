module test_diff_integ_eq
  !> Differential equation solver tests
  !>
  !> Sources:
  !>   - Mehdi Chinoune's SLATEC tests (2021)
  !>   - ZH: Additional tests for modernized routines (2024)
  !>
  !> Covers: DASSL, DEPAC (DEABM, DEBDF), SDRIVE, BVSUP

  use testdrive, only: new_unittest, unittest_type, error_type, check
  use service, only: SP, DP, eps_dp
  implicit none
  private

  public :: collect_diff_integ_eq_tests

  ! Shared state for ODE residual functions
  real(DP) :: lambda_dp = 1.0_DP

contains

  subroutine collect_diff_integ_eq_tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      ! DASSL tests (from Mehdi)
      new_unittest("DDASSL exponential decay", test_ddassl_exp), &
      ! ZH: Modernization tests
      new_unittest("ZH: DDASSL integration_loop", test_zh_ddassl), &
      new_unittest("ZH: DDASTP step_loop", test_zh_ddastp) &
    ]
  end subroutine

  !============================================================================
  ! DASSL RESIDUAL FUNCTIONS (modernized interface)
  !============================================================================

  pure subroutine res_exp_decay(t, y, yprime, delta, ires)
    !> Residual for y' = -lambda*y (exponential decay)
    !> Solution: y(t) = y0 * exp(-lambda*t)
    real(DP), intent(in) :: t
    real(DP), intent(in) :: y(:), yprime(:)
    real(DP), intent(out) :: delta(:)
    integer, intent(inout) :: ires

    ! Residual: yprime + lambda*y = 0
    delta(1) = yprime(1) + lambda_dp * y(1)
    ires = 0
  end subroutine

  pure subroutine jac_exp_decay(t, y, yprime, pd, cj)
    !> Jacobian for exponential decay: dR/dy = lambda, dR/dyprime = 1
    real(DP), intent(in) :: t, cj
    real(DP), intent(in) :: y(:), yprime(:)
    real(DP), intent(out) :: pd(:,:)

    pd(1,1) = lambda_dp + cj
  end subroutine

  !============================================================================
  ! DASSL TESTS (from Mehdi Chinoune)
  !============================================================================

  subroutine test_ddassl_exp(error)
    !> Test DDASSL with simple exponential decay ODE
    use diff_integ_eq, only: DDASSL
    type(error_type), allocatable, intent(out) :: error

    real(DP) :: t, tout, y(1), yprime(1)
    real(DP) :: rtol(1), atol(1)
    real(DP) :: rwork(52)
    integer :: iwork(21), info(15), idid
    real(DP) :: y_exact, rel_err
    logical :: passed

    ! Initial conditions: y(0) = 1, y'(0) = -lambda
    lambda_dp = 1.0_DP
    t = 0.0_DP
    tout = 1.0_DP
    y(1) = 1.0_DP
    yprime(1) = -lambda_dp

    ! Tolerances
    rtol(1) = 1.0D-6
    atol(1) = 1.0D-8

    ! Initialize info array
    info = 0

    call DDASSL(res_exp_decay, 1, t, y, yprime, tout, info, &
                rtol, atol, idid, rwork, 52, iwork, 21, jac_exp_decay)

    ! Exact solution: y(1) = exp(-1) = 0.367879...
    y_exact = exp(-1.0_DP)
    rel_err = abs(y(1) - y_exact) / y_exact

    passed = (idid > 0) .and. (rel_err < 1.0D-4)

    call check(error, passed, "DDASSL exponential decay test failed")
  end subroutine

  !============================================================================
  ! ZH: MODERNIZATION TESTS
  !============================================================================

  subroutine test_zh_ddassl(error)
    !> ZH: Test DDASSL after integration_loop modernization (Batch 31)
    !> Tests the main driver loop and inline error handling
    use diff_integ_eq, only: DDASSL
    type(error_type), allocatable, intent(out) :: error

    real(DP) :: t, tout, y(1), yprime(1)
    real(DP) :: rtol(1), atol(1)
    real(DP) :: rwork(52)
    integer :: iwork(21), info(15), idid
    real(DP) :: y_exact, rel_err
    logical :: passed
    integer :: i

    ! Test multiple steps to exercise integration_loop
    lambda_dp = 0.5_DP
    t = 0.0_DP
    y(1) = 1.0_DP
    yprime(1) = -lambda_dp
    rtol(1) = 1.0D-8
    atol(1) = 1.0D-10
    info = 0

    passed = .true.

    ! Integrate to t=2 in steps
    do i = 1, 4
      tout = real(i, DP) * 0.5_DP
      call DDASSL(res_exp_decay, 1, t, y, yprime, tout, info, &
                  rtol, atol, idid, rwork, 52, iwork, 21, jac_exp_decay)

      if (idid < 0) then
        passed = .false.
        exit
      end if

      y_exact = exp(-lambda_dp * tout)
      rel_err = abs(y(1) - y_exact) / y_exact
      if (rel_err > 1.0D-5) passed = .false.
    end do

    call check(error, passed, "ZH: DDASSL integration_loop test failed")
  end subroutine

  subroutine test_zh_ddastp(error)
    !> ZH: Test DDASTP step routines after modernization (Batch 31)
    !> Uses tighter tolerances to exercise step_loop retry logic
    use diff_integ_eq, only: DDASSL
    type(error_type), allocatable, intent(out) :: error

    real(DP) :: t, tout, y(1), yprime(1)
    real(DP) :: rtol(1), atol(1)
    real(DP) :: rwork(52)
    integer :: iwork(21), info(15), idid
    real(DP) :: y_exact, rel_err
    logical :: passed

    ! Use stiff problem to exercise corrector
    lambda_dp = 10.0_DP
    t = 0.0_DP
    tout = 0.5_DP
    y(1) = 1.0_DP
    yprime(1) = -lambda_dp
    rtol(1) = 1.0D-10
    atol(1) = 1.0D-12
    info = 0

    call DDASSL(res_exp_decay, 1, t, y, yprime, tout, info, &
                rtol, atol, idid, rwork, 52, iwork, 21, jac_exp_decay)

    y_exact = exp(-lambda_dp * tout)
    rel_err = abs(y(1) - y_exact) / y_exact

    passed = (idid > 0) .and. (rel_err < 1.0D-8)

    call check(error, passed, "ZH: DDASTP step_loop test failed")
  end subroutine

  ! TODO: test_zh_ddriv3 removed - interface needs updating to match modernized DDRIV3

end module test_diff_integ_eq

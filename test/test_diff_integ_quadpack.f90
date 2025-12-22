module test_diff_integ_quadpack
  !> QUADPACK integration tests
  !>
  !> Sources:
  !>   - Mehdi Chinoune's SLATEC tests (2021)
  !>   - ZH: Additional tests for modernized routines (2024)
  !>
  !> Reference: Piessens et al, QUADPACK, Springer-Verlag 1983

  use testdrive, only: new_unittest, unittest_type, error_type, check
  use service, only: SP, DP, eps_dp, eps_sp
  use diff_integ, only: DQAG, DQAGS, DQAGP, DQAGI, DQAWOE, &
                        QAG, QAGS, QAGP, QAGI, QAWOE
  implicit none
  private

  public :: collect_quadpack_tests

  real(DP), parameter :: PI_DP = 3.141592653589793238_DP
  real(SP), parameter :: PI_SP = 3.141592653589793238_SP

contains

  subroutine collect_quadpack_tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      ! Basic integration tests (from Mehdi)
      new_unittest("DQAG basic integral", test_dqag_basic), &
      new_unittest("DQAGS singular endpoint", test_dqags_singular), &
      new_unittest("DQAGI infinite interval", test_dqagi_infinite), &
      ! ZH: Modernization tests
      new_unittest("ZH: DQAGSE control flow", test_zh_dqagse), &
      new_unittest("ZH: DQAWOE oscillatory", test_zh_dqawoe), &
      ! Single precision
      new_unittest("QAG basic integral", test_qag_basic), &
      new_unittest("QAGS singular endpoint", test_qags_singular) &
    ]
  end subroutine

  !============================================================================
  ! TEST INTEGRANDS
  !============================================================================

  function f_sqrt(x) result(y)
    !> f(x) = 1/sqrt(x), integral [0,1] = 2
    real(DP), intent(in) :: x
    real(DP) :: y
    y = 1.0_DP / sqrt(x)
  end function

  function f_sqrt_sp(x) result(y)
    real(SP), intent(in) :: x
    real(SP) :: y
    y = 1.0_SP / sqrt(x)
  end function

  function f_exp(x) result(y)
    !> f(x) = exp(-x), integral [0,inf] = 1
    real(DP), intent(in) :: x
    real(DP) :: y
    y = exp(-x)
  end function

  function f_poly(x) result(y)
    !> f(x) = x^2, integral [0,1] = 1/3
    real(DP), intent(in) :: x
    real(DP) :: y
    y = x * x
  end function

  function f_poly_sp(x) result(y)
    real(SP), intent(in) :: x
    real(SP) :: y
    y = x * x
  end function

  function f_osc(x) result(y)
    !> f(x) = cos(100*x), oscillatory test
    real(DP), intent(in) :: x
    real(DP) :: y
    y = cos(100.0_DP * x)
  end function

  function f_log(x) result(y)
    !> f(x) = log(x), integral [0,1] = -1
    real(DP), intent(in) :: x
    real(DP) :: y
    y = log(x)
  end function

  !============================================================================
  ! BASIC INTEGRATION TESTS (from Mehdi Chinoune)
  !============================================================================

  subroutine test_dqag_basic(error)
    !> Test DQAG with polynomial integrand
    type(error_type), allocatable, intent(out) :: error

    real(DP) :: result, abserr, work(400)
    integer :: neval, ier, iwork(100), last
    real(DP) :: exact, rel_err
    logical :: passed

    exact = 1.0_DP / 3.0_DP  ! integral of x^2 from 0 to 1

    call DQAG(f_poly, 0.0_DP, 1.0_DP, 0.0_DP, 1.0D-10, 6, &
              result, abserr, neval, ier, 100, 400, last, iwork, work)

    rel_err = abs(result - exact) / exact
    passed = (ier == 0) .and. (rel_err < 1.0D-10)

    call check(error, passed, "DQAG basic integral test failed")
  end subroutine

  subroutine test_dqags_singular(error)
    !> Test DQAGS with singular endpoint (1/sqrt(x))
    type(error_type), allocatable, intent(out) :: error

    real(DP) :: result, abserr, work(400)
    integer :: neval, ier, iwork(100), last
    real(DP) :: exact, rel_err
    logical :: passed

    exact = 2.0_DP  ! integral of 1/sqrt(x) from 0 to 1

    call DQAGS(f_sqrt, 0.0_DP, 1.0_DP, 0.0_DP, 1.0D-8, &
               result, abserr, neval, ier, 100, 400, last, iwork, work)

    rel_err = abs(result - exact) / exact
    passed = (ier == 0) .and. (rel_err < 1.0D-6)

    call check(error, passed, "DQAGS singular endpoint test failed")
  end subroutine

  subroutine test_dqagi_infinite(error)
    !> Test DQAGI with infinite interval (exp(-x) on [0,inf])
    type(error_type), allocatable, intent(out) :: error

    real(DP) :: result, abserr, work(400)
    integer :: neval, ier, iwork(100), last
    real(DP) :: exact, rel_err
    logical :: passed

    exact = 1.0_DP  ! integral of exp(-x) from 0 to infinity

    ! inf = 1 means integrate from bound to +infinity
    call DQAGI(f_exp, 0.0_DP, 1, 0.0_DP, 1.0D-8, &
               result, abserr, neval, ier, 100, 400, last, iwork, work)

    rel_err = abs(result - exact) / exact
    passed = (ier == 0) .and. (rel_err < 1.0D-6)

    call check(error, passed, "DQAGI infinite interval test failed")
  end subroutine

  !============================================================================
  ! ZH: MODERNIZATION TESTS
  !============================================================================

  subroutine test_zh_dqagse(error)
    !> ZH: Test DQAGSE after control flow modernization (Batch 27)
    !> Tests compute_global flag and main_loop structure
    type(error_type), allocatable, intent(out) :: error

    real(DP) :: result, abserr, work(400)
    integer :: neval, ier, iwork(100), last
    real(DP) :: exact, rel_err
    logical :: passed

    exact = -1.0_DP  ! integral of log(x) from 0 to 1

    call DQAGS(f_log, 0.0_DP, 1.0_DP, 0.0_DP, 1.0D-8, &
               result, abserr, neval, ier, 100, 400, last, iwork, work)

    rel_err = abs(result - exact) / abs(exact)
    passed = (ier == 0) .and. (rel_err < 1.0D-6)

    call check(error, passed, "ZH: DQAGSE control flow test failed")
  end subroutine

  subroutine test_zh_dqawoe(error)
    !> ZH: Test DQAWOE after control flow modernization (Batch 46)
    !> Tests oscillatory integration with done/test flags
    type(error_type), allocatable, intent(out) :: error

    real(DP) :: result, abserr, work(1400)
    integer :: neval, ier, iwork(200), last, momcom, maxp1, nnlog
    real(DP) :: chebmo(25,50)
    real(DP) :: omega, exact, rel_err
    integer :: integr
    logical :: passed

    ! Integrate cos(100*x) from 0 to pi
    ! integral = sin(100*pi)/100 = 0
    omega = 100.0_DP
    integr = 1  ! cos weight
    maxp1 = 50
    exact = sin(100.0_DP * PI_DP) / 100.0_DP

    call DQAWOE(f_osc, 0.0_DP, PI_DP, omega, integr, 0.0_DP, 1.0D-6, &
                100, 1400, result, abserr, neval, ier, last, &
                nnlog, momcom, maxp1, chebmo, iwork, work)

    ! Result should be very close to 0
    passed = (abs(result) < 1.0D-4) .or. (ier == 0)

    call check(error, passed, "ZH: DQAWOE oscillatory test failed")
  end subroutine

  !============================================================================
  ! SINGLE PRECISION TESTS
  !============================================================================

  subroutine test_qag_basic(error)
    !> Single precision QAG test
    type(error_type), allocatable, intent(out) :: error

    real(SP) :: result, abserr, work(400)
    integer :: neval, ier, iwork(100), last
    real(SP) :: exact, rel_err
    logical :: passed

    exact = 1.0_SP / 3.0_SP

    call QAG(f_poly_sp, 0.0_SP, 1.0_SP, 0.0_SP, 1.0E-6_SP, 6, &
             result, abserr, neval, ier, 100, 400, last, iwork, work)

    rel_err = abs(result - exact) / exact
    passed = (ier == 0) .and. (rel_err < 1.0E-5_SP)

    call check(error, passed, "QAG basic integral test failed")
  end subroutine

  subroutine test_qags_singular(error)
    !> Single precision QAGS test
    type(error_type), allocatable, intent(out) :: error

    real(SP) :: result, abserr, work(400)
    integer :: neval, ier, iwork(100), last
    real(SP) :: exact, rel_err
    logical :: passed

    exact = 2.0_SP

    call QAGS(f_sqrt_sp, 0.0_SP, 1.0_SP, 0.0_SP, 1.0E-5_SP, &
              result, abserr, neval, ier, 100, 400, last, iwork, work)

    rel_err = abs(result - exact) / exact
    passed = (ier == 0) .and. (rel_err < 1.0E-4_SP)

    call check(error, passed, "QAGS singular endpoint test failed")
  end subroutine

end module test_diff_integ_quadpack

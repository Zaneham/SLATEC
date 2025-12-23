module test_diff_integ_quadpack
  !> QUADPACK integration tests - Golden regression tests
  !>
  !> Sources:
  !>   - Mehdi Chinoune's SLATEC tests (2021)
  !>   - ZH: Golden tests from Piessens et al, QUADPACK, Springer-Verlag 1983
  !>
  !> All test values are exact analytical results, not library outputs.

  use testdrive, only: new_unittest, unittest_type, error_type, check
  use service, only: SP, DP
  use diff_integ, only: DQAG, DQAGS, DQAGI, QAG, QAGS, QAGI
  implicit none
  private

  public :: collect_quadpack_tests

contains

  subroutine collect_quadpack_tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      ! Mehdi: Basic integration tests
      new_unittest("DQAG x^2 [0,1] = 1/3", test_dqag_poly), &
      new_unittest("DQAGS 1/sqrt(x) [0,1] = 2", test_dqags_sqrt), &
      new_unittest("DQAGI exp(-x) [0,inf] = 1", test_dqagi_exp), &
      ! ZH: Golden tests from QUADPACK book Table 5.1
      new_unittest("ZH: DQAGS log(x) [0,1] = -1", test_zh_log), &
      ! Single precision
      new_unittest("QAG x^2 [0,1] = 1/3 (SP)", test_qag_poly_sp), &
      new_unittest("QAGS 1/sqrt(x) [0,1] = 2 (SP)", test_qags_sqrt_sp) &
    ]
  end subroutine

  !============================================================================
  ! TEST INTEGRANDS (all PURE for modernized interfaces)
  !============================================================================

  pure function f_sqrt_dp(x) result(y)
    !> f(x) = 1/sqrt(x), integral [0,1] = 2 (exact)
    real(DP), intent(in) :: x
    real(DP) :: y
    y = 1.0_DP / sqrt(x)
  end function

  pure function f_sqrt_sp(x) result(y)
    real(SP), intent(in) :: x
    real(SP) :: y
    y = 1.0_SP / sqrt(x)
  end function

  pure function f_exp_dp(x) result(y)
    !> f(x) = exp(-x), integral [0,inf] = 1 (exact)
    real(DP), intent(in) :: x
    real(DP) :: y
    y = exp(-x)
  end function

  pure function f_poly_dp(x) result(y)
    !> f(x) = x^2, integral [0,1] = 1/3 (exact)
    real(DP), intent(in) :: x
    real(DP) :: y
    y = x * x
  end function

  pure function f_poly_sp(x) result(y)
    real(SP), intent(in) :: x
    real(SP) :: y
    y = x * x
  end function

  pure function f_log_dp(x) result(y)
    !> f(x) = log(x), integral [0,1] = -1 (exact)
    !> Reference: QUADPACK Table 5.1, Problem 1
    real(DP), intent(in) :: x
    real(DP) :: y
    y = log(x)
  end function

  !============================================================================
  ! GOLDEN REGRESSION TESTS - Double Precision
  ! All expected values are exact analytical results
  !============================================================================

  subroutine test_dqag_poly(error)
    !> Mehdi: integral of x^2 from 0 to 1 = 1/3 (exact)
    type(error_type), allocatable, intent(out) :: error

    real(DP) :: result, abserr, work(400)
    integer :: neval, ier, iwork(100), last
    real(DP), parameter :: EXACT = 1.0_DP / 3.0_DP
    real(DP) :: rel_err

    call DQAG(f_poly_dp, 0.0_DP, 1.0_DP, 0.0_DP, 1.0D-10, 6, &
              result, abserr, neval, ier, 100, 400, last, iwork, work)

    rel_err = abs(result - EXACT) / EXACT
    call check(error, ier == 0 .and. rel_err < 1.0D-10, &
               "DQAG: integral of x^2 should equal 1/3")
  end subroutine

  subroutine test_dqags_sqrt(error)
    !> Mehdi: integral of 1/sqrt(x) from 0 to 1 = 2 (exact)
    !> Tests handling of integrable singularity at x=0
    type(error_type), allocatable, intent(out) :: error

    real(DP) :: result, abserr, work(400)
    integer :: neval, ier, iwork(100), last
    real(DP), parameter :: EXACT = 2.0_DP
    real(DP) :: rel_err

    call DQAGS(f_sqrt_dp, 0.0_DP, 1.0_DP, 0.0_DP, 1.0D-8, &
               result, abserr, neval, ier, 100, 400, last, iwork, work)

    rel_err = abs(result - EXACT) / EXACT
    call check(error, ier == 0 .and. rel_err < 1.0D-6, &
               "DQAGS: integral of 1/sqrt(x) should equal 2")
  end subroutine

  subroutine test_dqagi_exp(error)
    !> Mehdi: integral of exp(-x) from 0 to infinity = 1 (exact)
    type(error_type), allocatable, intent(out) :: error

    real(DP) :: result, abserr, work(400)
    integer :: neval, ier, iwork(100), last
    real(DP), parameter :: EXACT = 1.0_DP
    real(DP) :: rel_err

    call DQAGI(f_exp_dp, 0.0_DP, 1, 0.0_DP, 1.0D-8, &
               result, abserr, neval, ier, 100, 400, last, iwork, work)

    rel_err = abs(result - EXACT) / EXACT
    call check(error, ier == 0 .and. rel_err < 1.0D-6, &
               "DQAGI: integral of exp(-x) [0,inf] should equal 1")
  end subroutine

  subroutine test_zh_log(error)
    !> ZH: integral of log(x) from 0 to 1 = -1 (exact)
    !> Reference: QUADPACK book Table 5.1, Problem 1
    !> Tests DQAGSE control flow after GOTO elimination (Batch 27)
    type(error_type), allocatable, intent(out) :: error

    real(DP) :: result, abserr, work(400)
    integer :: neval, ier, iwork(100), last
    real(DP), parameter :: EXACT = -1.0_DP
    real(DP) :: rel_err

    call DQAGS(f_log_dp, 0.0_DP, 1.0_DP, 0.0_DP, 1.0D-8, &
               result, abserr, neval, ier, 100, 400, last, iwork, work)

    rel_err = abs(result - EXACT) / abs(EXACT)
    call check(error, ier == 0 .and. rel_err < 1.0D-6, &
               "DQAGS: integral of log(x) [0,1] should equal -1")
  end subroutine

  !============================================================================
  ! GOLDEN REGRESSION TESTS - Single Precision
  !============================================================================

  subroutine test_qag_poly_sp(error)
    !> Single precision: integral of x^2 from 0 to 1 = 1/3
    type(error_type), allocatable, intent(out) :: error

    real(SP) :: result, abserr, work(400)
    integer :: neval, ier, iwork(100), last
    real(SP), parameter :: EXACT = 1.0_SP / 3.0_SP
    real(SP) :: rel_err

    call QAG(f_poly_sp, 0.0_SP, 1.0_SP, 0.0_SP, 1.0E-6_SP, 6, &
             result, abserr, neval, ier, 100, 400, last, iwork, work)

    rel_err = abs(result - EXACT) / EXACT
    call check(error, ier == 0 .and. rel_err < 1.0E-5_SP, &
               "QAG (SP): integral of x^2 should equal 1/3")
  end subroutine

  subroutine test_qags_sqrt_sp(error)
    !> Single precision: integral of 1/sqrt(x) from 0 to 1 = 2
    type(error_type), allocatable, intent(out) :: error

    real(SP) :: result, abserr, work(400)
    integer :: neval, ier, iwork(100), last
    real(SP), parameter :: EXACT = 2.0_SP
    real(SP) :: rel_err

    call QAGS(f_sqrt_sp, 0.0_SP, 1.0_SP, 0.0_SP, 1.0E-5_SP, &
              result, abserr, neval, ier, 100, 400, last, iwork, work)

    rel_err = abs(result - EXACT) / EXACT
    call check(error, ier == 0 .and. rel_err < 1.0E-4_SP, &
               "QAGS (SP): integral of 1/sqrt(x) should equal 2")
  end subroutine

end module test_diff_integ_quadpack

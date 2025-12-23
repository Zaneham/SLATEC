program quadpack_deviation
  !> ZH: QUADPACK floating-point deviation tests
  !>
  !> This program tests QUADPACK integration routines against known exact values
  !> and documents any deviations across different compiler optimisation levels.
  !>
  !> Reference integrals:
  !>   - ∫₀¹ x² dx = 1/3 (exactly)
  !>   - ∫₀^π sin(x) dx = 2 (exactly)
  !>   - ∫₀¹ exp(x) dx = e - 1
  !>   - ∫₀¹ 1/(1+x²) dx = π/4 (arctan)
  !>   - ∫₀^∞ exp(-x²) dx = √π/2 (Gaussian, via DQAGI)

  use service, only: SP, DP
  use diff_integ, only: DQAGS, QAGS, DQAGI, QAGI, DQNG, QNG
  implicit none

  real(DP), parameter :: PI_DP = 3.14159265358979323846264338327950288_DP
  real(SP), parameter :: PI_SP = 3.14159265358979323846264338327950288_SP
  real(DP), parameter :: E_DP = 2.71828182845904523536028747135266250_DP
  real(SP), parameter :: E_SP = 2.71828182845904523536028747135266250_SP

  print *, "========================================"
  print *, "QUADPACK Floating-Point Deviation Report"
  print *, "========================================"
  print *

  print *, "Compiler: gfortran (assumed)"
  print *

  call test_polynomial_integral()
  call test_trig_integral()
  call test_exp_integral()
  call test_arctan_integral()
  call test_gaussian_integral()
  call test_qng_integrals()

contains

  !> ∫₀¹ x² dx = 1/3 exactly
  subroutine test_polynomial_integral()
    real(DP) :: result_dp, abserr_dp, expected_dp
    real(SP) :: result_sp, abserr_sp, expected_sp
    real(DP) :: work_dp(400)
    real(SP) :: work_sp(400)
    integer :: iwork(100), neval, ier, last

    expected_dp = 1.0_DP / 3.0_DP
    expected_sp = 1.0_SP / 3.0_SP

    print *, "=== Polynomial: ∫₀¹ x² dx = 1/3 ==="
    print *

    ! Double precision
    call DQAGS(f_x2_dp, 0.0_DP, 1.0_DP, 1.0D-12, 1.0D-12, &
               result_dp, abserr_dp, neval, ier, 100, 400, last, iwork, work_dp)

    call print_result("x² (DP)", expected_dp, result_dp, abserr_dp, neval, ier)

    ! Single precision (tolerance must be >= 50*machine_epsilon ~ 6E-6 for SP)
    call QAGS(f_x2_sp, 0.0_SP, 1.0_SP, 1.0E-5_SP, 1.0E-5_SP, &
              result_sp, abserr_sp, neval, ier, 100, 400, last, iwork, work_sp)

    call print_result_sp("x² (SP)", expected_sp, result_sp, abserr_sp, neval, ier)
    print *
  end subroutine

  !> ∫₀^π sin(x) dx = 2 exactly
  subroutine test_trig_integral()
    real(DP) :: result_dp, abserr_dp, expected_dp
    real(SP) :: result_sp, abserr_sp, expected_sp
    real(DP) :: work_dp(400)
    real(SP) :: work_sp(400)
    integer :: iwork(100), neval, ier, last

    expected_dp = 2.0_DP
    expected_sp = 2.0_SP

    print *, "=== Trigonometric: ∫₀^π sin(x) dx = 2 ==="
    print *

    ! Double precision
    call DQAGS(f_sin_dp, 0.0_DP, PI_DP, 1.0D-12, 1.0D-12, &
               result_dp, abserr_dp, neval, ier, 100, 400, last, iwork, work_dp)

    call print_result("sin(x) (DP)", expected_dp, result_dp, abserr_dp, neval, ier)

    ! Single precision (tolerance must be >= 50*machine_epsilon ~ 6E-6 for SP)
    call QAGS(f_sin_sp, 0.0_SP, PI_SP, 1.0E-5_SP, 1.0E-5_SP, &
              result_sp, abserr_sp, neval, ier, 100, 400, last, iwork, work_sp)

    call print_result_sp("sin(x) (SP)", expected_sp, result_sp, abserr_sp, neval, ier)
    print *
  end subroutine

  !> ∫₀¹ exp(x) dx = e - 1
  subroutine test_exp_integral()
    real(DP) :: result_dp, abserr_dp, expected_dp
    real(SP) :: result_sp, abserr_sp, expected_sp
    real(DP) :: work_dp(400)
    real(SP) :: work_sp(400)
    integer :: iwork(100), neval, ier, last

    expected_dp = E_DP - 1.0_DP
    expected_sp = E_SP - 1.0_SP

    print *, "=== Exponential: ∫₀¹ exp(x) dx = e - 1 ==="
    print *

    ! Double precision
    call DQAGS(f_exp_dp, 0.0_DP, 1.0_DP, 1.0D-12, 1.0D-12, &
               result_dp, abserr_dp, neval, ier, 100, 400, last, iwork, work_dp)

    call print_result("exp(x) (DP)", expected_dp, result_dp, abserr_dp, neval, ier)

    ! Single precision (tolerance must be >= 50*machine_epsilon ~ 6E-6 for SP)
    call QAGS(f_exp_sp, 0.0_SP, 1.0_SP, 1.0E-5_SP, 1.0E-5_SP, &
              result_sp, abserr_sp, neval, ier, 100, 400, last, iwork, work_sp)

    call print_result_sp("exp(x) (SP)", expected_sp, result_sp, abserr_sp, neval, ier)
    print *
  end subroutine

  !> ∫₀¹ 1/(1+x²) dx = π/4
  subroutine test_arctan_integral()
    real(DP) :: result_dp, abserr_dp, expected_dp
    real(SP) :: result_sp, abserr_sp, expected_sp
    real(DP) :: work_dp(400)
    real(SP) :: work_sp(400)
    integer :: iwork(100), neval, ier, last

    expected_dp = PI_DP / 4.0_DP
    expected_sp = PI_SP / 4.0_SP

    print *, "=== Arctan: ∫₀¹ 1/(1+x²) dx = π/4 ==="
    print *

    ! Double precision
    call DQAGS(f_arctan_dp, 0.0_DP, 1.0_DP, 1.0D-12, 1.0D-12, &
               result_dp, abserr_dp, neval, ier, 100, 400, last, iwork, work_dp)

    call print_result("1/(1+x²) (DP)", expected_dp, result_dp, abserr_dp, neval, ier)

    ! Single precision
    call QAGS(f_arctan_sp, 0.0_SP, 1.0_SP, 1.0E-5_SP, 1.0E-5_SP, &
              result_sp, abserr_sp, neval, ier, 100, 400, last, iwork, work_sp)

    call print_result_sp("1/(1+x²) (SP)", expected_sp, result_sp, abserr_sp, neval, ier)
    print *
  end subroutine

  !> ∫₀^∞ exp(-x²) dx = √π/2 (Gaussian integral, tests DQAGI)
  subroutine test_gaussian_integral()
    real(DP) :: result_dp, abserr_dp, expected_dp
    real(SP) :: result_sp, abserr_sp, expected_sp
    real(DP) :: work_dp(400)
    real(SP) :: work_sp(400)
    integer :: iwork(100), neval, ier, last

    expected_dp = sqrt(PI_DP) / 2.0_DP
    expected_sp = sqrt(PI_SP) / 2.0_SP

    print *, "=== Gaussian: ∫₀^∞ exp(-x²) dx = √π/2 ==="
    print *

    ! Double precision (inf=1 means integrate from bound to +∞)
    call DQAGI(f_gaussian_dp, 0.0_DP, 1, 1.0D-12, 1.0D-12, &
               result_dp, abserr_dp, neval, ier, 100, 400, last, iwork, work_dp)

    call print_result("exp(-x²) (DP)", expected_dp, result_dp, abserr_dp, neval, ier)

    ! Single precision
    call QAGI(f_gaussian_sp, 0.0_SP, 1, 1.0E-5_SP, 1.0E-5_SP, &
              result_sp, abserr_sp, neval, ier, 100, 400, last, iwork, work_sp)

    call print_result_sp("exp(-x²) (SP)", expected_sp, result_sp, abserr_sp, neval, ier)
    print *
  end subroutine

  !> Test QNG (non-adaptive Gauss-Kronrod) for comparison
  subroutine test_qng_integrals()
    real(DP) :: result_dp, abserr_dp, expected_dp
    real(SP) :: result_sp, abserr_sp, expected_sp
    integer :: neval, ier

    expected_dp = 1.0_DP / 3.0_DP
    expected_sp = 1.0_SP / 3.0_SP

    print *, "=== QNG (Non-adaptive): ∫₀¹ x² dx = 1/3 ==="
    print *

    ! Double precision
    call DQNG(f_x2_dp, 0.0_DP, 1.0_DP, 1.0D-12, 1.0D-12, &
              result_dp, abserr_dp, neval, ier)

    call print_result("QNG x² (DP)", expected_dp, result_dp, abserr_dp, neval, ier)

    ! Single precision
    call QNG(f_x2_sp, 0.0_SP, 1.0_SP, 1.0E-5_SP, 1.0E-5_SP, &
             result_sp, abserr_sp, neval, ier)

    call print_result_sp("QNG x² (SP)", expected_sp, result_sp, abserr_sp, neval, ier)
    print *
  end subroutine

  !--- Integrand functions (must be PURE) ---

  pure real(DP) function f_x2_dp(x)
    real(DP), intent(in) :: x
    f_x2_dp = x * x
  end function

  pure real(SP) function f_x2_sp(x)
    real(SP), intent(in) :: x
    f_x2_sp = x * x
  end function

  pure real(DP) function f_sin_dp(x)
    real(DP), intent(in) :: x
    f_sin_dp = sin(x)
  end function

  pure real(SP) function f_sin_sp(x)
    real(SP), intent(in) :: x
    f_sin_sp = sin(x)
  end function

  pure real(DP) function f_exp_dp(x)
    real(DP), intent(in) :: x
    f_exp_dp = exp(x)
  end function

  pure real(SP) function f_exp_sp(x)
    real(SP), intent(in) :: x
    f_exp_sp = exp(x)
  end function

  pure real(DP) function f_arctan_dp(x)
    real(DP), intent(in) :: x
    f_arctan_dp = 1.0_DP / (1.0_DP + x * x)
  end function

  pure real(SP) function f_arctan_sp(x)
    real(SP), intent(in) :: x
    f_arctan_sp = 1.0_SP / (1.0_SP + x * x)
  end function

  pure real(DP) function f_gaussian_dp(x)
    real(DP), intent(in) :: x
    f_gaussian_dp = exp(-x * x)
  end function

  pure real(SP) function f_gaussian_sp(x)
    real(SP), intent(in) :: x
    f_gaussian_sp = exp(-x * x)
  end function

  !--- Output routines ---

  subroutine print_result(name, expected, computed, abserr, neval, ier)
    character(*), intent(in) :: name
    real(DP), intent(in) :: expected, computed, abserr
    integer, intent(in) :: neval, ier
    real(DP) :: rel_err, abs_err

    abs_err = abs(computed - expected)
    if (expected /= 0.0_DP) then
      rel_err = abs_err / abs(expected)
    else
      rel_err = abs_err
    end if

    print "(A20,A)", name, ":"
    print "(A,F25.16)", "   Expected:        ", expected
    print "(A,F25.16)", "   Computed:        ", computed
    print "(A,Z16)",    "   Hex:             ", computed
    print "(A,ES12.5)", "   Abs error:       ", abs_err
    print "(A,ES12.5)", "   Rel error:       ", rel_err
    print "(A,ES12.5)", "   QUADPACK abserr: ", abserr
    print "(A,I8)",     "   Function evals:  ", neval
    print "(A,I8)",     "   Return code:     ", ier
  end subroutine

  subroutine print_result_sp(name, expected, computed, abserr, neval, ier)
    character(*), intent(in) :: name
    real(SP), intent(in) :: expected, computed, abserr
    integer, intent(in) :: neval, ier
    real(SP) :: rel_err, abs_err

    abs_err = abs(computed - expected)
    if (expected /= 0.0_SP) then
      rel_err = abs_err / abs(expected)
    else
      rel_err = abs_err
    end if

    print "(A20,A)", name, ":"
    print "(A,F18.10)", "   Expected:        ", expected
    print "(A,F18.10)", "   Computed:        ", computed
    print "(A,Z8)",     "   Hex:             ", computed
    print "(A,ES12.5)", "   Abs error:       ", abs_err
    print "(A,ES12.5)", "   Rel error:       ", rel_err
    print "(A,ES12.5)", "   QUADPACK abserr: ", abserr
    print "(A,I8)",     "   Function evals:  ", neval
    print "(A,I8)",     "   Return code:     ", ier
  end subroutine

end program quadpack_deviation

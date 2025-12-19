
module test_utils
  use service, only: SP, DP
  implicit none
  private

  public :: drift_report_sp, drift_report_dp
  public :: check_drift_sp, check_drift_dp
  public :: relative_error_sp, relative_error_dp
  public :: ulp_error_sp, ulp_error_dp, ulp_sp, ulp_dp, check_ulp_sp, check_ulp_dp

  type, public :: ulp_tolerances
    integer :: j0_dp=3, j0_sp=2, j1_dp=4, j1_sp=4, jn_dp=8, jn_sp=7
    integer :: y0_dp=3, y0_sp=3, y1_dp=5, y1_sp=4, yn_dp=6, yn_sp=7
    integer :: i_dp=4, i_sp=4, k_dp=4, k_sp=4
  end type

  type(ulp_tolerances), parameter, public :: GLIBC_ULP_TOL = ulp_tolerances()

  !> SLATEC-specific ULP tolerances
  !> SLATEC targets 10-14 significant digits, not machine epsilon
  !> These are more practical for testing SLATEC implementations
  type, public :: slatec_ulp_tolerances
    !> Double precision: allow 100 ULPs (~2e-14 relative)
    integer :: j_dp = 100
    integer :: i_dp = 100
    integer :: k_dp = 100000
    integer :: y_dp = 100
    !> Single precision: allow 200 ULPs (~2e-5 relative)
    integer :: j_sp = 200
    integer :: i_sp = 200
    integer :: k_sp = 1000
    integer :: y_sp = 200
  end type

  type(slatec_ulp_tolerances), parameter, public :: SLATEC_ULP_TOL = slatec_ulp_tolerances()

  type, public :: drift_stats_sp
    real(SP) :: max_abs_error=0.0_SP, max_rel_error=0.0_SP, mean_abs_error=0.0_SP, max_ulp_error=0.0_SP
    integer :: n_comparisons=0, n_failures=0
  end type

  type, public :: drift_stats_dp
    real(DP) :: max_abs_error=0.0_DP, max_rel_error=0.0_DP, mean_abs_error=0.0_DP, max_ulp_error=0.0_DP
    integer :: n_comparisons=0, n_failures=0
  end type

contains

  pure function ulp_sp(x) result(u)
    real(SP), intent(in) :: x
    real(SP) :: u
    if (abs(x) < tiny(1.0_SP)) then
      u = tiny(1.0_SP)
    else
      u = epsilon(1.0_SP) * abs(x)
    end if
  end function

  pure function ulp_dp(x) result(u)
    real(DP), intent(in) :: x
    real(DP) :: u
    if (abs(x) < tiny(1.0_DP)) then
      u = tiny(1.0_DP)
    else
      u = epsilon(1.0_DP) * abs(x)
    end if
  end function

  pure function ulp_error_sp(computed, reference) result(ulps)
    real(SP), intent(in) :: computed, reference
    real(SP) :: ulps, u
    u = ulp_sp(reference)
    if (u > 0.0_SP) then
      ulps = abs(computed - reference) / u
    else
      ulps = 0.0_SP
    end if
  end function

  pure function ulp_error_dp(computed, reference) result(ulps)
    real(DP), intent(in) :: computed, reference
    real(DP) :: ulps, u
    u = ulp_dp(reference)
    if (u > 0.0_DP) then
      ulps = abs(computed - reference) / u
    else
      ulps = 0.0_DP
    end if
  end function

  subroutine check_ulp_sp(reference, computed, max_ulps, stats, passed)
    real(SP), intent(in) :: reference(:), computed(:)
    integer, intent(in) :: max_ulps
    type(drift_stats_sp), intent(inout) :: stats
    logical, intent(out) :: passed
    integer :: i, n
    real(SP) :: ulps
    n = min(size(reference), size(computed))
    passed = .true.
    do i = 1, n
      ulps = ulp_error_sp(computed(i), reference(i))
      stats%n_comparisons = stats%n_comparisons + 1
      if (ulps > stats%max_ulp_error) stats%max_ulp_error = ulps
      if (ulps > real(max_ulps, SP)) then
        stats%n_failures = stats%n_failures + 1
        passed = .false.
      end if
    end do
  end subroutine

  subroutine check_ulp_dp(reference, computed, max_ulps, stats, passed)
    real(DP), intent(in) :: reference(:), computed(:)
    integer, intent(in) :: max_ulps
    type(drift_stats_dp), intent(inout) :: stats
    logical, intent(out) :: passed
    integer :: i, n
    real(DP) :: ulps
    n = min(size(reference), size(computed))
    passed = .true.
    do i = 1, n
      ulps = ulp_error_dp(computed(i), reference(i))
      stats%n_comparisons = stats%n_comparisons + 1
      if (ulps > stats%max_ulp_error) stats%max_ulp_error = ulps
      if (ulps > real(max_ulps, DP)) then
        stats%n_failures = stats%n_failures + 1
        passed = .false.
      end if
    end do
  end subroutine

  pure function relative_error_sp(computed, reference) result(rel_err)
    real(SP), intent(in) :: computed, reference
    real(SP) :: rel_err
    if (abs(reference) > 1.0E-38_SP) then
      rel_err = abs(computed - reference) / abs(reference)
    else
      rel_err = abs(computed - reference)
    end if
  end function

  pure function relative_error_dp(computed, reference) result(rel_err)
    real(DP), intent(in) :: computed, reference
    real(DP) :: rel_err
    if (abs(reference) > 1.0D-308) then
      rel_err = abs(computed - reference) / abs(reference)
    else
      rel_err = abs(computed - reference)
    end if
  end function

  subroutine check_drift_sp(original, modern, tol, stats, passed)
    real(SP), intent(in) :: original(:), modern(:)
    real(SP), intent(in) :: tol
    type(drift_stats_sp), intent(inout) :: stats
    logical, intent(out) :: passed
    integer :: i, n
    real(SP) :: abs_err, rel_err, ulps
    n = min(size(original), size(modern))
    passed = .true.
    do i = 1, n
      abs_err = abs(modern(i) - original(i))
      rel_err = relative_error_sp(modern(i), original(i))
      ulps = ulp_error_sp(modern(i), original(i))
      stats%n_comparisons = stats%n_comparisons + 1
      stats%mean_abs_error = stats%mean_abs_error + abs_err
      if (abs_err > stats%max_abs_error) stats%max_abs_error = abs_err
      if (rel_err > stats%max_rel_error) stats%max_rel_error = rel_err
      if (ulps > stats%max_ulp_error) stats%max_ulp_error = ulps
      if (rel_err > tol .and. abs_err > tol) then
        stats%n_failures = stats%n_failures + 1
        passed = .false.
      end if
    end do
    if (stats%n_comparisons > 0) stats%mean_abs_error = stats%mean_abs_error / stats%n_comparisons
  end subroutine

  subroutine check_drift_dp(original, modern, tol, stats, passed)
    real(DP), intent(in) :: original(:), modern(:)
    real(DP), intent(in) :: tol
    type(drift_stats_dp), intent(inout) :: stats
    logical, intent(out) :: passed
    integer :: i, n
    real(DP) :: abs_err, rel_err, ulps
    n = min(size(original), size(modern))
    passed = .true.
    do i = 1, n
      abs_err = abs(modern(i) - original(i))
      rel_err = relative_error_dp(modern(i), original(i))
      ulps = ulp_error_dp(modern(i), original(i))
      stats%n_comparisons = stats%n_comparisons + 1
      stats%mean_abs_error = stats%mean_abs_error + abs_err
      if (abs_err > stats%max_abs_error) stats%max_abs_error = abs_err
      if (rel_err > stats%max_rel_error) stats%max_rel_error = rel_err
      if (ulps > stats%max_ulp_error) stats%max_ulp_error = ulps
      if (rel_err > tol .and. abs_err > tol) then
        stats%n_failures = stats%n_failures + 1
        passed = .false.
      end if
    end do
    if (stats%n_comparisons > 0) stats%mean_abs_error = stats%mean_abs_error / stats%n_comparisons
  end subroutine

  subroutine drift_report_sp(label, stats, unit)
    character(len=*), intent(in) :: label
    type(drift_stats_sp), intent(in) :: stats
    integer, intent(in), optional :: unit
    integer :: u
    u = 6
    if (present(unit)) u = unit
    write(u, '(a)') "=== Drift Report: " // trim(label) // " ==="
    write(u, '(a, i0)') "  Comparisons: ", stats%n_comparisons
    write(u, '(a, es12.5)') "  Max absolute error: ", stats%max_abs_error
    write(u, '(a, es12.5)') "  Max relative error: ", stats%max_rel_error
    write(u, '(a, f8.2, a)') "  Max ULP error: ", stats%max_ulp_error, " ULPs"
    write(u, '(a, es12.5)') "  Mean absolute error: ", stats%mean_abs_error
    write(u, '(a, i0)') "  Failures: ", stats%n_failures
    write(u, '(a)') ""
  end subroutine

  subroutine drift_report_dp(label, stats, unit)
    character(len=*), intent(in) :: label
    type(drift_stats_dp), intent(in) :: stats
    integer, intent(in), optional :: unit
    integer :: u
    u = 6
    if (present(unit)) u = unit
    write(u, '(a)') "=== Drift Report: " // trim(label) // " ==="
    write(u, '(a, i0)') "  Comparisons: ", stats%n_comparisons
    write(u, '(a, es12.5)') "  Max absolute error: ", stats%max_abs_error
    write(u, '(a, es12.5)') "  Max relative error: ", stats%max_rel_error
    write(u, '(a, f8.2, a)') "  Max ULP error: ", stats%max_ulp_error, " ULPs"
    write(u, '(a, es12.5)') "  Mean absolute error: ", stats%mean_abs_error
    write(u, '(a, i0)') "  Failures: ", stats%n_failures
    write(u, '(a)') ""
  end subroutine

end module test_utils

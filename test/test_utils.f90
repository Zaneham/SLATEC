module test_utils
  !> Utility routines for numeric drift analysis and test comparison
  use service, only: SP, DP
  implicit none
  private

  public :: drift_report_sp, drift_report_dp
  public :: check_drift_sp, check_drift_dp
  public :: relative_error_sp, relative_error_dp

  type, public :: drift_stats_sp
    real(SP) :: max_abs_error = 0.0_SP
    real(SP) :: max_rel_error = 0.0_SP
    real(SP) :: mean_abs_error = 0.0_SP
    integer :: n_comparisons = 0
    integer :: n_failures = 0
  end type

  type, public :: drift_stats_dp
    real(DP) :: max_abs_error = 0.0_DP
    real(DP) :: max_rel_error = 0.0_DP
    real(DP) :: mean_abs_error = 0.0_DP
    integer :: n_comparisons = 0
    integer :: n_failures = 0
  end type

contains

  pure function relative_error_sp(computed, reference) result(rel_err)
    real(SP), intent(in) :: computed, reference
    real(SP) :: rel_err
    real(SP), parameter :: tiny_val = 1.0E-38_SP
    if (abs(reference) > tiny_val) then
      rel_err = abs(computed - reference) / abs(reference)
    else
      rel_err = abs(computed - reference)
    end if
  end function

  pure function relative_error_dp(computed, reference) result(rel_err)
    real(DP), intent(in) :: computed, reference
    real(DP) :: rel_err
    real(DP), parameter :: tiny_val = 1.0D-308
    if (abs(reference) > tiny_val) then
      rel_err = abs(computed - reference) / abs(reference)
    else
      rel_err = abs(computed - reference)
    end if
  end function

  subroutine check_drift_sp(original, modern, tol, stats, passed)
    !> Compare original and modern values, update drift statistics
    real(SP), intent(in) :: original(:), modern(:)
    real(SP), intent(in) :: tol
    type(drift_stats_sp), intent(inout) :: stats
    logical, intent(out) :: passed
    integer :: i, n
    real(SP) :: abs_err, rel_err

    n = min(size(original), size(modern))
    passed = .true.

    do i = 1, n
      abs_err = abs(modern(i) - original(i))
      rel_err = relative_error_sp(modern(i), original(i))

      stats%n_comparisons = stats%n_comparisons + 1
      stats%mean_abs_error = stats%mean_abs_error + abs_err

      if (abs_err > stats%max_abs_error) stats%max_abs_error = abs_err
      if (rel_err > stats%max_rel_error) stats%max_rel_error = rel_err

      if (rel_err > tol .and. abs_err > tol) then
        stats%n_failures = stats%n_failures + 1
        passed = .false.
      end if
    end do

    if (stats%n_comparisons > 0) then
      stats%mean_abs_error = stats%mean_abs_error / stats%n_comparisons
    end if
  end subroutine

  subroutine check_drift_dp(original, modern, tol, stats, passed)
    !> Compare original and modern values, update drift statistics
    real(DP), intent(in) :: original(:), modern(:)
    real(DP), intent(in) :: tol
    type(drift_stats_dp), intent(inout) :: stats
    logical, intent(out) :: passed
    integer :: i, n
    real(DP) :: abs_err, rel_err

    n = min(size(original), size(modern))
    passed = .true.

    do i = 1, n
      abs_err = abs(modern(i) - original(i))
      rel_err = relative_error_dp(modern(i), original(i))

      stats%n_comparisons = stats%n_comparisons + 1
      stats%mean_abs_error = stats%mean_abs_error + abs_err

      if (abs_err > stats%max_abs_error) stats%max_abs_error = abs_err
      if (rel_err > stats%max_rel_error) stats%max_rel_error = rel_err

      if (rel_err > tol .and. abs_err > tol) then
        stats%n_failures = stats%n_failures + 1
        passed = .false.
      end if
    end do

    if (stats%n_comparisons > 0) then
      stats%mean_abs_error = stats%mean_abs_error / stats%n_comparisons
    end if
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
    write(u, '(a, es12.5)') "  Mean absolute error: ", stats%mean_abs_error
    write(u, '(a, i0)') "  Failures: ", stats%n_failures
    write(u, '(a)') ""
  end subroutine

end module test_utils

program deviation_report
  !> ZH: Deviation documentation program - EXPANDED
  !> Outputs exact computed values at maximum precision for cross-platform comparison
  !> This is NOT a pass/fail test - it documents what each configuration produces

  use service, only: SP, DP
  use special_functions, only: BESJ, DBESJ, BESI, DBESI, BESK, DBESK, DBESY, BESY
  implicit none

  print *, "============================================================"
  print *, "SLATEC-MODERN DEVIATION REPORT (EXPANDED)"
  print *, "============================================================"
  print *, ""
#if defined(__GFORTRAN__)
  print "(A,I0,A,I0,A,I0)", " Compiler: gfortran ", __GNUC__, ".", __GNUC_MINOR__, ".", __GNUC_PATCHLEVEL__
#elif defined(__INTEL_COMPILER)
  print "(A,I0)", " Compiler: Intel ", __INTEL_COMPILER
#elif defined(__flang__)
  print *, " Compiler: Flang"
#else
  print *, " Compiler: unknown"
#endif
  print "(A,I0)", " SP precision (digits): ", precision(1.0_SP)
  print "(A,I0)", " DP precision (digits): ", precision(1.0_DP)
  print "(A,ES10.3)", " SP epsilon: ", epsilon(1.0_SP)
  print "(A,ES10.3)", " DP epsilon: ", epsilon(1.0_DP)
  print *, ""

  ! Bessel functions - expanded
  call report_bessel_j_expanded()
  call report_bessel_i_expanded()
  call report_bessel_k_expanded()
  call report_bessel_y_expanded()

  print *, ""
  print *, "============================================================"
  print *, "END REPORT"
  print *, "============================================================"

contains

  subroutine print_dp_value(name, ref, val)
    character(*), intent(in) :: name, ref
    real(DP), intent(in) :: val
    print *, name
    print "(A,A)", "  Reference:        ", ref
    print "(A,F25.16)", "   DP computed:     ", val
    print "(A,Z16)",    "   DP hex:          ", val
    print *, ""
  end subroutine

  subroutine print_sp_dp_value(name, ref, sp_val, dp_val)
    character(*), intent(in) :: name, ref
    real(SP), intent(in) :: sp_val
    real(DP), intent(in) :: dp_val
    print *, name
    print "(A,A)", "  Reference:        ", ref
    print "(A,F25.16)", "   SP computed:     ", sp_val
    print "(A,F25.16)", "   DP computed:     ", dp_val
    print "(A,Z8)",     "   SP hex:          ", sp_val
    print "(A,Z16)",    "   DP hex:          ", dp_val
    print *, ""
  end subroutine

  subroutine report_bessel_j_expanded()
    real(SP) :: ys(10)
    real(DP) :: yd(10)
    integer :: nz

    print *, "================================================================"
    print *, "--- Bessel J (First Kind) - EXPANDED ---"
    print *, "================================================================"
    print *, ""

    print *, "=== Standard test points (A&S Table 9.1) ==="
    print *, ""

    call BESJ(1.0_SP, 0.0_SP, 1, ys, nz)
    call DBESJ(1.0_DP, 0.0_DP, 1, yd, nz)
    call print_sp_dp_value("J_0(1.0):", "0.7651976865579666", ys(1), yd(1))

    call BESJ(2.0_SP, 0.0_SP, 1, ys, nz)
    call DBESJ(2.0_DP, 0.0_DP, 1, yd, nz)
    call print_sp_dp_value("J_0(2.0):", "0.2238907791412357", ys(1), yd(1))

    call BESJ(5.0_SP, 0.0_SP, 1, ys, nz)
    call DBESJ(5.0_DP, 0.0_DP, 1, yd, nz)
    call print_sp_dp_value("J_0(5.0):", "-0.1775967713143383", ys(1), yd(1))

    call BESJ(1.0_SP, 0.0_SP, 2, ys, nz)
    call DBESJ(1.0_DP, 0.0_DP, 2, yd, nz)
    call print_sp_dp_value("J_1(1.0):", "0.4400505857449335", ys(2), yd(2))

    call BESJ(2.0_SP, 0.0_SP, 2, ys, nz)
    call DBESJ(2.0_DP, 0.0_DP, 2, yd, nz)
    call print_sp_dp_value("J_1(2.0):", "0.5767248077568734", ys(2), yd(2))

    print *, "=== Near zeros of J_0 (sensitive points) ==="
    print *, ""

    ! First zero of J_0 is at x = 2.4048255576957728...
    call DBESJ(2.4048255576957728_DP, 0.0_DP, 1, yd, nz)
    call print_dp_value("J_0(2.404825557695773) [1st zero]:", "0.0", yd(1))

    ! Second zero at x = 5.5200781102863106...
    call DBESJ(5.5200781102863106_DP, 0.0_DP, 1, yd, nz)
    call print_dp_value("J_0(5.520078110286311) [2nd zero]:", "0.0", yd(1))

    ! Third zero at x = 8.6537279129110122...
    call DBESJ(8.6537279129110122_DP, 0.0_DP, 1, yd, nz)
    call print_dp_value("J_0(8.653727912911012) [3rd zero]:", "0.0", yd(1))

    print *, "=== Near-zero argument (series expansion region) ==="
    print *, ""

    call DBESJ(1.0D-10, 0.0_DP, 1, yd, nz)
    call print_dp_value("J_0(1E-10):", "1.0 (exact)", yd(1))

    call DBESJ(1.0D-15, 0.0_DP, 1, yd, nz)
    call print_dp_value("J_0(1E-15):", "1.0 (exact)", yd(1))

    call DBESJ(1.0D-10, 0.0_DP, 2, yd, nz)
    call print_dp_value("J_1(1E-10):", "5E-11 (x/2)", yd(2))

    print *, "=== Large argument (asymptotic region) ==="
    print *, ""

    call BESJ(50.0_SP, 0.0_SP, 1, ys, nz)
    call DBESJ(50.0_DP, 0.0_DP, 1, yd, nz)
    call print_sp_dp_value("J_0(50.0):", "0.0558123276692332 (NIST)", ys(1), yd(1))

    call BESJ(100.0_SP, 0.0_SP, 1, ys, nz)
    call DBESJ(100.0_DP, 0.0_DP, 1, yd, nz)
    call print_sp_dp_value("J_0(100.0):", "0.0199858503042232 (NIST)", ys(1), yd(1))

    call DBESJ(500.0_DP, 0.0_DP, 1, yd, nz)
    call print_dp_value("J_0(500.0):", "0.0179753822938868 (NIST)", yd(1))

    call DBESJ(1000.0_DP, 0.0_DP, 1, yd, nz)
    call print_dp_value("J_0(1000.0):", "-0.0246711376936846 (NIST)", yd(1))

    print *, "=== Higher orders ==="
    print *, ""

    call DBESJ(5.0_DP, 0.0_DP, 5, yd, nz)
    call print_dp_value("J_4(5.0):", "0.3912057568779922 (A&S)", yd(5))

    call DBESJ(10.0_DP, 0.0_DP, 6, yd, nz)
    call print_dp_value("J_5(10.0):", "-0.2340615281867937 (NIST)", yd(6))

    call DBESJ(20.0_DP, 0.0_DP, 10, yd, nz)
    call print_dp_value("J_9(20.0):", "0.2453202954891665 (NIST)", yd(10))

  end subroutine

  subroutine report_bessel_i_expanded()
    real(SP) :: ys(5)
    real(DP) :: yd(5)
    integer :: nz

    print *, "================================================================"
    print *, "--- Bessel I (Modified First Kind) - EXPANDED ---"
    print *, "================================================================"
    print *, ""

    print *, "=== Standard test points (A&S Table 9.8) ==="
    print *, ""

    call BESI(1.0_SP, 0.0_SP, 1, 1, ys, nz)
    call DBESI(1.0_DP, 0.0_DP, 1, 1, yd, nz)
    call print_sp_dp_value("I_0(1.0):", "1.2660658777520084", ys(1), yd(1))

    call BESI(2.0_SP, 0.0_SP, 1, 1, ys, nz)
    call DBESI(2.0_DP, 0.0_DP, 1, 1, yd, nz)
    call print_sp_dp_value("I_0(2.0):", "2.2795853023360673", ys(1), yd(1))

    call BESI(1.0_SP, 0.0_SP, 1, 2, ys, nz)
    call DBESI(1.0_DP, 0.0_DP, 1, 2, yd, nz)
    call print_sp_dp_value("I_1(1.0):", "0.5651591039924851", ys(2), yd(2))

    call BESI(2.0_SP, 0.0_SP, 1, 2, ys, nz)
    call DBESI(2.0_DP, 0.0_DP, 1, 2, yd, nz)
    call print_sp_dp_value("I_1(2.0):", "1.5906368546373291", ys(2), yd(2))

    print *, "=== Near-zero argument ==="
    print *, ""

    call DBESI(1.0D-10, 0.0_DP, 1, 1, yd, nz)
    call print_dp_value("I_0(1E-10):", "1.0 (exact)", yd(1))

    call DBESI(1.0D-10, 0.0_DP, 1, 2, yd, nz)
    call print_dp_value("I_1(1E-10):", "5E-11 (x/2)", yd(2))

    print *, "=== Large argument (exponential growth) ==="
    print *, ""

    call DBESI(10.0_DP, 0.0_DP, 1, 1, yd, nz)
    call print_dp_value("I_0(10.0):", "2815.7166284662544 (NIST)", yd(1))

    call DBESI(50.0_DP, 0.0_DP, 1, 1, yd, nz)
    call print_dp_value("I_0(50.0):", "2.9325537838493547E+20 (NIST)", yd(1))

    print *, "=== Higher orders ==="
    print *, ""

    call DBESI(5.0_DP, 0.0_DP, 1, 3, yd, nz)
    call print_dp_value("I_2(5.0):", "17.505614966624236 (NIST)", yd(3))

  end subroutine

  subroutine report_bessel_k_expanded()
    real(SP) :: ys(5)
    real(DP) :: yd(5)
    integer :: nz

    print *, "================================================================"
    print *, "--- Bessel K (Modified Second Kind) - EXPANDED ---"
    print *, "*** THIS FUNCTION SHOWS OPTIMISATION-DEPENDENT DEVIATIONS ***"
    print *, "================================================================"
    print *, ""

    print *, "=== Standard test points (A&S Table 9.8) ==="
    print *, ""

    call BESK(1.0_SP, 0.0_SP, 1, 1, ys, nz)
    call DBESK(1.0_DP, 0.0_DP, 1, 1, yd, nz)
    call print_sp_dp_value("K_0(1.0): *** KNOWN DEVIATION ***", "0.4210244382407084", ys(1), yd(1))

    call BESK(2.0_SP, 0.0_SP, 1, 1, ys, nz)
    call DBESK(2.0_DP, 0.0_DP, 1, 1, yd, nz)
    call print_sp_dp_value("K_0(2.0): *** KNOWN DEVIATION ***", "0.1138938727495334", ys(1), yd(1))

    call BESK(1.0_SP, 0.0_SP, 1, 2, ys, nz)
    call DBESK(1.0_DP, 0.0_DP, 1, 2, yd, nz)
    call print_sp_dp_value("K_1(1.0):", "0.6019072301972346", ys(2), yd(2))

    call BESK(2.0_SP, 0.0_SP, 1, 2, ys, nz)
    call DBESK(2.0_DP, 0.0_DP, 1, 2, yd, nz)
    call print_sp_dp_value("K_1(2.0):", "0.1398658818165224", ys(2), yd(2))

    print *, "=== Near-zero argument (logarithmic singularity) ==="
    print *, ""

    call DBESK(0.1_DP, 0.0_DP, 1, 1, yd, nz)
    call print_dp_value("K_0(0.1):", "2.4270690247020166 (NIST)", yd(1))

    call DBESK(0.01_DP, 0.0_DP, 1, 1, yd, nz)
    call print_dp_value("K_0(0.01):", "4.7212359359101134 (NIST)", yd(1))

    call DBESK(1.0D-5, 0.0_DP, 1, 1, yd, nz)
    call print_dp_value("K_0(1E-5):", "~11.629 (-ln(x/2)-gamma)", yd(1))

    print *, "=== Large argument (exponential decay) ==="
    print *, ""

    call DBESK(10.0_DP, 0.0_DP, 1, 1, yd, nz)
    call print_dp_value("K_0(10.0):", "1.7780062316167652E-05 (NIST)", yd(1))

    call DBESK(50.0_DP, 0.0_DP, 1, 1, yd, nz)
    call print_dp_value("K_0(50.0):", "3.4101677497894956E-23 (NIST)", yd(1))

    print *, "=== Additional K_0 points for deviation tracking ==="
    print *, ""

    call BESK(0.5_SP, 0.0_SP, 1, 1, ys, nz)
    call DBESK(0.5_DP, 0.0_DP, 1, 1, yd, nz)
    call print_sp_dp_value("K_0(0.5):", "0.9244190712276659 (NIST)", ys(1), yd(1))

    call BESK(1.5_SP, 0.0_SP, 1, 1, ys, nz)
    call DBESK(1.5_DP, 0.0_DP, 1, 1, yd, nz)
    call print_sp_dp_value("K_0(1.5):", "0.2138055626545264 (NIST)", ys(1), yd(1))

    call BESK(3.0_SP, 0.0_SP, 1, 1, ys, nz)
    call DBESK(3.0_DP, 0.0_DP, 1, 1, yd, nz)
    call print_sp_dp_value("K_0(3.0):", "0.0347395043862552 (NIST)", ys(1), yd(1))

    call BESK(5.0_SP, 0.0_SP, 1, 1, ys, nz)
    call DBESK(5.0_DP, 0.0_DP, 1, 1, yd, nz)
    call print_sp_dp_value("K_0(5.0):", "0.0036910982867735 (NIST)", ys(1), yd(1))

  end subroutine

  subroutine report_bessel_y_expanded()
    real(SP) :: ys(5)
    real(DP) :: yd(5)

    print *, "================================================================"
    print *, "--- Bessel Y (Second Kind) - EXPANDED ---"
    print *, "================================================================"
    print *, ""

    print *, "=== Standard test points (A&S Table 9.1) ==="
    print *, ""

    call BESY(1.0_SP, 0.0_SP, 1, ys)
    call DBESY(1.0_DP, 0.0_DP, 1, yd)
    call print_sp_dp_value("Y_0(1.0):", "0.0882569642156769", ys(1), yd(1))

    call BESY(2.0_SP, 0.0_SP, 1, ys)
    call DBESY(2.0_DP, 0.0_DP, 1, yd)
    call print_sp_dp_value("Y_0(2.0):", "0.5103756726497451", ys(1), yd(1))

    call BESY(5.0_SP, 0.0_SP, 1, ys)
    call DBESY(5.0_DP, 0.0_DP, 1, yd)
    call print_sp_dp_value("Y_0(5.0):", "-0.3085176252490338", ys(1), yd(1))

    call BESY(1.0_SP, 0.0_SP, 2, ys)
    call DBESY(1.0_DP, 0.0_DP, 2, yd)
    call print_sp_dp_value("Y_1(1.0):", "-0.7812128213002887", ys(2), yd(2))

    print *, "=== Near zeros of Y_0 ==="
    print *, ""

    ! First zero of Y_0 at x = 0.8935769662791675...
    call DBESY(0.8935769662791675_DP, 0.0_DP, 1, yd)
    call print_dp_value("Y_0(0.893576966...) [1st zero]:", "0.0", yd(1))

    ! Second zero at x = 3.9576784193148578...
    call DBESY(3.9576784193148578_DP, 0.0_DP, 1, yd)
    call print_dp_value("Y_0(3.957678419...) [2nd zero]:", "0.0", yd(1))

    print *, "=== Large argument ==="
    print *, ""

    call BESY(50.0_SP, 0.0_SP, 1, ys)
    call DBESY(50.0_DP, 0.0_DP, 1, yd)
    call print_sp_dp_value("Y_0(50.0):", "-0.0560404718523358 (NIST)", ys(1), yd(1))

    call DBESY(100.0_DP, 0.0_DP, 1, yd)
    call print_dp_value("Y_0(100.0):", "-0.0772433752531550 (NIST)", yd(1))

  end subroutine

end program deviation_report

program standalone_portability
  !> ZH: Standalone portability test - simulates different hardware by
  !> comparing results at extreme ranges and between SP/DP
  
  use service, only: SP, DP
  use special_functions, only: BESJ, DBESJ, BESI, DBESI, BESK, DBESK
  implicit none
  
  print *, "=== Bessel Portability Tests (Standalone) ==="
  print *, "Simulating hardware variations via extreme cases"
  print *
  
  call test_sp_vs_dp()
  call test_large_argument()
  call test_small_argument()
  
  print *
  print *, "Portability tests complete."
  
contains

  subroutine test_sp_vs_dp()
    !> Compare single and double precision - should agree to SP precision
    real(SP) :: ys(1)
    real(DP) :: yd(1)
    integer :: nz
    real(DP) :: diff, sp_eps
    
    sp_eps = epsilon(1.0_SP)  ! ~1.2E-7
    
    print *, "Single vs Double precision agreement (SP eps ~1.2E-7):"
    
    ! J_0(1)
    call BESJ(1.0_SP, 0.0_SP, 1, ys, nz)
    call DBESJ(1.0_DP, 0.0_DP, 1, yd, nz)
    diff = abs(real(ys(1), DP) - yd(1))
    print "(A,F12.9,A,F18.15,A,ES10.3,A,F6.1,A)", &
      "  J_0(1): SP=", real(ys(1),DP), "  DP=", yd(1), &
      "  diff=", diff, " (", diff/sp_eps, " SP ulps)"
    
    ! I_0(2)
    call BESI(2.0_SP, 0.0_SP, 1, 1, ys, nz)
    call DBESI(2.0_DP, 0.0_DP, 1, 1, yd, nz)
    diff = abs(real(ys(1), DP) - yd(1))
    print "(A,F12.9,A,F18.15,A,ES10.3,A,F6.1,A)", &
      "  I_0(2): SP=", real(ys(1),DP), "  DP=", yd(1), &
      "  diff=", diff, " (", diff/sp_eps, " SP ulps)"
    
    ! K_0(1)
    call BESK(1.0_SP, 0.0_SP, 1, 1, ys, nz)
    call DBESK(1.0_DP, 0.0_DP, 1, 1, yd, nz)
    diff = abs(real(ys(1), DP) - yd(1))
    print "(A,F12.9,A,F18.15,A,ES10.3,A,F6.1,A)", &
      "  K_0(1): SP=", real(ys(1),DP), "  DP=", yd(1), &
      "  diff=", diff, " (", diff/sp_eps, " SP ulps)"
  end subroutine

  subroutine test_large_argument()
    !> Large arguments stress FPU - tests range reduction algorithms
    real(DP) :: y(1)
    integer :: nz
    
    print *
    print *, "Large argument (FPU range reduction stress test):"
    
    ! J_0(50) - oscillatory, tests trig range reduction
    call DBESJ(50.0_DP, 0.0_DP, 1, y, nz)
    print "(A,F20.16)", "  J_0(50)  = ", y(1)
    
    ! J_0(100) - more stress on range reduction
    call DBESJ(100.0_DP, 0.0_DP, 1, y, nz)
    print "(A,F20.16)", "  J_0(100) = ", y(1)
    
    ! I_0(700) - near overflow boundary
    call DBESI(700.0_DP, 0.0_DP, 1, 1, y, nz)
    print "(A,ES23.16,A,I0)", "  I_0(700) = ", y(1), "  (nz=", nz, ")"
    
    ! K_0(700) - underflow region
    call DBESK(700.0_DP, 0.0_DP, 1, 1, y, nz)
    print "(A,ES23.16,A,I0)", "  K_0(700) = ", y(1), "  (nz=", nz, ")"
  end subroutine

  subroutine test_small_argument()
    !> Small arguments test series expansion - different FPU rounding
    real(DP) :: y(1)
    integer :: nz
    
    print *
    print *, "Small argument (series expansion, rounding sensitivity):"
    
    ! J_0 near zero -> 1
    call DBESJ(1.0D-10, 0.0_DP, 1, y, nz)
    print "(A,F23.20,A,ES10.3)", "  J_0(1E-10) = ", y(1), "  err from 1: ", abs(y(1)-1.0_DP)
    
    ! I_0 near zero -> 1
    call DBESI(1.0D-10, 0.0_DP, 1, 1, y, nz)
    print "(A,F23.20,A,ES10.3)", "  I_0(1E-10) = ", y(1), "  err from 1: ", abs(y(1)-1.0_DP)
    
    ! K_0 near zero -> -ln(x/2) - gamma ~ 11.5 + 0.577
    call DBESK(1.0D-5, 0.0_DP, 1, 1, y, nz)
    print "(A,F20.16)", "  K_0(1E-5) = ", y(1)
  end subroutine

end program standalone_portability

program detailed_compare
  use service, only: SP, DP
  use special_functions, only: DBESJ, DBESI, DBESK
  implicit none
  
  real(DP) :: y(1)
  integer :: nz
  
  print *, "Double precision values at high precision:"
  
  call DBESJ(1.0_DP, 0.0_DP, 1, y, nz)
  print "(A,F25.22)", "  J_0(1)   = ", y(1)
  
  call DBESJ(50.0_DP, 0.0_DP, 1, y, nz)
  print "(A,F25.22)", "  J_0(50)  = ", y(1)
  
  call DBESJ(100.0_DP, 0.0_DP, 1, y, nz)
  print "(A,F25.22)", "  J_0(100) = ", y(1)
  
  call DBESK(1.0_DP, 0.0_DP, 1, 1, y, nz)
  print "(A,F25.22)", "  K_0(1)   = ", y(1)
  
  call DBESI(2.0_DP, 0.0_DP, 1, 1, y, nz)
  print "(A,F25.22)", "  I_0(2)   = ", y(1)

end program detailed_compare

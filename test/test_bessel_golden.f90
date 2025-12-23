program test_bessel_golden
  !> ZH: Golden regression tests for Bessel functions
  !>
  !> Reference values from:
  !>   - Abramowitz & Stegun, Handbook of Mathematical Functions (1964)
  !>     Table 9.1 (J Bessel), Table 9.8 (I,K Bessel)
  !>   - NIST Digital Library of Mathematical Functions (dlmf.nist.gov)
  !>
  !> These are standalone tests that do not require external BLAS/LAPACK.

  use service, only: DP
  use special_functions, only: DBESJ, DBESI, DBESK
  implicit none

  integer :: nfail, npass
  
  nfail = 0
  npass = 0

  print *, "=== Bessel Function Golden Regression Tests ==="
  print *, "Reference: Abramowitz & Stegun Tables 9.1, 9.8"
  print *

  ! J Bessel tests (A&S Table 9.1)
  call test_j0_values()
  call test_j1_values()
  
  ! I Bessel tests (A&S Table 9.8)
  call test_i0_values()
  call test_i1_values()
  
  ! K Bessel tests (A&S Table 9.8)
  call test_k0_values()
  call test_k1_values()

  print *
  print *, "=== Summary ==="
  print "(A,I0,A,I0)", " Passed: ", npass, "  Failed: ", nfail
  
  if (nfail > 0) then
    error stop "Some tests failed"
  else
    print *, "All tests passed!"
  end if

contains

  subroutine test_j0_values()
    !> J_0(x) values from A&S Table 9.1
    real(DP) :: y(1)
    integer :: nz
    logical :: ok
    
    print *, "Testing J_0(x)..."
    
    ! J_0(0) = 1.0 exactly
    call DBESJ(1.0D-15, 0.0_DP, 1, y, nz)
    call check_value("J_0(0)", y(1), 1.0_DP, 1.0D-10, ok)
    
    ! J_0(1) = 0.7651976865579666 (A&S)
    call DBESJ(1.0_DP, 0.0_DP, 1, y, nz)
    call check_value("J_0(1)", y(1), 0.7651976865579666_DP, 1.0D-12, ok)
    
    ! J_0(2) = 0.2238907791412357 (A&S)
    call DBESJ(2.0_DP, 0.0_DP, 1, y, nz)
    call check_value("J_0(2)", y(1), 0.2238907791412357_DP, 1.0D-12, ok)
    
    ! J_0(5) = -0.1775967713143383 (A&S)
    call DBESJ(5.0_DP, 0.0_DP, 1, y, nz)
    call check_value("J_0(5)", y(1), -0.1775967713143383_DP, 1.0D-12, ok)
  end subroutine

  subroutine test_j1_values()
    !> J_1(x) values from A&S Table 9.1
    real(DP) :: y(2)
    integer :: nz
    logical :: ok
    
    print *, "Testing J_1(x)..."
    
    ! J_1(1) = 0.4400505857449335 (A&S)
    call DBESJ(1.0_DP, 0.0_DP, 2, y, nz)
    call check_value("J_1(1)", y(2), 0.4400505857449335_DP, 1.0D-12, ok)
    
    ! J_1(2) = 0.5767248077568734 (A&S)
    call DBESJ(2.0_DP, 0.0_DP, 2, y, nz)
    call check_value("J_1(2)", y(2), 0.5767248077568734_DP, 1.0D-12, ok)
  end subroutine

  subroutine test_i0_values()
    !> I_0(x) values from A&S Table 9.8
    real(DP) :: y(1)
    integer :: nz
    logical :: ok
    
    print *, "Testing I_0(x)..."
    
    ! I_0(1) = 1.2660658777520084 (A&S)
    call DBESI(1.0_DP, 0.0_DP, 1, 1, y, nz)
    call check_value("I_0(1)", y(1), 1.2660658777520084_DP, 1.0D-12, ok)
    
    ! I_0(2) = 2.2795853023360673 (A&S)
    call DBESI(2.0_DP, 0.0_DP, 1, 1, y, nz)
    call check_value("I_0(2)", y(1), 2.2795853023360673_DP, 1.0D-12, ok)
  end subroutine

  subroutine test_i1_values()
    !> I_1(x) values from A&S Table 9.8
    real(DP) :: y(2)
    integer :: nz
    logical :: ok
    
    print *, "Testing I_1(x)..."
    
    ! I_1(1) = 0.5651591039924851 (A&S)
    call DBESI(1.0_DP, 0.0_DP, 1, 2, y, nz)
    call check_value("I_1(1)", y(2), 0.5651591039924851_DP, 1.0D-12, ok)
    
    ! I_1(2) = 1.5906368546373291 (A&S)
    call DBESI(2.0_DP, 0.0_DP, 1, 2, y, nz)
    call check_value("I_1(2)", y(2), 1.5906368546373291_DP, 1.0D-12, ok)
  end subroutine

  subroutine test_k0_values()
    !> K_0(x) values from A&S Table 9.8
    real(DP) :: y(1)
    integer :: nz
    logical :: ok
    
    print *, "Testing K_0(x)..."
    
    ! K_0(1) = 0.4210244382407084 (A&S)
    call DBESK(1.0_DP, 0.0_DP, 1, 1, y, nz)
    call check_value("K_0(1)", y(1), 0.4210244382407084_DP, 1.0D-12, ok)
    
    ! K_0(2) = 0.1138938727495334 (A&S)
    call DBESK(2.0_DP, 0.0_DP, 1, 1, y, nz)
    call check_value("K_0(2)", y(1), 0.1138938727495334_DP, 1.0D-12, ok)
  end subroutine

  subroutine test_k1_values()
    !> K_1(x) values from A&S Table 9.8
    real(DP) :: y(2)
    integer :: nz
    logical :: ok
    
    print *, "Testing K_1(x)..."
    
    ! K_1(1) = 0.6019072301972346 (A&S)
    call DBESK(1.0_DP, 0.0_DP, 1, 2, y, nz)
    call check_value("K_1(1)", y(2), 0.6019072301972346_DP, 1.0D-12, ok)
    
    ! K_1(2) = 0.1398658818165224 (A&S)
    call DBESK(2.0_DP, 0.0_DP, 1, 2, y, nz)
    call check_value("K_1(2)", y(2), 0.1398658818165224_DP, 1.0D-12, ok)
  end subroutine

  subroutine check_value(name, computed, expected, tol, ok)
    character(*), intent(in) :: name
    real(DP), intent(in) :: computed, expected, tol
    logical, intent(out) :: ok
    real(DP) :: rel_err, abs_err
    character(4) :: status

    abs_err = abs(computed - expected)
    if (expected /= 0.0_DP) then
      rel_err = abs_err / abs(expected)
    else
      rel_err = abs_err
    end if

    ok = (rel_err < tol)

    if (ok) then
      status = "PASS"
      npass = npass + 1
    else
      status = "FAIL"
      nfail = nfail + 1
    end if

    ! Always print: name, A&S reference, computed, absolute error, relative error
    print "(A4,1X,A8,A,F19.16,A,F19.16,A,ES10.3,A,ES10.3)", &
      status, name, &
      "  ref=", expected, "  got=", computed, &
      "  abs=", abs_err, "  rel=", rel_err
  end subroutine

end program test_bessel_golden

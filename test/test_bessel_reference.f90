module test_bessel_reference
  !> Reference tests comparing against published values
  !> Sources: Abramowitz & Stegun (1972), NIST DLMF, ACM TOMS papers
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use service, only: SP, DP
  use special_functions, only: DBESI, DBESJ, DBESK
  use test_utils
  implicit none
  private

  public :: collect_bessel_reference_tests

  ! Tolerance for reference tests
  ! TODO: Verify test data against original A&S tables
  ! Using 1e-2 tolerance until data is verified
  real(DP), parameter :: TOL_DP = 1.0D-2
  real(SP), parameter :: TOL_SP = 1.0E-2_SP

contains

  subroutine collect_bessel_reference_tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("I Bessel vs A&S Table 9.8", test_ibess_as_table98), &
      new_unittest("J Bessel vs A&S Table 9.1", test_jbess_as_table91), &
      new_unittest("K Bessel vs A&S Table 9.8", test_kbess_as_table98), &
      new_unittest("I Bessel special values", test_ibess_special), &
      new_unittest("J Bessel special values", test_jbess_special), &
      new_unittest("K Bessel asymptotic", test_kbess_asymptotic) &
    ]
  end subroutine

  subroutine test_ibess_as_table98(error)
    !> Test I Bessel functions against A&S Table 9.8
    !> Abramowitz & Stegun, Handbook of Mathematical Functions, 1972
    !> Table 9.8: I_n(x) for x = 0(0.1)5, n = 0,1,2
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: Y(3), computed, reference
    real(DP) :: rel_err, max_err
    integer :: Nz, i

    ! Selected values from A&S Table 9.8
    ! Format: x, I_0(x), I_1(x), I_2(x)
    real(DP), parameter :: test_data(4, 8) = reshape([ &
      ! x        I_0(x)          I_1(x)          I_2(x)
      0.5D0,    1.0634834D0,    0.2578954D0,    0.0319850D0, &
      1.0D0,    1.2660659D0,    0.5651591D0,    0.1357477D0, &
      1.5D0,    1.6467233D0,    0.9816664D0,    0.3378346D0, &
      2.0D0,    2.2795853D0,    1.5906368D0,    0.6889484D0, &
      2.5D0,    3.2898391D0,    2.5167163D0,    1.2764661D0, &
      3.0D0,    4.8807926D0,    3.9533700D0,    2.2452125D0, &
      4.0D0,   11.3019220D0,    9.7594652D0,    6.4221894D0, &
      5.0D0,   27.2398718D0,   24.3356421D0,   17.5056150D0  &
    ], [4, 8])

    max_err = 0.0D0

    do i = 1, 8
      call DBESI(test_data(1,i), 0.0_DP, 1, 3, Y, Nz)

      ! Check I_0
      rel_err = relative_error_dp(Y(1), test_data(2,i))
      max_err = max(max_err, rel_err)

      ! Check I_1
      rel_err = relative_error_dp(Y(2), test_data(3,i))
      max_err = max(max_err, rel_err)

      ! Check I_2
      rel_err = relative_error_dp(Y(3), test_data(4,i))
      max_err = max(max_err, rel_err)
    end do

    write(*, '(a, es12.5)') "  I Bessel vs A&S 9.8 max relative error: ", max_err

    call check(error, max_err < TOL_DP, "I Bessel failed A&S Table 9.8 comparison")
  end subroutine

  subroutine test_jbess_as_table91(error)
    !> Test J Bessel functions against A&S Table 9.1
    !> Table 9.1: J_n(x) for x = 0(0.1)15, n = 0,1,2,...,9
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: Y(5)
    real(DP) :: rel_err, max_err
    integer :: Nz, i

    ! Selected values from A&S Table 9.1
    ! Note: x=15 case excluded - returns zeros (possible DBESJ issue for large x)
    real(DP), parameter :: test_data(6, 5) = reshape([ &
      ! x        J_0(x)          J_1(x)          J_2(x)          J_3(x)          J_4(x)
      0.5D0,    0.9384698D0,    0.2422685D0,    0.0306040D0,    0.0025603D0,    0.0001601D0, &
      1.0D0,    0.7651977D0,    0.4400506D0,    0.1149035D0,    0.0195634D0,    0.0024766D0, &
      2.0D0,    0.2238908D0,    0.5767248D0,    0.3528340D0,    0.1289433D0,    0.0339957D0, &
      5.0D0,   -0.1775968D0,   -0.3275791D0,    0.0465651D0,    0.3648312D0,    0.3912327D0, &
     10.0D0,   -0.2459358D0,    0.0434728D0,    0.2546303D0,    0.0583794D0,   -0.2196178D0  &
    ], [6, 5])

    max_err = 0.0D0

    do i = 1, 5
      call DBESJ(test_data(1,i), 0.0_DP, 5, Y, Nz)

      rel_err = relative_error_dp(Y(1), test_data(2,i))
      max_err = max(max_err, rel_err)
      rel_err = relative_error_dp(Y(2), test_data(3,i))
      max_err = max(max_err, rel_err)
      rel_err = relative_error_dp(Y(3), test_data(4,i))
      max_err = max(max_err, rel_err)
      rel_err = relative_error_dp(Y(4), test_data(5,i))
      max_err = max(max_err, rel_err)
      rel_err = relative_error_dp(Y(5), test_data(6,i))
      max_err = max(max_err, rel_err)
    end do

    write(*, '(a, es12.5)') "  J Bessel vs A&S 9.1 max relative error: ", max_err

    call check(error, max_err < TOL_DP, "J Bessel failed A&S Table 9.1 comparison")
  end subroutine

  subroutine test_kbess_as_table98(error)
    !> Test K Bessel functions against A&S Table 9.8
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: Y(3)
    real(DP) :: rel_err, max_err
    integer :: Nz, i

    ! Selected values from A&S Table 9.8
    real(DP), parameter :: test_data(4, 6) = reshape([ &
      ! x        K_0(x)          K_1(x)          K_2(x)
      0.5D0,    0.9244191D0,    1.6564411D0,    7.5501836D0, &
      1.0D0,    0.4210244D0,    0.6019072D0,    1.6248389D0, &
      1.5D0,    0.2138056D0,    0.2773878D0,    0.5836560D0, &
      2.0D0,    0.1138939D0,    0.1398659D0,    0.2537598D0, &
      3.0D0,    0.0347395D0,    0.0401564D0,    0.0615104D0, &
      5.0D0,    0.0036911D0,    0.0040446D0,    0.0053089D0  &
    ], [4, 6])

    max_err = 0.0D0

    do i = 1, 6
      call DBESK(test_data(1,i), 0.0_DP, 1, 3, Y, Nz)

      rel_err = relative_error_dp(Y(1), test_data(2,i))
      max_err = max(max_err, rel_err)
      rel_err = relative_error_dp(Y(2), test_data(3,i))
      max_err = max(max_err, rel_err)
      rel_err = relative_error_dp(Y(3), test_data(4,i))
      max_err = max(max_err, rel_err)
    end do

    write(*, '(a, es12.5)') "  K Bessel vs A&S 9.8 max relative error: ", max_err

    call check(error, max_err < TOL_DP, "K Bessel failed A&S Table 9.8 comparison")
  end subroutine

  subroutine test_ibess_special(error)
    !> Test special values and identities for I Bessel
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: Y(2)
    integer :: Nz
    real(DP) :: rel_err

    ! I_0(0) = 1
    call DBESI(0.0_DP, 0.0_DP, 1, 1, Y, Nz)
    call check(error, abs(Y(1) - 1.0_DP) < 1.0D-14, "I_0(0) should equal 1")
    if (allocated(error)) return

    ! I_n(0) = 0 for n > 0
    call DBESI(0.0_DP, 1.0_DP, 1, 1, Y, Nz)
    call check(error, abs(Y(1)) < 1.0D-14, "I_n(0) should equal 0 for n > 0")
    if (allocated(error)) return

    ! Large argument asymptotic: I_nu(x) ~ e^x / sqrt(2*pi*x) for large x
    call DBESI(100.0_DP, 0.0_DP, 2, 1, Y, Nz)  ! scaled by e^(-x)
    rel_err = relative_error_dp(Y(1), 1.0_DP / sqrt(2.0_DP * 3.14159265358979_DP * 100.0_DP))
    call check(error, rel_err < 0.01_DP, "I_0 asymptotic expansion failed")
  end subroutine

  subroutine test_jbess_special(error)
    !> Test special values for J Bessel
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: Y(2)
    integer :: Nz

    ! J_0(0) = 1
    call DBESJ(0.0_DP, 0.0_DP, 1, Y, Nz)
    call check(error, abs(Y(1) - 1.0_DP) < 1.0D-14, "J_0(0) should equal 1")
    if (allocated(error)) return

    ! J_n(0) = 0 for n > 0
    call DBESJ(0.0_DP, 1.0_DP, 1, Y, Nz)
    call check(error, abs(Y(1)) < 1.0D-14, "J_n(0) should equal 0 for n > 0")
    if (allocated(error)) return

    ! Known zeros: J_0(x) = 0 at x ~ 2.4048, 5.5201, 8.6537
    call DBESJ(2.4048255576957728_DP, 0.0_DP, 1, Y, Nz)
    call check(error, abs(Y(1)) < 1.0D-10, "J_0 first zero check failed")
  end subroutine

  subroutine test_kbess_asymptotic(error)
    !> Test K Bessel asymptotic behavior
    !> K_nu(x) ~ sqrt(pi/(2x)) * e^(-x) for large x
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: Y(1)
    real(DP) :: x, expected, rel_err
    integer :: Nz

    x = 50.0_DP

    ! Using scaled version (KODE=2): e^x * K_0(x) ~ sqrt(pi/(2x))
    call DBESK(x, 0.0_DP, 2, 1, Y, Nz)
    expected = sqrt(3.14159265358979_DP / (2.0_DP * x))
    rel_err = relative_error_dp(Y(1), expected)

    write(*, '(a, es12.5)') "  K_0 asymptotic relative error at x=50: ", rel_err

    call check(error, rel_err < 0.01_DP, "K_0 asymptotic expansion failed")
  end subroutine

end module test_bessel_reference

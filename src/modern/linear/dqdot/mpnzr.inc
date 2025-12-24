!** MPNZR
SUBROUTINE MPNZR(Rs,Re,Z,Trunc)
  !> Subsidiary to DQDOTA and DQDOTI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (MPNZR-A)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !  Modified for use with BLAS.  Blank COMMON changed to named COMMON.
  !  Assumes long (i.e. (t+4)-DIGIT) fraction in R, sign = RS, exponent
  !  = RE.  Normalizes, and returns 'mp' result in Z. Integer arguments
  !  RS and RE are not preserved. R*-rounding is used if TRUNC=0
  !
  !  The argument Z(*) and the variable R in COMMON are INTEGER arrays
  !  of size 30.  See the comments in the routine MPBLAS for the reason
  !  for this choice.
  !
  !***
  ! **See also:**  DQDOTA, DQDOTI, MPBLAS
  !***
  ! **Routines called:**  MPERR, MPOVFL, MPUNFL
  !***
  ! COMMON BLOCKS    MPCOM

  !* REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
  USE MPCOM, ONLY : b_com, lun_com, m_com, t_com, r_com, mxr_com

  INTEGER, INTENT(IN) :: Rs, Trunc
  INTEGER, INTENT(INOUT) :: Re, Z(mxr_com)
  INTEGER :: i, i2, i2m, i2p, is, it, j, k, b2
  LOGICAL :: found_nonzero, do_round
  !* FIRST EXECUTABLE STATEMENT  MPNZR
  !
  ! Modernisation notes (202412):
  !   Eliminated GOTOs using logical flags.
  !   Original F77 used labels 100/150/200.
  !   Algorithm unchanged - Brent multiple-precision arithmetic.
  !
  i2 = t_com + 4
  found_nonzero = .FALSE.
  is = 0
  !
  IF( Rs/=0 ) THEN
    ! CHECK THAT SIGN = +-1
    IF( ABS(Rs)<=1 ) THEN
      ! LOOK FOR FIRST NONZERO DIGIT
      DO i = 1, i2
        is = i - 1
        IF( r_com(i)>0 ) THEN
          found_nonzero = .TRUE.
          EXIT
        END IF
        ! FRACTION ZERO
      END DO
    ELSE
      ERROR STOP ' *** SIGN NOT 0, +1 OR -1 IN CALL TO MPNZR, POSSIBLE OVERWRITING PROBLEM ***'
    END IF
  END IF
  !
  IF( .NOT. found_nonzero ) THEN
    ! STORE ZERO IN Z
    Z(1) = 0
    RETURN
  END IF
  !
  IF( is/=0 ) THEN
    ! NORMALIZE
    Re = Re - is
    i2m = i2 - is
    DO j = 1, i2m
      k = j + is
      r_com(j) = r_com(k)
    END DO
    i2p = i2m + 1
    DO j = i2p, i2
      r_com(j) = 0
    END DO
  END IF
  !
  ! CHECK TO SEE IF TRUNCATION IS DESIRED
  IF( Trunc==0 ) THEN
    ! SEE IF ROUNDING NECESSARY
    ! TREAT EVEN AND ODD BASES DIFFERENTLY
    b2 = b_com/2
    do_round = .FALSE.
    IF( (2*b2)/=b_com ) THEN
      ! ODD BASE, ROUND IF R(T+1)... > 1/2
      DO i = 1, 4
        it = t_com + i
        IF( r_com(it)<b2 ) EXIT
        IF( r_com(it)/=b2 ) THEN
          do_round = .TRUE.
          EXIT
        END IF
      END DO
    ELSE
      ! B EVEN.  ROUND IF R(T+1)>=B2 UNLESS R(T) ODD AND ALL ZEROS
      ! AFTER R(T+2).
      IF( r_com(t_com+1)>=b2 ) THEN
        do_round = .TRUE.
        IF( r_com(t_com+1)==b2 ) THEN
          IF( MOD(r_com(t_com),2)/=0 ) THEN
            IF( (r_com(t_com+2)+r_com(t_com+3)+r_com(t_com+4))==0 ) do_round = .FALSE.
          END IF
        END IF
      END IF
    END IF
    !
    IF( do_round ) THEN
      ! ROUND
      DO j = 1, t_com
        i = t_com + 1 - j
        r_com(i) = r_com(i) + 1
        IF( r_com(i)<b_com ) EXIT
        r_com(i) = 0
        IF( j==t_com ) THEN
          ! EXCEPTIONAL CASE, ROUNDED UP TO .10000...
          Re = Re + 1
          r_com(1) = 1
        END IF
      END DO
    END IF
  END IF
  !
  ! CHECK FOR OVERFLOW
  IF( Re>m_com ) THEN
    WRITE (lun_com,99002)
    99002 FORMAT (' *** OVERFLOW OCCURRED IN MPNZR ***')
    CALL MPOVFL(Z)
    RETURN
    ! CHECK FOR UNDERFLOW
  ELSEIF( Re<(-m_com) ) THEN
    ! UNDERFLOW HERE
    CALL MPUNFL(Z)
    RETURN
  END IF
  !
  ! STORE RESULT IN Z
  Z(1) = Rs
  Z(2) = Re
  DO i = 1, t_com
    Z(i+2) = r_com(i)
  END DO
  !
END SUBROUTINE MPNZR
!** MPDIVI
SUBROUTINE MPDIVI(X,Iy,Z)
  !> Subsidiary to DQDOTA and DQDOTI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (MPDIVI-A)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !  Divides 'mp' X by the single-precision integer IY giving 'mp' Z.
  !  This is much faster than division by an 'mp' number.
  !
  !  The arguments X(*) and Z(*), and the variable R in COMMON are all
  !  INTEGER arrays of size 30.  See the comments in the routine MPBLAS
  !  for the reason for this choice.
  !
  !***
  ! **See also:**  DQDOTA, DQDOTI, MPBLAS
  !***
  ! **Routines called:**  MPCHK, MPERR, MPNZR, MPSTR, MPUNFL
  !***
  ! COMMON BLOCKS    MPCOM

  !* REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   ??????  Modified for use with BLAS.  Blank COMMON changed to named
  !           COMMON.  R given dimension 12.
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
  !   202412  Eliminated GOTOs using DO loops and logical flags.  (modernization)
  USE MPCOM, ONLY : b_com, lun_com, m_com, t_com, r_com, mxr_com

  INTEGER, INTENT(IN) :: Iy
  INTEGER, INTENT(IN) :: X(mxr_com)
  INTEGER, INTENT(INOUT) :: Z(mxr_com)
  INTEGER :: i, i2, iq, iqj, ir, j, j1, j11, j2, k, kh, rs, re, r1, c, c2, b2
  LOGICAL :: overflow_error, found_digit, do_large_case
  !* FIRST EXECUTABLE STATEMENT  MPDIVI
  ! Modernized: replaced GOTOs with DO WHILE loops and flags
  overflow_error = .FALSE.
  found_digit = .FALSE.
  do_large_case = .FALSE.

  rs = X(1)
  j = Iy
  IF( j<0 ) THEN
    j = -j
    rs = -rs
  ELSEIF( j==0 ) THEN
    WRITE (lun_com,99001)
    99001 FORMAT (' *** ATTEMPTED DIVISION BY ZERO IN CALL TO MPDIVI ***')
    ERROR STOP ' *** INTEGER OVERFLOW IN MPDIVI, b_com TOO LARGE ***'
    Z(1) = 0
    RETURN
  END IF
  re = X(2)
  ! CHECK FOR ZERO DIVIDEND
  IF( rs/=0 ) THEN
    ! CHECK FOR DIVISION BY B
    IF( j==b_com ) THEN
      CALL MPSTR(X,Z)
      IF( re<=(-m_com) ) THEN
        ! UNDERFLOW HERE
        CALL MPUNFL(Z)
        RETURN
      ELSE
        Z(1) = rs
        Z(2) = re - 1
        RETURN
      END IF
      ! CHECK FOR DIVISION BY 1 OR -1
    ELSEIF( j/=1 ) THEN
      c = 0
      i2 = t_com + 4
      i = 0
      ! IF J*B NOT REPRESENTABLE AS AN INTEGER HAVE TO SIMULATE
      ! LONG DIVISION.   ASSUME AT LEAST 16-BIT WORD.
      b2 = MAX(8*b_com,32767/b_com)
      IF( j>=b2 ) THEN
        ! HERE NEED SIMULATED DOUBLE-PRECISION DIVISION
        c2 = 0
        j1 = j/b_com
        j2 = j - j1*b_com
        j11 = j1 + 1
        find_digit_large: DO
          ! LOOK FOR FIRST NONZERO DIGIT
          i = i + 1
          c = b_com*c + c2
          c2 = 0
          IF( i<=t_com ) c2 = X(i+2)
          IF( c<j1 ) THEN
            CYCLE find_digit_large
          ELSEIF( c==j1 ) THEN
            IF( c2>=j2 ) THEN
              found_digit = .TRUE.
              do_large_case = .TRUE.
              EXIT find_digit_large
            END IF
          ELSE
            found_digit = .TRUE.
            do_large_case = .TRUE.
            EXIT find_digit_large
          END IF
        END DO find_digit_large
      ELSE
        find_digit_small: DO
          ! LOOK FOR FIRST NONZERO DIGIT IN QUOTIENT
          i = i + 1
          c = b_com*c
          IF( i<=t_com ) c = c + X(i+2)
          r1 = c/j
          IF( r1<0 ) THEN
            overflow_error = .TRUE.
            EXIT find_digit_small
          END IF
          IF( r1/=0 ) THEN
            ! ADJUST EXPONENT AND GET T+4 DIGITS IN QUOTIENT
            re = re + 1 - i
            r_com(1) = r1
            c = b_com*(c-j*r1)
            kh = 2
            IF( i<t_com ) THEN
              kh = 1 + t_com - i
              DO k = 2, kh
                i = i + 1
                c = c + X(i+2)
                r_com(k) = c/j
                c = b_com*(c-j*r_com(k))
              END DO
              IF( c<0 ) THEN
                overflow_error = .TRUE.
                EXIT find_digit_small
              END IF
              kh = kh + 1
            END IF
            DO k = kh, i2
              r_com(k) = c/j
              c = b_com*(c-j*r_com(k))
            END DO
            IF( c<0 ) overflow_error = .TRUE.
            EXIT find_digit_small
          END IF
        END DO find_digit_small
      END IF
    ELSE
      CALL MPSTR(X,Z)
      Z(1) = rs
      RETURN
    END IF
  END IF

  IF( overflow_error ) THEN
    ERROR STOP ' *** INTEGER OVERFLOW IN MPDIVI, b_com TOO LARGE ***'
    Z(1) = 0
    RETURN
  END IF

  ! Handle large case (was label 200-300 loop)
  IF( do_large_case ) THEN
    re = re + 1 - i
    k = 1
    ! Main loop for large divisor case
    large_loop: DO
      ! GET APPROXIMATE QUOTIENT FIRST
      ir = c/j11
      ! NOW REDUCE SO OVERFLOW DOES NOT OCCUR
      iq = c - ir*j1
      IF( iq>=b2 ) THEN
        ! HERE IQ*B WOULD POSSIBLY OVERFLOW SO INCREASE IR
        ir = ir + 1
        iq = iq - j1
      END IF
      iq = iq*b_com - ir*j2
      IF( iq<0 ) THEN
        ! HERE IQ NEGATIVE SO IR WAS TOO LARGE
        ir = ir - 1
        iq = iq + j
      END IF
      IF( i<=t_com ) iq = iq + X(i+2)
      iqj = iq/j
      ! R(K) = QUOTIENT, C = REMAINDER
      r_com(k) = iqj + ir
      c = iq - j*iqj
      IF( c<0 ) THEN
        ! CARRY NEGATIVE SO OVERFLOW MUST HAVE OCCURRED
        ERROR STOP ' *** INTEGER OVERFLOW IN MPDIVI, b_com TOO LARGE ***'
        Z(1) = 0
        RETURN
      END IF
      ! MAIN LOOP FOR LARGE ABS(IY) CASE
      k = k + 1
      IF( k>i2 ) EXIT large_loop
      i = i + 1
    END DO large_loop
  END IF

  ! NORMALIZE AND ROUND RESULT
  CALL MPNZR(rs,re,Z,0)
  RETURN
END SUBROUTINE MPDIVI